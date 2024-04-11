use std::collections::{HashMap, HashSet};
use std::io::Read as IoRead;
use std::mem::{self, swap};
use std::path::Path;

use clap::Parser;
use csv::ReaderBuilder;
use edit_distance::edit_distance;
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::record::{Genotype, GenotypeAllele};
use rust_htslib::bcf::{Format, Header, Read, Reader, Writer};

/// Filter variants based on expected haplotypes
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// CSV containing the inheritance vectors
    #[arg(short, long)]
    inheritance: String,

    /// VCF file
    #[arg(short, long)]
    vcf: String,

    #[arg(short, long)]
    prefix: String,

    #[arg(short, long)]
    father: String,

    #[arg(short, long)]
    mother: String,
}

const POSSIBLE_PATTERNS: [[i32; 4]; 15] = [
    //[ 0, 0, 0, 0], // dad homref mom homref
    [0, 0, 0, 1], // dad homref mom het (1)
    [0, 0, 1, 0], // dad homref mom het (2)
    [0, 0, 1, 1], // dad homref mom homalt
    [0, 1, 0, 0], // dad het mom homref (1)
    [1, 0, 0, 0], // dad het mom homref (2)
    [1, 0, 1, 0], // dad het mom het (1)
    [0, 1, 0, 1], // dad het mom het (2)
    [1, 0, 0, 1], // dad het mom het (3)
    [0, 1, 1, 0], // dad het mom het (4)
    [0, 1, 1, 1], // dad het mom homalt (1)
    [1, 0, 1, 1], // dad het mom homalt (2)
    [1, 1, 0, 0], // dad is homalt mom is homref
    [1, 1, 0, 1], // dad is homalt mom is het (1)
    [1, 1, 1, 0], // dad is homalt mom is het (2)
    [1, 1, 1, 1], // dad is homalt mom is homalt
];

fn geno_conversion(geno: String) -> i32 {
    match geno.as_str() {
        "0/0" => 0,
        "0|0" => 0,
        "0/1" => 1,
        "1/0" => 1,
        "0|1" => 1,
        "1|0" => 1,
        "1/1" => 2,
        "1|1" => 2,
        "./." => 3,
        ".|." => 3,
        _ => 3,
    }
}

fn allele_conversion(allele: char) -> usize {
    match allele {
        'A' => 0,
        'B' => 1,
        'C' => 2,
        'D' => 3,
        _ => 4,
    }
}

struct InheritanceBlock {
    chrom: String,
    start: i32,
    end: i32,
    passing_count: i32,
    failing_count: i32,
    samples: Vec<String>,
    sample_lookups: HashMap<String, usize>,
    parental_hap: Vec<String>,
    patterns: HashMap<String, Vec<[i32; 2]>>,
    inherited_haps: HashSet<char>,
}
// returns a HashMap were the key is a simplified genotype string (0=homref,1=het,2=homalt), and the corresponding phase
fn okay_genotypes(iht: &InheritanceBlock) -> HashMap<String, Vec<[i32; 2]>> {
    let mut lookup: HashMap<String, Vec<[i32; 2]>> = HashMap::new();
    // possible genotypes
    for i in POSSIBLE_PATTERNS {
        let mut onevec: String = "".to_string();
        let mut phased_genotypes: Vec<[i32; 2]> = Vec::new();

        // possible haplotype configurations
        for j in &iht.parental_hap {
            let possible_geno = format!(
                "{}|{}",
                i[allele_conversion(j.chars().nth(0).unwrap())],
                i[allele_conversion(j.chars().nth(1).unwrap())]
            );
            phased_genotypes.push([
                i[allele_conversion(j.chars().nth(0).unwrap())],
                i[allele_conversion(j.chars().nth(1).unwrap())],
            ]);

            onevec = format!("{}{}", onevec, geno_conversion(possible_geno))
        }

        if lookup.contains_key(&onevec) {
            let seen_phase = lookup.get(&onevec).unwrap();
            for (g1, g2) in seen_phase.iter().zip(phased_genotypes.iter()) {
                if g1 != g2 {
                    println!(
                        "WARNING disagreement in phase: {} {:?} {:?}",
                        onevec, seen_phase, phased_genotypes
                    );
                }
            }
        } else {
            lookup.insert(onevec, phased_genotypes);
        }
    }
    return lookup;
}

impl std::fmt::Display for InheritanceBlock {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "inheritance-block: {} {} {} {:?} {:?} {:?}",
            self.chrom, self.start, self.end, self.samples, self.parental_hap, self.patterns
        )
    }
}

fn parse_inht(inht_fn: String) -> Vec<InheritanceBlock> {
    let mut inht_info = Vec::new();
    let mut reader = ReaderBuilder::new().from_path(inht_fn);

    let header = reader
        .as_mut()
        .expect("Error reading inheritance CSV header")
        .headers()
        .expect("Error reading inheritance CSV header")
        .clone();

    for record in reader.expect("Error reading inheritance CSV.").records() {
        let mut ihtblock = InheritanceBlock {
            chrom: record.as_ref().unwrap()[0].to_string().clone(),
            start: record.as_ref().unwrap()[1].parse::<i32>().unwrap().clone(),
            end: record.as_ref().unwrap()[2].parse::<i32>().unwrap().clone(),
            passing_count: 0,
            failing_count: 0,
            samples: Vec::new(),
            sample_lookups: HashMap::new(),
            parental_hap: Vec::new(),
            patterns: HashMap::new(),
            inherited_haps: HashSet::new(),
        };

        let mut one_parent = false;

        let mut sidx: usize = 0;
        for i in 3..header.len() {
            // Some inheritance vectors do not contain both parental marker, for example you might see `A` rather than `AB`.
            // We want to skip these sites as the expected genotypes are unknown.
            if record.as_ref().unwrap()[i].to_string().len() == 1 {
                one_parent = true;
            }
            // putting the sample names in the header into the inheritance block for easy lookup
            ihtblock
                .sample_lookups
                .insert((&header[i]).to_string(), sidx);
            ihtblock.samples.push((&header[i]).to_string());
            sidx += 1;
            // geno
            ihtblock
                .parental_hap
                .push(record.as_ref().unwrap()[i].to_string());
        }
        // counting up haplotypes seen in children
        for i in 0..ihtblock.parental_hap.len() - 2 {
            let geno = ihtblock.parental_hap.get(i).unwrap();
            ihtblock.inherited_haps.insert(geno.as_bytes()[0].into());
            ihtblock.inherited_haps.insert(geno.as_bytes()[1].into());
        }

        if one_parent {
            println!("Warning skipping block missing both parents {}", ihtblock);
            continue;
        }
        ihtblock.patterns = okay_genotypes(&ihtblock);
        println!("{}", ihtblock);
        inht_info.push(ihtblock);
    }
    return inht_info;
}

fn header_to_idx(h: &HeaderView) -> HashMap<String, usize> {
    let mut lookup = HashMap::new();
    for (i, mut x) in (*h).samples().into_iter().enumerate() {
        let mut s = String::new();
        x.read_to_string(&mut s).expect("converting name to idx"); // Read sample name in to `s`
        println!("idx:{} sample:{}", i, s); // output sample name
        lookup.insert(s, i);
    }
    return lookup;
}

fn get_iht_block<'a>(
    ihts: &'a Vec<InheritanceBlock>,
    chr: &str,
    pos: i32,
    current: &mut usize,
) -> Option<&'a InheritanceBlock> {
    for i in *current..ihts.len() {
        if (ihts[i].chrom == chr) && (pos >= ihts[i].start) && (pos <= ihts[i].end) {
            return Some(&ihts[i]);
        }
        *current += 1;
    }
    return None;
}

fn concordant(
    parents: [[GenotypeAllele; 4]; 4],
    haps: Vec<String>,
    genos: Vec<GenotypeAllele>,
) -> i8 {
    for (i, c) in parents.iter().enumerate() {
        let mut genovec: Vec<GenotypeAllele> = Vec::new();
        for p in &haps {
            let mut first_allele = c[allele_conversion(p.chars().nth(0).unwrap())];
            let mut second_allele = c[allele_conversion(p.chars().nth(1).unwrap())];

            if first_allele.index() > second_allele.index() {
                swap(&mut first_allele, &mut second_allele);
            }
            genovec.push(first_allele);
            genovec.push(second_allele);
        }
        if genovec == genos {
            return i as i8;
        }
    }
    return -1;
}

fn build_phased_haplotypes(
    parents: [GenotypeAllele; 4],
    haps: Vec<String>,
    sample_names: Vec<String>,
) -> HashMap<String, [GenotypeAllele; 2]> {
    let mut sample_to_alleles: HashMap<String, [GenotypeAllele; 2]> = HashMap::new();

    for (i, p) in haps.iter().enumerate() {
        sample_to_alleles.insert(
            sample_names[i].clone(),
            [
                parents[allele_conversion(p.chars().nth(0).unwrap())],
                parents[allele_conversion(p.chars().nth(1).unwrap())],
            ],
        );
    }

    return sample_to_alleles;
}

fn main() {
    let args = Args::parse();
    let mut bcf = Reader::from_path(args.vcf).expect("Error opening vcf file.");
    let header = bcf.header().clone();
    let wheader: Header = Header::from_template(&header);
    let mut inheritance = parse_inht(args.inheritance);
    let ped_idx_lookup = header_to_idx(&header);
    let mut outvcf = Writer::from_path(
        Path::new(&format!("{}.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();

    let mut outfailvcf = Writer::from_path(
        Path::new(&format!("{}.failedsites.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();

    let mut block_fail: i64 = 0;
    let mut failed: i64 = 0;
    let mut passed: i64 = 0;
    let mut all_het: i64 = 0;
    let mut fail_counts: HashMap<String, i64> = HashMap::new();

    for (_i, record_result) in bcf.records().enumerate() {
        let mut record = record_result.expect("Fail to read record");

        let gts = record.genotypes().expect("Error reading genotypes");
        let chr = std::str::from_utf8(header.rid2name(record.rid().unwrap()).expect("Invalid rid"))
            .unwrap();
        let mut current_block_idx: usize = 0;

        let block = get_iht_block(
            &inheritance,
            chr,
            record.pos().try_into().unwrap(),
            &mut current_block_idx,
        );
        match block {
            None => {
                println!(
                    "Warning: skipping variant as it is not in a block {} {}",
                    chr,
                    record.pos()
                );
                block_fail += 1;
                continue;
            }
            Some(_) => {}
        }
        let mut onevec: String = "".to_string();
        let mut genovec: Vec<GenotypeAllele> = Vec::new();

        let mut failed_site = false;
        let mut alt_count: usize = 0;
        let mut het_count: usize = 0;

        // samples are ordered by inheritance vector header
        for s in block.unwrap().samples.iter() {
            onevec = format!(
                "{}{}",
                onevec,
                geno_conversion(gts.get(ped_idx_lookup[s]).to_string())
            );
            let gt = gts.get(ped_idx_lookup[s]);
            let mut first_allele = gt[0];
            let mut second_allele = gt[1];

            if first_allele == GenotypeAllele::UnphasedMissing {
                failed_site = true;
                break;
            }
            // this will blow up unless the previous if statement is removed
            if first_allele.index().unwrap() > 0 || second_allele.index().unwrap() > 0 {
                alt_count += 1;
            }

            if first_allele.index().unwrap() != second_allele.index().unwrap() {
                het_count += 1;
            }

            // ensuring that the alleles are ordered allowing us to test if we can get the right configuration.
            if first_allele.index() > second_allele.index() {
                swap(&mut first_allele, &mut second_allele);
            }
            genovec.push(first_allele);
            genovec.push(second_allele);
        }

        // if the alt count is zero, we have a hom-ref site.
        if alt_count == 0 {
            failed_site = true;
        }

        let mother_gt = gts.get(ped_idx_lookup[&args.mother]);
        let father_gt = gts.get(ped_idx_lookup[&args.father]);

        let mut parent_allele_count: HashSet<GenotypeAllele> = HashSet::new();

        let configurations = [
            [father_gt[0], father_gt[1], mother_gt[0], mother_gt[1]],
            [father_gt[1], father_gt[0], mother_gt[0], mother_gt[1]],
            [father_gt[0], father_gt[1], mother_gt[1], mother_gt[0]],
            [father_gt[1], father_gt[0], mother_gt[1], mother_gt[0]],
        ];

        parent_allele_count.insert(father_gt[0]);
        parent_allele_count.insert(father_gt[1]);
        parent_allele_count.insert(mother_gt[0]);
        parent_allele_count.insert(mother_gt[1]);

        // if everything is het we check that the parent alleles are the same too.
        if het_count == block.unwrap().samples.len() && parent_allele_count.len() == 2 {
            all_het += 1;
            failed_site = true;
        }

        let con = concordant(configurations, block.unwrap().parental_hap.clone(), genovec);

        // unable to find a genotype configuration that matches the inheritance vector
        if con == -1 {
            failed_site = true;
        }

        if failed_site {
            failed += 1;
            inheritance[current_block_idx].failing_count += 1;
            *fail_counts.entry(onevec.clone()).or_insert(0) += 1;

            outfailvcf.write(&record).unwrap();

            continue;
        }

        let mut new_gts: Vec<GenotypeAllele> = Vec::new();
        let sample_count = usize::try_from(record.sample_count()).unwrap();

        let phased = build_phased_haplotypes(
            configurations[con as usize],
            block.unwrap().parental_hap.clone(),
            block.unwrap().samples.clone(),
        );

        for sample_index in 0..sample_count {
            let sample_name = String::from_utf8(header.samples()[sample_index].to_vec()).unwrap();
            let gt = gts.get(sample_index);

            if phased.contains_key(&sample_name) {
                let pg = phased.get(&sample_name).unwrap();

                let a0: u32 = pg[0].index().unwrap();
                let a1: u32 = pg[1].index().unwrap();

                new_gts.push(GenotypeAllele::Unphased(a0 as i32));
                new_gts.push(GenotypeAllele::Phased(a1 as i32));
            } else {
                for g in gt.iter() {
                    new_gts.push(g.clone());
                }
            }
        }

        record.push_genotypes(&new_gts).unwrap();

        passed += 1;
        outvcf.write(&record).unwrap();
        inheritance[current_block_idx].passing_count += 1;
    }

    for (ov, c) in &fail_counts {
        println!("failed_reason\t{}\t{}", ov, c);
    }

    for b in &inheritance {
        println!(
            "block_stats\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            args.prefix,
            b.chrom,
            b.start,
            b.end,
            b.passing_count,
            b.failing_count,
            b.passing_count + b.failing_count,
            b.inherited_haps.len(),
        );
    }

    println!(
        "not in block: {} passed: {} failed: {} all-het: {} ",
        block_fail, passed, failed, all_het
    );
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_geno_hom() {
        assert_eq!(geno_conversion("0/0".to_string()), 0);
    }
    #[test]
    fn test_geno_hom_phase() {
        assert_eq!(geno_conversion("0|0".to_string()), 0);
    }
    #[test]
    fn test_geno_het1() {
        assert_eq!(geno_conversion("0/1".to_string()), 1);
    }
    #[test]
    fn test_geno_het2() {
        assert_eq!(geno_conversion("1/0".to_string()), 1);
    }
    #[test]
    fn test_geno_het3() {
        assert_eq!(geno_conversion("1|0".to_string()), 1);
    }
    #[test]
    fn test_geno_het4() {
        assert_eq!(geno_conversion("0|1".to_string()), 1);
    }
    #[test]
    fn test_geno_nc1() {
        assert_eq!(geno_conversion("./1".to_string()), 3);
    }
    #[test]
    fn test_geno_nc2() {
        assert_eq!(geno_conversion("./.".to_string()), 3);
    }
    #[test]
    fn test_geno_nc3() {
        assert_eq!(geno_conversion("./0".to_string()), 3);
    }
}
