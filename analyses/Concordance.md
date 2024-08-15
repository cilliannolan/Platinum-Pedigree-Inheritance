# Pedigree filtering and phasing using inheritance vectors

 The concordance program accepts a VCF (as detailed in variant calling) and the transmission vectors (as detailed in defining haplotype structure). The concordance code generates all possible combinations of valid genotypes over a haplotype (region) and tests if the observed genotypes at a site match a valid combination. Sites that pass the inheritance filtering are kept if they contain a genotype call for every individual (excluding any site with nocalls). Using the inheritance vectors we can directly phase genotypes for bi-allelic and multi-allelic sites. There is a bijective mapping of observed genotypes in the VCF with the haplotype structure with the exception of all heterozygous genotypes, where the phase is ambiguous. We exclude all heterozygous sites from the truth set. The command used to filter and phase the variants is:

```
concordance \
--father NA12877 \
--mother NA12878 \
--vcf {}.vcf.gz \
--inheritance data/ceph.GRCh38.viterbi.csv \
--prefix {} > std.out
```
