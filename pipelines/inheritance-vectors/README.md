# Inheritance vectors

## Workflow steps

1. Prepare SNPs (`prepare_snps.py`)

This step assigns haplotype identity for each SNP in the cohort vcf, by taking the identity from grandparents (paternal: A, B, maternal: C, D).

It is run as follows:

```bash
python {input.script} \
    --cohort-calls {input.cohort_vcf} \
    --pedigree {input.cohort_ped} \
    --dad-sample {params.dad} \
    --mom-sample {params.mom} \
    --subset-children {params.children} \
    --filtered {output.filtered} \
    > {output.bed}
```

It ouputs a TSV file with the following structure:

```
| #CHROM | start | end  | REF | ALT | called_parent | grandparent | phase | children_calls |
|--------|-------|------|-----|-----|---------------|-------------|-------|----------------|
| chr1   | 4736  | 4736 | C   | A   | NA12886       | NA12878     | B     | 200102         |
| chr1   | 5324  | 5324 | C   | A   | NA12886       | NA12878     | B     | 200102         |
| chr1   | 6137  | 6137 | G   | A   | NA12886       | NA12878     | B     | 200102         |
```

2. Viterbi
3. Output

## Run the workflow

Set up and activate your snakemake environment using conda to run the workflow:

```bash
conda create -n snakemake -f envs/snakemake.yaml
conda activate snakemake
```

Run a test dataset using just chromosome 22:

```bash
VCF="data/chr22.dv.all.vcf.gz"

snakemake \
  -s inhHMM.smk \
  --config \
  "input_vcf=${VCF}" \
  "mom=NA12878" \
  "dad=NA12877" \
  "children=NA12879,NA12881,NA12882,NA12883,NA12884,NA12885,NA12886,NA12887" \
  "male_children=NA12882,NA12883,NA12884,NA12886" \
  "chromosomes=chr22" \
  --use-conda \
  -p
```

Final inheritance vectors will be output to `output/viterbi/{vcf_prefix}.inht_vectors.csv`

To run on the entire genome just remove the `chromosome` parameter from the above command, and update your vcf path, ie:

```bash
VCF="data/cohort.vcf.gz"

snakemake \
  -s inhHMM.smk \
  --config \
  "input_vcf=${VCF}" \
  "mom=NA12878" \
  "dad=NA12877" \
  "children=NA12879,NA12881,NA12882,NA12883,NA12884,NA12885,NA12886,NA12887" \
  "male_children=NA12882,NA12883,NA12884,NA12886" \
  --use-conda \
  -p
```