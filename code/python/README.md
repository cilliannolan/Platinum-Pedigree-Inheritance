# SV merging

## VCF normalization

Prior to SV merging it is important to handle the differences of SV callers by normalizing the VCF files. Most notably this includes dropping INFO/FORMAT fields that are not shared across different callers, reformatting fields, and prefixing SV caller id to variant IDs. This can be done using `strip_vcf.py`, given an input json file containing VCF file paths (*.pedfilt.vcf files), their caller identity and (optional) prefix to add to each variant ID (see `strip_vcf.json` for an example):

```
{
"caller": "sawfish",
"path": "/.../sawfish.pedfilt.vcf",
"prefix_tag": "sawfish"
}
```
and an output path. Which will write a normalized VCF file (*.strip.pedfilt.vcf) for each input. Supported SV callers: sawfish, sniffles, pav, pggb, pbsv.

```
python strip_vcf.py -i strip_vcf.json -o ../sv_vcfs/out/
```

## VCF merging

After normalization the VCF files can be merged using `interval_tree_merge.py` given the space separated paths of normalized VCFs which also determines the merging ordering and a output file. Optional CLI `--flank_len` (default = 200), allows for setting the size of flanks when building intervals.

The below command (assuming each SV VCF file is present), first merge pbsv into sawfish, then sniffles in that result, and finally do the same with PAV.

```
python interval_tree_merge.py --vcf_files ../sv_vcfs/out/sawfish.vcf ../sv_vcfs/out/pbsv.vcf ../sv_vcfs/out/sniffles.vcf ../sv_vcfs/out/pav.vcf -o merged.vcf 
```