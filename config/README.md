# Configuration options

A description of every valid option in `config.yaml`.

```
# If you're performing SBS96 signature analysis, sample name should have _SBS96 suffix.
# If you're performing DBS78 signature analysis, sample name should have _DBS78 suffix.
# If you're performing ID83 signature analysis, sample name should have _ID83 suffix.
# I know. This is not ideal. But it works for the time being. Note again this is not for VCFs based on haploid genome alignment.

samples:
  <your_sample_name_prefix>_SBS96:
    vcf: "<path.to.your.snv.vcf[.gz]>"
    signatures: ["SBS96"]

  <your_sample_name_prefix>_ID83:
    vcf: "path.to.your.indels.vcf[.gz]>"
    signatures: ["ID83"]

reference_genome: "path/to/your/reference/genome.fasta"

output_dir: "results"
```
