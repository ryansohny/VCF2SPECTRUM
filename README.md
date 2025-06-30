[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/CI/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
# Make Mutation Spectrum from VCF file (SNV/Indels)
Currently support SBS96, DBS78 and ID83 signatures.
> Note:
> SBS96 is based on VCF obtained from diploid genome alignment.
> If your VCF is generated using haploid genome (e.g., hg38/CHM13-T2T) alignment, this will undercount homozygous alt alleles.

This is a Snakemake project template. The `Snakefile` is under `workflow`.

[Slides](https://mrvollger.github.io/SmkTemplate/slides) describing and justifying the use of this template.

## Install

Please start by installing [pixi](https://pixi.sh/latest/) which handles the environment of this Snakemake workflow.

You can then install the `pixi` environment by cloning this repository and running:

```bash
pixi install
```

## Usage

`pixi` handles the execution of the Snakemake workflows:

```bash
pixi run snakemake --cores 4 --use-conda --conda-prefix /mmfs1/gscratch/stergachislab/mhsohny/Miniconda/envs/sigprofiler
```
