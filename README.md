# HERMES2 GWAS

This repository contains code to reproduce key analyses in GWAS of heart failure subtypes and dilated cardiomyopathy from the HERMES Consortium.

### :warning: This repository is currently under active development. Contents will be frequently updated.

The repository is structured as follows (details below):

```
.
├── data
├── results
├── scripts
└── workflow
```

## `data`

The `data` directory contains data required to perform some the analyses, including study-level data for QC and meta-analysis and
software-specific datasets (e.g. reference genotype to construct LD matrix).

Due to size and sharing policy, the actual data used in the analysis are *not* stored in this repository but some can be obtained from public repositories.

## `results`

The `results` directory 

## `scripts`

The `scripts` directory contains custom software / scripts required to perform some analyses.
Some of these can be directly downloaded from public repositories e.g. [`ldsc`](https://github.com/bulik/ldsc) and [`polyfun`](https://github.com/omerwe/polyfun)
and can be symbolic links (e.g. from existing copies)

## `workflow`

The `workflow` directory contains a [`snakemake`](https://snakemake.readthedocs.io/en/stable/index.html) pipeline to perform most of the analyses.
Each key analysis is represented by a set of rules, as listed on `workflow/rules/` directory.
For reproduciblity purposes, most rules can be run in an isolated `conda` environment or `singularity` container using `--use-conda` and `--use-singularity` arguments, respectively.
Please consult the [`snakemake` documentation](https://snakemake.readthedocs.io/en/stable/index.html) for guidance and additional options.




