# HERMES 2 GWAS

This repository contains code and computational pipeline
to reproduce key analyses performed in the genome-wide association study (GWAS)
of heart failure subtypes and dilated cardiomyopathy manuscripts (2024)
from the HERMES Consortium.

The analyses are wrapped in a [*snakemake*](https://snakemake.github.io/) workflow management system
as illustrated below:

![Workflow rule graph](workflow/dag/test.svg)

The figure above is a toy example which performed the following:
* Two-step quality control (QC) of GWAS summary statistics of a phenotype labelled as `Pheno1_EUR` 
from study A and B following the procedure
described in [Winkler T, *et al.* 2014](https://www.nature.com/articles/nprot.2014.071) 
* Fixed-effect meta-analysis of GWAS of `Pheno1_EUR` using [METAL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922887/) 
* Conditional-joint association analysis using GWAS summary statistics from the meta-analysis with (GCTA-CoJo)[https://yanglab.westlake.edu.cn/software/gcta/#COJO]
* SNP Heritability assessment with [LDAK](https://dougspeed.com/snp-heritability/)  

Details of each step of the analysis, as represented by a *rule* (polygon in the figure),
are described on the [Rules](#rules) section

## Rules

| Rule                        | Description                                                                                                                                                                                                  |
| :-------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `qc_step1`                  | **Step 1 QC**: sanity check, create unique variant ID, harmonise allele (see [Note](#sec-qc-note) )                                                                                                          |
| `qc_step2`                  | **Step 2 QC**: QC based on allele comparison with reference panel (see [Note](#sec-qc-note) ). For I/O efficiency, this step also makes an *AFCHECK* plot (allele frequency comparison with reference panel) |
| `plot_qq`                   | Make *QQPLOT* (observed vs. expected log P-value)                                                                                                                                                            |
| `plot_pz`                   | Make *PZPLOT* (reported vs. calculated P-value)                                                                                                                                                              |
| `metal_make_config`         | Create config file for meta-analysis with METAL                                                                                                                                                              |
| `metal_run`                 | Run meta-analysis with METAL                                                                                                                                                                                 |
| `extract_N_meta`            | Extract N sample and other metadata from METAL log file                                                                                                                                                      |
| `format_metal`              | Reformat METAL output                                                                                                                                                                                        |
| `dl_ldak_tagging_file`      | Download pre-computed [LDAK tagging file](https://dougspeed.com/pre-computed-tagging-files/)                                                                                                                 |
| `ldak_format_sumstats_meta` | Reformat GWAS summary statistics for LDAK input                                                                                                                                                              |
| `ldak_sumher_meta`          | Perform LDAK SumHer SNP heritability estimation                                                                                                                                                              |
| `gcta_slct_aggregate`       | Perform chromosome-wide conditional analysis with GCTA-CoJo                                                                                                                                                  |
| `all`                       | Dummy rule to get target outputs (adjust in `workflow/Snakefile`)                                                                                                                                            |

**Note** {#sec-qc-note}

**Step 1 QC**

-   For sanity check, the `qc_step1` rule **excludes** variants with **any** of the following criteria:
    -   beta \> 10
    -   standard error \> 10
    -   *P* value outside 0-1 range
    -   imputation (INFO) score outside 0-1 range
    -   Allele frequencies outside 0-1 range
    -   N effective \< 50
    -   imputation (INFO) score \< 0.6
    -   If INFO score is missing & N effective cannot be calculated:
        -   Minor allele frequency (MAF) \< 0.01
-   N effective (effective sample size) is calculated as $N_{eff} = 2 \times MAF \times (1-MAF) \times N_{total} \times INFO$
-   Each variant will be assigned a unique ID in the format **chr:pos:A1_A2**, where **chr:pos** refers to chromosome and base pair position according to the [NCBI GRCh37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) genome assembly, and **A1**\_**A2** refers to allele 1 (effect allele) and allele 2 (other allele) in alphabetical order.
-   Accordingly, the regression coefficient (i.e. beta / log odds) of each variant is harmonised to reflect the effect allele (A1)

**Step 2 QC**

-   Based on allele comparison with reference panel, the `qc_step2` rule further **excludes** variants with **any** of the following criteria:
    -   unique variant ID not found in the reference panel (as the unique variant ID is constructed using genomic position and allele information, this will exclude any mismatch on those)
    -   MAF difference with reference panel \> 0.2


## Input data

This repository includes dummy study-level GWAS summary statistics and a small subset of reference genotype from 1000G project
that are stored in the `data` directory. 
**These data are for illustration only and are not suitable for running actual analysis**

To reuse the pipeline with other data, please follow the structure within the `data` directory

## Executing the workflow

The worfklow can be executed by using snakemake CLI e.g. 

```
snakemake -c all
```

This will generate files stored on the `results` directory using the given test dataset

To solve package dependencies, the workflow can also utilize [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)
and / or [Apptainer](https://apptainer.org/docs/user/latest/index.html)(formerly Singularity)
e.g. by running

```
snakemake -c all --use-conda --use-singularity
```

Please consult the [`snakemake` documentation](https://snakemake.readthedocs.io/en/stable/index.html) for guidance and additional options.


