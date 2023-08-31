# Rules to perform common tasks & prepare data for downstream analysis

# Get reference genotype file from 1000G

rule get_ref_1000G:
    output:
        "data/geno_ref/1000GP_Phase3_combined.legend.gz"
    shell:
        "wget https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz"

rule format_1000G:
    input:
        "data/geno_ref/1000GP_Phase3_combined.legend.gz"
    output:
        bgz = temp(expand("data/geno_ref/GENO-REF_1000G_{ancestry}.tsv.gz",
                          ancestry = ["EUR", "SAS", "EAS", "AFR", "AMR"])),
        tbi = temp(expand("data/geno_ref/GENO-REF_1000G_{ancestry}.tsv.gz.tbi",
                          ancestry = ["EUR", "SAS", "EAS", "AFR", "AMR"]))
    conda:
        "envs/common.yaml"
    script:
        "scripts/common/format_1000G.R"
