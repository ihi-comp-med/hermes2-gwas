# -----------------------------------
# Snakemake rules to perform GWAS QC
# -----------------------------------

rule qc_step1:
    """
    Implements Step 1 of the GWAS QC pipeline.
    Step 1 includes sanity checks, variant filtering, and genomic control
    """
    input:
        data = "data/studies/{study}/sumstats/GWAS_{pheno}_{ancestry}.tsv.gz"
    output:
        bgz = "results/qc/step1/{study}/GWAS-QC1_{pheno}_{ancestry}.tsv.gz",
        tbi = "results/qc/step1/{study}/GWAS-QC1_{pheno}_{ancestry}.tsv.gz.tbi",
        exclude = "results/qc/step1/{study}/EXCLUDED-VARS-QC1_{pheno}_{ancestry}.tsv.gz",
        varstats = "results/qc/step1/{study}/VAR-STATS-QC1_{pheno}_{ancestry}.tsv"
    # group: "qc_format"
    threads: 8
    params:
        output_prefix = lambda x, output: Path(str(output.bgz)).with_suffix('').as_posix(),
        max_beta = config["qc_step1"]["max_beta"],
        max_se = config["qc_step1"]["max_se"],
        # min_n_events = 200,
        min_maf = config["qc_step1"]["min_maf"],
        min_imp_score = config["qc_step1"]["min_imp_score"],
        min_n_eff = config["qc_step1"]["min_n_eff"],
        max_gc = config["qc_step1"]["max_gc"],
        # others = config["qc_step1"]["others"]
    conda:
        "envs/common.yaml"
    script:
        "scripts/qc/qc_step1.R"


rule qc_step2:
    input:
        gwas = rules.qc_step1.output.bgz,
        var_stats_qc1 = rules.qc_step1.output.varstats,
        ref = lambda x: config["qc_step2"]["afreq_ref"][x.ancestry]
    output:
        bgz = "results/qc/step2/{study}/GWAS-QC2_{pheno}_{ancestry}.tsv.gz",
        tbi = "results/qc/step2/{study}/GWAS-QC2_{pheno}_{ancestry}.tsv.gz.tbi",
        excluded = "results/qc/step2/{study}/EXCLUDED-VARS-QC2_{pheno}_{ancestry}.tsv.gz",
        var_stats = "results/qc/step2/{study}/VAR-STATS_{pheno}_{ancestry}.tsv",
        plot = "results/qc/step2/{study}/AFCHECK_{pheno}_{ancestry}.png"
    params:
        prefix_bgz = lambda wc, output: Path(str(output.bgz)).with_suffix('').as_posix(),
        color = {'main': "gray60",
                 'secondary': "#d62728",
                 'ref_line': "#006ba4"},
        decimal = 3,
        max_AF_diff = 0.2,
        # min_MAF_ref = 0.005,
        width = 2,
        height = 2,
        scaling = 0.5,
        device = "agg_png"
    threads: 8
    conda:
        "envs/common.yaml"
    script:
        "scripts/qc/qc_step2.R"

rule plot_pz:
    input:
        gwas = rules.qc_step2.output.bgz
    output:
        "results/qc/step2/{study}/PZPLOT_{pheno}_{ancestry}.png"
    conda:
        "envs/common.yaml"
    params:
        color = {'main': "#006ba4",
                 'ref_line': "gray60"},
        decimal = 3,
        width = 2,
        height = 2,
        scaling = 0.5,
        device = "agg_png"
    script:
        "scripts/qc/plot_pz.R"

rule plot_qq:
    input:
        gwas = rules.qc_step2.output.bgz
    output:
        "results/qc/step2/{study}/QQPLOT_{pheno}_{ancestry}.png"
    conda:
        "envs/common.yaml"
    params:
        color = {'main': "#006ba4",
                 'ref_line': "gray60"},
        decimal = 3,
        width = 2,
        height = 2,
        scaling = 0.5,
        device = "agg_png"
    script:
        "scripts/qc/plot_qq.R"

