# --------------------------------
# RULES for META-analysis
# --------------------------------

def get_metal_input(wc):
    df = pd.read_table(f"workflow/config/analysis_set/{wc.analysis_id}/meta_studies.tsv")
    cols = ['study', 'pheno', 'ancestry']
    files = df.apply(lambda x: f'results/qc/step2/{x["study"]}/GWAS-QC2_{x["pheno"]}_{x["ancestry"]}.tsv.gz', axis = 1).tolist()    
    files += df.apply(lambda x: f'results/qc/step2/{x["study"]}/QQPLOT_{x["pheno"]}_{x["ancestry"]}.png', axis = 1).tolist()    
    files += df.apply(lambda x: f'results/qc/step2/{x["study"]}/PZPLOT_{x["pheno"]}_{x["ancestry"]}.png', axis = 1).tolist()    
    return files


rule metal_make_config:
    input:
        get_metal_input
    output:
        "results/meta/{analysis_id}/METAL.config"
    params:
        stderr = True,
        avg_freq = True,
        minmax_freq = True,
        verbose = True,
        beta_precision = 12,
        se_precision = 12
    script:
        "scripts/meta/mk_metal_config.py"

rule metal_run:
    """Note: METAL needs to be installed and added to PATH"""
    input:
        rules.metal_make_config.output
    output:
        tsv = "results/meta/{analysis_id}/METAL-GWAS.tsv",
        log = "results/meta/{analysis_id}/METAL-GWAS.log",
        info = "results/meta/{analysis_id}/METAL-GWAS.tsv.info"
    threads: 8
    container:
        "docker://alhenry/docker-gwaskit"
    params:
        outfile = lambda wc,output: Path(str(output.tsv)).with_suffix('').as_posix()
    shell:
        """
        metal {input} | tee {output.log} && \
        mv {params.outfile}1.tsv {output.tsv} && \
        mv {params.outfile}1.tsv.info {output.info}
        """

# extract N based on metal data
rule extract_N_meta:
    input:
        metal_info = rules.metal_run.output.info,
        study_info = "workflow/config/analysis_set/{analysis_id}/meta_studies.tsv"
    output:
        "results/meta/{analysis_id}/METAL-GWAS.N.tsv"
    conda:
        "envs/common.yaml"
    script:
        "scripts/meta/extract_N.R"

rule format_metal:
    """
    Rule to format meta-analysis results from METAL
    & apply post-meta-analysis QC filters
    INPUT: metal -> metal output file
           metal_N -> metal_N output file
           map_rsid -> mapping file for rsID
            should contain at least 3 columns: chr, pos (same genome build), rsID
            such file can be created e.g. from https://db.systemsbiology.net/kaviar/
    OUTPUT: bgzipped & tabix indexed files, md5sum hash
    """
    input:
        metal = rules.metal_run.output.tsv,
        metal_N = rules.extract_N_meta.output,
        map_rsid = "data/utils/MAP_rsID.tsv",
    output:
        bgz = 'results/meta/{analysis_id}/FORMAT-GWAS.tsv.gz',
        tbi = 'results/meta/{analysis_id}/FORMAT-GWAS.tsv.gz.tbi',
        md5 = 'results/meta/{analysis_id}/FORMAT-GWAS.tsv.gz.md5'
    threads: 8
    conda:
        "envs/common.yaml"
    params:
        output_prefix = lambda x, output: Path(str(output[0])).with_suffix("").as_posix(),
        min_maf=0.01,
        min_study=2,
        min_N_prop=0.1 # minimum proportion of total N_effective = 4 / (1/N_case + 1/N_cont)
    script:
        'scripts/meta/format_metal.R'
