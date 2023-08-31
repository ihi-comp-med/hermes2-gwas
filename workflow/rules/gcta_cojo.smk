# Rules to identify  genomic susceptibility locus from GWAS summary statistics
# using GCTA-COJO

# -----------------------------------------
# GCTA-COJO: conditional/joint analysis
# -----------------------------------------
def gcta_get_n(wc):
    df = pd.read_table(config['dir_res_meta'] + f'{wc.pheno}_{wc.ancestry}/METAL.N.tsv')
    N_total = sum(df['N_case']) + sum(df['N_control'])
    return N_total

# separate GCTA sumstats per chr
checkpoint gcta_indep_loci_chr:
    input:
        sumstats = rules.format_metal.output.bgz,
        config = "workflow/config/analysis_set/{analysis_id}/locus_annot.yaml",
        # gcta_ma=config['dir_res_post'] + '{pheno}_{ancestry}/GCTA/sumstats.ma',
        # pThresh_info = config['dir_res_meta'] + '{pheno}_{ancestry}/INFO-pThresh.tsv'
    output:
        directory("results/gcta/{analysis_id}/indep_loci_chr")
    params:
        prefix_out = lambda wc,output: Path(str(output[0])).with_suffix('').as_posix(),
        min_pThresh = 1e-5
    run:
        import warnings
        # from scipy.stats import norm
        
        with open(input.config) as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        pThresh = float(config['gcta']['p_threshold'])

        df = pd.read_table(input.sumstats).query('pval < @pThresh')
        
        if len(df.index) == 0:
            warnings.warn(f"WARNING: no variants with p < {pThresh}")
        else:
            Path(str(output)).mkdir(exist_ok = True)
            for i in pd.unique(df.chr):
                outfile = Path(str(output)) / f"chr{i}_sumstats.ma"
                cols = ["SNP", "a1", "a2", "freq", "b", "se", "p", "N"]
                old_cols = ["#key", "A1", "A2", "A1_freq", "A1_beta", "se", "pval", "N_total"]
                d_cols = {k: v for k,v in zip(old_cols, cols)}
                df_chr = df[df.chr == i].rename(columns = d_cols)[cols]
                df_chr.to_csv(outfile, sep = " ", index = False)

def get_bfile_chr(wc):
    """get reference bfile set for gcta cojo"""
    config_file = f"workflow/config/analysis_set/{wc.analysis_id}/locus_annot.yaml"
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    prefix = config['gcta']['prefix_bfile_chr']
    return {'bfile': [prefix + f"{wc.chr}.{ext}" for ext in ['bed', 'bim', 'fam']],
            'p_threshold': config['gcta']['p_threshold'],
            'maf_threshold': config['gcta']['maf_threshold']}

rule gcta_slct_chr:
    input:
        gcta_ma = "results/gcta/{analysis_id}/indep_loci_chr/chr{chr}_sumstats.ma",
        bfile = lambda x: get_bfile_chr(x)['bfile']
    output:
        expand("results/gcta/{{analysis_id}}/indep_loci_chr/chr{{chr}}_slct.{ext}",
               ext = ['log', 'jma.cojo', 'cma.cojo', 'ldr.cojo'])
    threads: 8
    params:
        p_threshold = lambda x: get_bfile_chr(x)['p_threshold'],
        prefix_bfile = lambda x, input: Path(str(input.bfile[0])).with_suffix('').as_posix(),
        prefix_out = lambda x, output: Path(str(output[0])).with_suffix('').as_posix(),
        maf_threshold = lambda x: get_bfile_chr(x)['maf_threshold']
    run:
        import os, tempfile
        
        pThresh = params.p_threshold
        sumstats = pd.read_table(input.gcta_ma, sep = " ")
        bim = pd.read_table(input.bfile[1], header = 0, names = ['chr', 'SNP', 'cm', 'bp', 'a1', 'a2'])

        # check if snps overlap
        if not set(sumstats.SNP).isdisjoint(bim.SNP):
            with tempfile.NamedTemporaryFile() as tmp:
                sumstats.SNP.to_csv(tmp.name, header = False, index = False)
                args = ['gcta64',
                        f'--bfile {params.prefix_bfile}',
                        f'--cojo-file {input.gcta_ma}',
                        f'--cojo-slct',
                        f'--cojo-p {pThresh}',
                        f'--chr {wildcards.chr}',
                        f'--extract {tmp.name}',
                        f'--threads {threads}',
                        f'--maf {params.maf_threshold}',
                        f'--out {params.prefix_out}',
                        f'| tee {output[0]}']
                shell(" ".join(args))
        else:
            shell("touch {output}")

def gcta_slct_aggregate(wc):
    dir_indep_loci_chr = Path(checkpoints.gcta_indep_loci_chr.get(**wc).output[0])
    files_chr = list(dir_indep_loci_chr.glob("*_sumstats.ma"))
    chrs = [f.name.split("_")[0].replace("chr","") for f in files_chr]
    return [(dir_indep_loci_chr / f'chr{i}_slct.jma.cojo').as_posix() for i in chrs]

rule gcta_slct_aggregate:
    input: gcta_slct_aggregate
    output:
        'results/gcta/{analysis_id}/genome-wide_slct-cojo.tsv'
    params:
        bp_window = 500_000
    threads: 8
    script:
        "scripts/gcta/aggregate_gcta_slct.R"