# Snakemake rules for fine-mapping using polyfun and polyloc

rule fm_munge_sumstats:
    """munge sumstats using (modified) polyfun helper script"""
    input:
        sumstats = rules.format_metal.output.bgz,
        # slightly modified from the original to handle hermes2 column names
        script = "scripts/polyfun/munge_polyfun_sumstats.py"
    output:
        "results/finemap/{analysis_id}/sumstats.parquet"
    conda:
        "envs/polyfun.yaml"
    params:
        min_info = 0,
        min_maf = 0
    shell:
        """
        python {input.script} \
            --sumstats {input.sumstats} \
            --out {output} \
            --min-info {params.min_info} \
            --min-maf {params.min_maf} \
            --keep-hla
        """
# Approach 1: Using precomputed prior causal probabilities based on a meta-analysis of 15 UK Biobank traits
rule fm_extract_snpvar:
    input:
        sumstats = rules.fm_munge_sumstats.output,
        script = "scripts/polyfun/extract_snpvar.py"
    output:
        'results/finemap/{analysis_id}/snps_with_var.gz'
    conda:
        "envs/polyfun.yaml"
    threads: 8
    shell:
        "python {input.script} --allow-missing --sumstats {input.sumstats} --out {output}"

checkpoint fm_sep_locus:
    input:
        rules.gcta_slct_aggregate.output
    output:
        directory('results/finemap/{analysis_id}/loci')
    run:
        df = pd.read_table(str(input))
        Path(str(output)).mkdir(exist_ok = True)
        for i in pd.unique(df.locus):
            n_var = len(df[df.locus == i].index)
            chr,start,end = df[df.locus == i].reset_index() \
                            .loc[0,["Chr", "start_bp_locus", "end_bp_locus"]]
            outfile = (Path(str(output)) / f"locus{i}.txt").as_posix()
            df[df.locus == i].to_csv(outfile, sep = "\t", index = False)

rule fm_locus_susie:
    input:
        sumstats = rules.fm_extract_snpvar.output,
        locus = 'results/finemap/{analysis_id}/loci/locus{locus}.txt',
        bfiles = lambda wc: expand(load_config(wc)['polyfun']['prefix_bfile_chr'] + '{chr}.{ext}', 
                                   chr = range(1,23), ext = ['bed', 'bim', 'fam']),
        script = "scripts/polyfun/finemapper.py"
    output:
        "results/finemap/{analysis_id}/loci/locus{locus}_finemap.tsv.gz"
    conda:
        "envs/polyfun.yaml"
    threads: 8
    params:
        ld_cache = lambda wc,output: (Path(str(output)).parent / 'LD_cache').as_posix(),
        n_var = 5
    script:
        "scripts/finemap/run_polyfun_susie.py"

def get_baseline_annot(wc):
    finemap_file = f"results/finemap/{wc.analysis_id}/loci/locus{wc.locus}.txt"
    df_locus = pd.read_table(finemap_file)
    chr = df_locus.loc[0,"Chr"]
    annot_file = load_config(wc)['polyfun']['prefix_annot_chr'] + f"{chr}.annot.parquet"
    return annot_file

rule fm_annotate_locus:
    input:
        fm = rules.fm_locus_susie.output,
        annot = get_baseline_annot,
        script = "scripts/polyfun/extract_annotations.py"
    output:
        "results/finemap/{analysis_id}/loci/locus{locus}_annot.tsv.gz"
    conda:
        "envs/polyfun.yaml"
    params:
        pip_cutoff = 0
    shell:
        """
        python {input.script} \
            --pips {input.fm} \
            --annot {input.annot} \
            --pip-cutoff {params.pip_cutoff} \
            --allow-missing \
            --out {output}
        """

def aggregate_fm_locus(wc):
    dir_fm_loci = Path(checkpoints.fm_sep_locus.get(**wc).output[0])
    files = list(dir_fm_loci.glob("*locus*txt"))
    annot_files = [str(f.with_suffix('')) + "_annot.tsv.gz" for f in files]
    fm_files = [str(f.with_suffix('')) + "_finemap.tsv.gz" for f in files]
    return {"annot": annot_files, "fm": fm_files}

rule aggregate_fm_locus:
    input: unpack(aggregate_fm_locus)
    output: "results/finemap/{analysis_id}/credible_set.tsv"
    params:
        pip_percent_cs = lambda x: load_config(x)['polyfun']['pip_percent_cs']
    script: "scripts/finemap/aggregate_fm_locus.R"
