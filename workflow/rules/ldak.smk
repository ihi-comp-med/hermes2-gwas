# Run LDAK to calculate heritability
# Ref: https://dougspeed.com/ldak/

rule dl_ldak_tagging_file:
    """Download LDAK tagging file"""
    output:
        "data/ldak/bld.ldak.hapmap.gbr.tagging"
    shell:
        """
        mkdir -p data/ldak/ && \
        wget -P data/ldak https://genetics.ghpc.au.dk/doug/bld.ldak.hapmap.gbr.tagging.gz && \
        gunzip data/ldak/bld.ldak.hapmap.gbr.tagging.gz
        """



rule ldak_format_sumstats_meta:
    """Format meta-analysis summary stats for LDAK"""
    # use N effective as per https://www.sciencedirect.com/science/article/abs/pii/S0006322322013166?via%3Dihub#!
    input:
        rules.format_metal.output.bgz
    output:
        "results/ldak/{analysis_id}/ldak_sumstats.txt"
    run:
        base = ['A', 'C', 'T', 'G']
        df = pd.read_table(str(input), usecols = ['chr', 'pos_b37', 'A1', 'A2', 'A1_beta', 'pval', 'N_case', 'N_total']) \
               .query("A1 in @base and A2 in @base")
        df['Predictor'] = [':'.join([str(x),str(y)]) for x,y in zip(df['chr'], df['pos_b37'])]
        df.drop_duplicates(subset = ['Predictor'], keep = False, inplace = True)
        df['n'] = 4 / (1/df['N_case'] + 1/(df['N_total'] - df['N_case']))
        df.rename(columns = {'pval': 'P', 'A1_beta': 'Direction'}, inplace=True)
        df[['Predictor', 'A1', 'A2', 'n', 'Direction', 'P']] \
          .to_csv(str(output), index=False, sep=' ')


rule ldak_sumher_meta:
    """
    Wildcards: prev = population prevalence
    """
    input:
        sumstats = rules.ldak_format_sumstats_meta.output,
        tagging_files = lambda x: load_config(x)['ldak']['tagging_file'][x.model]
    output:
        "results/ldak/{analysis_id}/{model}.hers"
    params:
        ldak_exec = "ldak",
        output_prefix = lambda wc, output: strip_ext(wc, output),
        # use n_eff and assume sample prevalence of 0.5
        samp_prev = 0.5,
        pop_prev = lambda x: load_config(x)["ldak"]["pop_prev"]
    container:
        "docker://alhenry/gwaskit"
    shell:
        """
        {params.ldak_exec} \
            --sum-hers {params.output_prefix} \
            --tagfile {input.tagging_files} \
            --summary {input.sumstats} \
            --prevalence {params.pop_prev} \
            --ascertainment {params.samp_prev} \
            --check-sums NO \
            | tee {output}.log
        """
