# Run LDAK to calculate heritability
# Ref: https://dougspeed.com/ldak/

rule ldak_format_sumstats_cohort:
    """Calculate cohort-level summary stats"""
    input:
        rules.format_metal.output.bgz
    output:
        "results/ldak/{analysis_id}/ldak_sumstats.txt"
    shell:
        """
        gunzip -c {input} | awk \
            'BEGIN {{print "Predictor", "A1", "A2", "n", "Direction", "P"}} 
             NR > 1 && $4 ~ /^[ACTG]$/ && $5 ~ /^[ACTG]$/ && !a[$2":"$3]++ {{
                print $2":"$3, $4, $5, $15, $10, $12
                }}' \
             > {output}
        """

def ldak_get_tagging_file(wc):
    """Get tagging file for LDAK"""
    config_file = f"workflow/config/analysis_set/{wc.analysis_id}/locus_annot.yaml"
    with open(config_file, "r") as f:
        config = yaml.full_load(f)
    return config["ldak"][wc.model]

rule ldak_sumher_cohort:
    input:
        sumstats = rules.ldak_format_sumstats_cohort.output,
        tagging_files = ldak_get_tagging_file
    output:
        "results/study_level/{pheno}_{study}_{ref}_{ancestry}/LDAK/{model}-prevmult_{prevmult}.hers"
    params:
        ldak_exec = "ldak",
        output_prefix = lambda wc, output: strip_ext(wc, output),
        # these are defined in ldsc.smk
        samp_prev = 0.5,
        pop_prev = lambda wc: get_pop_prev(wc) * float(str(wc.prevmult))
    container:
        "docker://alhenry/docker-gwaskit:latest"
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

rule ldak_format_sumstats_meta:
    """Format meta-analysis summary stats for LDAK"""
    # use N effective as per https://www.sciencedirect.com/science/article/abs/pii/S0006322322013166?via%3Dihub#!
    input:
        "results/meta/{pheno}_{ancestry}/FORMAT-METAL_{pheno}_{ancestry}.tsv.gz"
    output:
        "results/post_meta/{pheno}_{ancestry}/LDAK/ldak_sumstats.txt"
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

ruleorder: ldak_format_sumstats_cohort > ldak_format_sumstats_meta

def calc_pop_prev_from_sample(wc, input):
    """Calculate population prevalence from sample"""
    df = pd.read_table(input.meta_n)
    return df.N_case.sum() / df.N_total.sum()

rule ldak_sumher_meta:
    """
    prevmult: multiplier for population prevalence(for sensitivity analysis)
    """
    input:
        sumstats = rules.ldak_format_sumstats_meta.output,
        meta_n = "results/meta/{pheno}_{ancestry}/METAL.N.tsv",
        tagging_files = ldak_get_tagging_file
    output:
        "results/post_meta/{pheno}_{ancestry}/LDAK/{model}-prevmult_{prevmult}.hers"
    params:
        ldak_exec = "ldak5.2.linux",
        output_prefix = lambda wc, output: strip_ext(wc, output),
        # these are defined in ldsc.smk - using best guestimate 
        # Doug suggested ascertainment (sample prev) = 0.5 as per https://www.biologicalpsychiatryjournal.com/article/S0006-3223(22)01316-6/fulltext
        samp_prev = 0.5,
        # pop_prev = lambda wc: get_pop_prev(wc) * float(str(wc.prevmult)),
        pop_prev = lambda wc,input: calc_pop_prev_from_sample(wc, input) * float(str(wc.prevmult))
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

ruleorder: ldak_sumher_cohort > ldak_sumher_meta


# rule ldak_gen_cor:
#     """
#     Genetic correlation with LDAK
#     """
#     input:
#         sumstats1 = "results/post_meta/{pheno1}_{ancestry1}/LDAK/ldak_sumstats.txt",
#         sumstats2 = "results/post_meta/{pheno2}_{ancestry2}/LDAK/ldak_sumstats.txt",
#         tagging_files = ldak_get_tagging_file
#     output:
#         "results/ldak/gencor/{pheno1}_{ancestry1}-{pheno2}_{ancestry2}.cors"
#     params:
#         ldak_exec = "ldak5.2.linux",
#         output_prefix = lambda wc, output: strip_ext(wc, output),
#         # these are defined in ldsc.smk
#         # Doug suggested ascertainment (sample prev) = 0.5
#         samp_prev = 0.5,
#         pop_prev = get_pop_prev
#     shell:
#         """
#         {params.ldak_exec} \
#             --sum-cors {params.output_prefix} \
#             --tagfile {input.tagging_files} \
#             --summary {input.sumstats1} \
#             --summary2 {input.sumstats2} \
#             --prevalence {params.pop_prev} \
#             --ascertainment {params.samp_prev} \
#             --allow-ambiguous YES \
#             --check-sums NO \
#             | tee {output}.log
#         """


# rule ldak_collect_results:
#     input:
#         sumstats = expand()