# ----------------
# RULES
# ----------------
# Calculate LDSC using random 10K ukb sample
rule make_symlink_UKB10K:
    input:
        multiext(config['dir_data_UKB10K'] + 'C{chr}_UKBrandom10K_v3_eur_indiv_var_qc_nodupvar',
                 ".bim", ".bed", ".fam")
    output:
        multiext("data/UKB10K_EUR/" + '{chr}', ".bim", ".bed", ".fam")
    run:
        for i, o in zip(input, output):
            Path(o).resolve().symlink_to(Path(i).resolve())

rule calc_ldsc_UKB10K:
    input:
        multiext("data/UKB10K_EUR/" + '{chr}', ".bim", ".bed", ".fam")
    output:
        multiext("data/UKB10K_EUR/" + '{chr}', ".log",
                 ".l2.M", ".l2.M_5_50", ".l2.ldscore.gz")
    conda:
        "envs/ldsc.yaml"
    params:
        prefix = "data/UKB10K_EUR/{chr}"
    shell:
        """
        ./workflow/rules/scripts/ldsc/ldsc.py \
            --bfile {params.prefix} \
            --l2 \
            --ld-wind-cm 1 \
            --out {params.prefix}
        """

# add rsID for study specific formatted results (e.g. UKBiobank)
rule add_rsID_study1:
    input:
        gwas = config['dir_data_qc_ref_check'] + \
            '{study}/REFCHECK-GWAS_{study}_{pheno}_{ref}_{ancestry}.tsv.gz',
        map = config['dir_data_ref'] + \
            'Kaviar-160204-Public-hg19-trim.rsID.vcf.gz'
    output:
        temp(config['dir_res_post'] +
             '{pheno}_{study}-{ref}-{ancestry}/RSID.tsv')
    shell:
        """
        tabix -R <(gunzip -c {input.gwas} | tail -n +2 | cut -f2,3) {input.map} > {output}
        """

rule add_rsID_study2:
    input:
        gwas = config['dir_data_qc_ref_check'] + \
            '{study}/REFCHECK-GWAS_{study}_{pheno}_{ref}_{ancestry}.tsv.gz',
        rsid = config['dir_res_post'] + \
            '{pheno}_{study}-{ref}-{ancestry}/RSID.tsv'
    output:
        temp(config['dir_res_post'] +
             '{pheno}_{study}-{ref}-{ancestry}/REFCHECK-GWAS-RSID.tsv')
    script:
        "scripts/add_rsID_study.R"

# rule to extract HapMap3 SNP for LDSC estimation
rule ldsc_extract_snps:
    input:
        sumstats = config['dir_res_meta'] + \
            '{pheno}_{ancestry}/FORMAT-METAL_{pheno}_{ancestry}.tsv.gz',
        snplist = "data/ldsc/eur_w_ld_chr_old/w_hm3.snplist"
    output:
        temp(config['dir_res_post'] + '{pheno}_{ancestry}/HM3_SNPs.txt')
    shell:
        """
        join --header -1 2 -2 1 -o 1.2,1.5,1.6,1.7,1.10 \
            <(gunzip -c {input.sumstats} | awk 'NR==1 {{print; next}} {{print | "sort -k2,2"}}') \
            <(awk 'NR==1 {{print; next}} {{print | "sort -k1,1"}}' {input.snplist}) > {output}
        """

def ldsc_get_n(wc):
    if wc.ancestry in ['EUR', 'ALL']:
        df = pd.read_table(config['dir_res_meta'] + f'{wc.pheno}_{wc.ancestry}/METAL.N.tsv')
        N_case = sum(df['N_case'])
        N_control = sum(df['N_control'])
    else:
        df = pd.read_table('data/misc/INFO_study_n.tsv')
        study = wc.ancestry.split("-")[0]
        ancestry = wc.ancestry.split("-")[2]
        df_query = df.query('study == @study & ancestry == @ancestry & pheno == @wc.pheno')
        N_case = df_query['N_case'].values[0]
        N_control = df_query['N_control'].values[0]
    return {'N_case': N_case, 'N_control': N_control}

def ldsc_get_sumstats(wc):
    if wc.ancestry in ['EUR', 'ALL']:
        return config['dir_res_meta'] + f'{wc.pheno}_{wc.ancestry}/FORMAT-METAL_{wc.pheno}_{wc.ancestry}.tsv.gz'
    else:
        x = wc.ancestry.split("-")
        return config['dir_res_post'] + f'{wc.pheno}_{x[0]}-{x[1]}-{x[2]}/REFCHECK-GWAS-RSID.tsv'

rule ldsc_munge_sumstats:
    input:
        sumstats = ldsc_get_sumstats
    output:
        config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}.sumstats.gz'
    conda:
        "envs/ldsc.yaml"
    params:
        prefix = config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}',
        N = ldsc_get_n,
        chunksize = 500000
    shell:
        """
        cmdArgs=('--sumstats' '{input.sumstats}' \
            '--out' '{params.prefix}' \
            '--snp' 'rsID' \
            '--a1' 'A1' \
            '--a2' 'A2' \
            '--p' 'pval' \
            '--N-cas' '{params.N[N_case]}' \
            '--N-con' '{params.N[N_control]}' \
            '--chunksize' '{params.chunksize}' \
            '--signed-sumstats' 'A1_beta,0')

        if [ '{wildcards.snpset}' == 'hm3' ]; then
            cmdArgs+=('--merge-allele' 'data/ldsc/eur_w_ld_chr_old/w_hm3.snplist')
        fi

        ./workflow/rules/scripts/ldsc/munge_sumstats.py "${{cmdArgs[@]}}"
        """

def calc_sample_prev(wc):
    if wc.ancestry in ['EUR', 'ALL']:
        file = config['dir_res_meta'] + f'{wc.pheno}_{wc.ancestry}/METAL.N.tsv'
        df = pd.read_table(file)
        sample_prev = sum(df['N_case']) / sum(df['N_total'])
    else:
        df_study = pd.read_table('data/misc/INFO_study_n.tsv')
        study = wc.ancestry.split('-')[0]
        ancestry = wc.ancestry.split('-')[2]
        df_query = df_study.query(
            'study == @study & ancestry == @ancestry & pheno == @wc.pheno')
        sample_prev = df_query['N_case'].values[0] / df_query['N_total'].values[0]
    return sample_prev

def get_pop_prev(wc):
    dict_pheno = {'HF': 'Pheno1', 'CMHF': 'Pheno2', 'CMHFrEF': 'Pheno3',
                  'CMHFpEF': 'Pheno4', 'DCM/HNDC': 'Pheno5', 'LVSD': 'Pheno6',
                  'DCM': 'Pheno5-DCM', 'DCM/HNDC_MTAG=': 'Pheno5-MTAG'}
    sheet_id = "1DfJkqXNOm0Ux_nBESGbzHdWAabbeziLXrjw3VaEMNuQ"
    sheet_name = "hermes2"
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    df = pd.read_csv(url)
    df['trait'].replace(dict_pheno, inplace=True)
    pheno = wc.pheno.split("-")[0]
    return df.pop_prev[df.trait == pheno].values[0]

def ldsc_get_ref_chr(wc):
    """get reference LD panel for ldsc calculation"""
    d = {'UKB10K-EUR': 'data/UKB10K_EUR',
         '1000G-EUR': 'data/ldsc/eur_w_ld_chr_old'}
    return d[wc.ldref]


rule ldsc_calc_h2:
    input:
        sumstats = config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}.sumstats.gz',
        dir_ld_chr = ldsc_get_ref_chr
    output:
        config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{ldref}_{snpset}_h2.log'
    params:
        prefix = lambda wc, output: Path(str(output)).with_suffix('').as_posix(),
        samp_prev = calc_sample_prev,
        pop_prev = get_pop_prev
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        ./workflow/rules/scripts/ldsc/ldsc.py \
            --h2 {input.sumstats} \
            --ref-ld-chr {input.dir_ld_chr}/ \
            --w-ld-chr {input.dir_ld_chr}/ \
            --pop-prev {params.pop_prev} \
            --samp-prev {params.samp_prev} \
            --out {params.prefix}
        """

# Reformat summary statistics for other traits for LDSC calculation
def get_source_sumstats(wc):
    dir_source = Path(config['dir_gwas_sumstats'])
    pattern = f'**/{wc.trait}_{wc.author}*.CLEAN.tsv.gz'
    return list(dir_source.glob(pattern))[0].as_posix()

rule make_symlink_sumstats:
    input:
        get_source_sumstats
    output:
        config['dir_data_gwas'] + 'SUMSTATS_{trait}_{author}.tsv.gz',
        config['dir_data_gwas'] + 'SUMSTATS_{trait}_{author}.tsv.gz.tbi'
    run:
        Path(str(output[0])).resolve().symlink_to(Path(str(input)).resolve())
        Path(str(output[1])).resolve().symlink_to(
            Path(str(input) + '.tbi').resolve())

rule ldsc_extract_snps_other:
    input:
        sumstats = config['dir_data_gwas'] + \
            'SUMSTATS_{trait}_{author}.tsv.gz',
        snplist = "data/ldsc/eur_w_ld_chr_old/w_hm3.snplist"
    output:
        temp(config['dir_data_gwas'] + '{trait}_{author}_HM3_SNPs.txt')
    shell:
        """
        join --header -1 3 -2 1 -o 1.3,1.4,1.5,1.6,1.8 \
            <(gunzip -c {input.sumstats} | awk 'NR==1 {{print; next}} {{print | "sort -k3,3"}}') \
            <(awk 'NR==1 {{print; next}} {{print | "sort -k1,1"}}' {input.snplist}) > {output}
        """

rule ldsc_interstudy_rg_format1:
    input:
        data = lambda wc: list((Path(config['dir_data_qc_ref_check']) / wc.study) \
                                .glob(f"REFCHECK-GWAS_{wc.study}_{wc.pheno}_*_{wc.ancestry}.tsv.gz"))[0],
        info = "data/misc/INFO_study_n.tsv",
        bim_files = expand("data/UKB10K_EUR/{chr}.bim", chr = range(1,23)),
        hm3_b37 = "data/reference/hapmap3_r2_b36_fwd.qc.poly/hapmap3_r3_b36_fwd.CEU.qc.poly.lift_overb37_key.tsv",
        hm3_b37_region = "data/reference/hm3_b37_tabix_query.txt"
    output: temp('data/study_ldsc/{study}/LDSC_{study}_{pheno}_{ancestry}.tsv.gz')
    script: "scripts/format_sumstats_ldsc_study.R"

rule ldsc_interstudy_rg_format2:
    input:
        file = rules.ldsc_interstudy_rg_format1.output,
        snplist = "data/ldsc/eur_w_ld_chr_old/w_hm3.snplist",
        munge_sumstats_script = "workflow/rules/scripts/ldsc/munge_sumstats.py"
    output: 'data/study_ldsc/{study}/LDSC_{study}_{pheno}_{ancestry}.sumstats.gz'
    conda: "envs/ldsc.yaml"
    params:
        prefix = 'data/study_ldsc/{study}/LDSC_{study}_{pheno}_{ancestry}',
        chunksize = 500000
    shell:
        """
        {input.munge_sumstats_script} \
            --sumstats {input.file} \
            --out {params.prefix} \
            --chunksize {params.chunksize} \
            --merge-allele {input.snplist}
        """


def ldsc_interstudy_rg_files(wc):
     min_N = {'Pheno1': 2000,
              'Pheno2': 1000,
              'Pheno3': 500,
              'Pheno4': 500}[wc.pheno]

     df_info = pd.read_table("data/misc/INFO_study_n.tsv") \
                 .query('ancestry == @wc.ancestry & \
                         pheno == @wc.pheno & \
                         N_case > @min_N')

     df_info['samp_prev'] = df_info.apply(lambda x: round(x['N_case'] / x['N_total'], 4), axis = 1)
     df_info['file'] = df_info.apply(lambda x: f'data/study_ldsc/{x.study}/LDSC_{x.study}_{x.pheno}_{x.ancestry}.sumstats.gz',
                          axis = 1)
     df = pd.concat((df_info.query('study == @wc.study'), df_info.query('study != @wc.study')))
     df.reset_index(drop = True, inplace = True)
     return df

rule ldsc_interstudy_rg:
    """
    Calculate genetic correlation between participating study
    """
    input:
        x = lambda wc: ldsc_interstudy_rg_files(wc).iloc[0,:]['file'],
        Y = lambda wc: ldsc_interstudy_rg_files(wc).iloc[1:,:]['file'].tolist(),
        dir_ld_chr = "data/UKB10K_EUR/",
        script = "scripts/ldsc/ldsc.py"
    output:
        log = config['dir_res_post'] + '{pheno}_{ancestry}/interstudy_rg/{study}.log',
        tsv = config['dir_res_post'] + '{pheno}_{ancestry}/interstudy_rg/{study}.tsv'
    params:
        prefix = config['dir_res_post'] + '{pheno}_{ancestry}/interstudy_rg/{study}',
        files = lambda wc, input: ",".join([input.x,*input.Y]),
        pop_prev = lambda wc: ",".join([str(get_pop_prev(wc))] * len(ldsc_interstudy_rg_files(wc).index)),
        samp_prev = lambda wc: ",".join(map(str, ldsc_interstudy_rg_files(wc).samp_prev.tolist()))
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        {input.script} \
            --rg {params.files} \
            --ref-ld-chr {input.dir_ld_chr}/ \
            --w-ld-chr {input.dir_ld_chr}/ \
            --pop-prev {params.pop_prev} \
            --samp-prev {params.samp_prev} \
            --out {params.prefix} && \
        sed -e '1,/Summary of Genetic Correlation/d' \
            {output.log} | head -n -3 > {output.tsv}
        """

def get_n_gwas(wc):
    gsheet_url = "https://docs.google.com/spreadsheets/d/1DfJkqXNOm0Ux_nBESGbzHdWAabbeziLXrjw3VaEMNuQ/"
    df = pd.read_csv(gsheet_url + "export?gid=495134986&format=csv")
    return df.query('trait == @wc.trait & GWAS_author == @wc.author')['n_total'].values[0]

# convert DTAdb to LDSC format
rule ldsc_munge_sumstats_dtadb:
    input:
        # defined in common.smk
        sumstats = get_gwas_sumstats,
        snplist = "data/ldsc/eur_w_ld_chr_old/w_hm3_uni_id.snplist",
        # this is to convert uni_id to rsID
        hm3_vcf = "data/ldsc/eur_w_ld_chr_old/w_hm3_annot.recode.vcf"
    output: config['dir_data_gwas'] + 'LDSC_{trait}_{author}.sumstats.gz'
    conda: "envs/ldsc.yaml"
    params:
        prefix = config['dir_data_gwas'] + 'LDSC_{trait}_{author}',
        N = get_n_gwas,
        chunksize = 500_000,
        snp_col_dtadb = 'uni_id'
    shell:
        """
        ./workflow/rules/scripts/ldsc/munge_sumstats_ext.py \
            --sumstats {input.sumstats} \
            --out {params.prefix} \
            --snp {params.snp_col_dtadb} \
            --a1 effect_allele \
            --a2 other_allele \
            --p pvalue \
            --p-is-logged \
            --N {params.N} \
            --chunksize {params.chunksize} \
            --signed-sumstats effect_size,0 \
            --merge-allele {input.snplist}

        TMP=$(mktemp)
        trap "rm -f $TMP" 0 2 3 15        
        awk -v OFS="\\t" '
        NR == FNR {{
            a1 = $4 < $5 ? $4 : $5
            a2 = $4 < $5 ? $5 : $4
            uni_id = $1"_"$2"_"a1"_"a2
            a[uni_id] = $3
            next }}
        FNR == 1 {{print; next}}
        $1 in a {{print a[$1], $2, $3, $4, $5}}     
        ' <(egrep -v "^#" {input.hm3_vcf}) <(zcat {output}) | \
            gzip -c > $TMP && \
        mv $TMP {output}
        """

rule ldsc_munge_sumstats_other:
    input:
        sumstats = config['dir_data_gwas'] + '{trait}_{author}_HM3_SNPs.txt',
        snplist = "data/ldsc/eur_w_ld_chr_old/w_hm3.snplist"
    output:
        config['dir_data_gwas'] + 'LDSC_{trait}_{author}.sumstats.gz'
    conda:
        "envs/ldsc.yaml"
    params:
        prefix = config['dir_data_gwas'] + 'LDSC_{trait}_{author}',
        N = get_n_gwas,
        chunksize = 500000
    shell:
        """
        ./workflow/rules/scripts/ldsc/munge_sumstats.py \
            --sumstats {input.sumstats} \
            --out {params.prefix} \
            --snp varID \
            --a1 A1 \
            --a2 A2 \
            --p P \
            --N {params.N} \
            --chunksize {params.chunksize} \
            --signed-sumstats beta,0 \
            --merge-allele {input.snplist}
        """

ruleorder: ldsc_munge_sumstats_dtadb > ldsc_munge_sumstats_other

def calc_prev_trait(wc):
    df = pd.read_table(config['dir_data_gwas'] + "GWAS_N.tsv")
    df_trait = df.query('trait == @wc.trait & GWAS_author == @wc.author')
    if df_trait['type'].values[0] == 'quantitative':
        d = {'samp': 'nan', 'pop': 'nan'}
    else:
        d = {'samp': df_trait['sample_prev'].values[0],
             'pop': df_trait['pop_prev'].values[0]}
    return d

def get_trait_sumstats(wc):
    import yaml
    if wc.ancestry in ['EUR', 'ALL']:
        dict_sample_prev = {}
        for f in list(Path(config['dir_res_meta']).glob(f"*_{wc.ancestry}/METAL.N.tsv")):
            pheno = f.parent.stem.split("_")[0]
            df = pd.read_table(f)
            sample_prev = sum(df['N_case']) / sum(df['N_total'])
            dict_sample_prev[pheno] = sample_prev
    else:
        # this is assuming estimation with study data (e.g. UKB)
        df_study = pd.read_table('data/misc/INFO_study_n.tsv')
        study = wc.ancestry.split('-')[0]
        df_query = df_study.query('study == @study')
        dict_sample_prev = dict(zip(df_query['pheno'],
                                  df_query['N_case'] / df_query['N_total']))

    dict_pheno = {'HF': 'Pheno1', 'CMHF': 'Pheno2', 'CMHFrEF': 'Pheno3',
                  'CMHFpEF': 'Pheno4', 'DCM': 'Pheno5', 'LVSD': 'Pheno6'}
    gsheet_url = "https://docs.google.com/spreadsheets/d/1DfJkqXNOm0Ux_nBESGbzHdWAabbeziLXrjw3VaEMNuQ/"
    df_info = pd.read_csv(gsheet_url + "export?gid=495134986&format=csv")
    df_info['trait'].replace(dict_pheno, inplace=True)
    df_info['sample_prev'] = df_info['sample_prev'].fillna(
        df_info['trait'].map(dict_sample_prev))

    # df_info['file'] = [config['dir_res_post'] + f'{x}_{wc.ancestry}/LDSC/{wc.snpset}.sumstats.gz'
    #                    if y == 'HERMES2.0'
    #                    else config['dir_data_gwas'] + f'LDSC_{x}_{y}.sumstats.gz'
    #                    for x, y in zip(df_info['trait'], df_info['GWAS_author'])]

    list_files = [config['dir_res_post'] + f'{wc.pheno}_{wc.ancestry}/LDSC/{wc.snpset}.sumstats.gz']

    # add hermes phenotypes which are not included in the present analysis
    traits = [wc.pheno] + [f'Pheno{x}' for x in range(1,6) if f'Pheno{x}' != wc.pheno]
    list_files = [config['dir_res_post'] + f'{t}_{wc.ancestry}/LDSC/{wc.snpset}.sumstats.gz'
                   for t in traits]
    with open("workflow/config/target_ldsc_traits.yaml") as f:
        target_pheno = yaml.full_load(f)['target_pheno']
        traits += target_pheno
        list_files += [config['dir_data_gwas'] + f"LDSC_{t}.sumstats.gz"
                       for t in target_pheno]

    df_info['trait'] = [f'{x}_{y}' if y != 'HERMES2.0' else x
                        for x,y in zip(df_info['trait'], df_info['GWAS_author'])]
    df_info.set_index('trait', inplace=True)
    df_target = df_info.reindex(traits)
    return {'traits': ','.join(list_files),
            'pop_prev': ','.join([str(round(x, 4)) for x in df_target.pop_prev]),
            'samp_prev': ','.join([str(round(x, 4)) for x in df_target.sample_prev])}


rule ldsc_calc_rg:
    input:
        files = lambda wc: get_trait_sumstats(wc)['traits'].split(','),
        target_sumstats = config['dir_res_post'] + \
            '{pheno}_{ancestry}/LDSC/{snpset}.sumstats.gz',
        dir_ld_chr = "data/UKB10K_EUR/",
        script = "scripts/ldsc/ldsc.py"
    output:
        config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}_rg.log'
    params:
        prefix = config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}_rg',
        trait_info = get_trait_sumstats
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        {input.script} \
            --rg {params.trait_info[traits]} \
            --ref-ld-chr {input.dir_ld_chr}/ \
            --w-ld-chr {input.dir_ld_chr}/ \
            --pop-prev {params.trait_info[pop_prev]} \
            --samp-prev {params.trait_info[samp_prev]} \
            --out {params.prefix}
        """

rule ldsc_get_summary_rg:
    input:
        config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}_rg.log'
    output:
        config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}_rg_summary.tsv'
    shell:
        "sed -e '1,/Summary of Genetic Correlation/d' {input} | head -n -3 > {output}"


# Rules to run LDSC partitioned heritability for cell-type / tissue enrichment analysis
# See https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses for example
rule ldsc_cts_analysis:
    input:
        sumstats = config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}.sumstats.gz',
        dir_ld_chr = "data/ldsc/baseline_v1.2",
        dir_cts_ldscore = "data/ldsc/{cts}_1000Gv3_ldscores",
        ldcts = "data/ldsc/{cts}.ldcts",
        dir_weights = "data/ldsc/1000G_Phase3_weights_hm3_no_MHC"
    output:
        log = config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}-{cts}.log',
        results = config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{snpset}-{cts}.cell_type_results.txt'
    params:
        prefix = lambda wc, output: Path(str(output.log)).with_suffix('').as_posix(),
        dir_cts_ldscore_prefix = lambda wc, input: Path(str(input.dir_cts_ldscore)).name
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        ./workflow/rules/scripts/ldsc/ldsc.py \
            --h2-cts {input.sumstats} \
            --ref-ld-chr {input.dir_ld_chr}/baseline. \
            --out {params.prefix} \
            --ref-ld-chr-cts <(sed s@{params.dir_cts_ldscore_prefix}@{input.dir_cts_ldscore}@g {input.ldcts}) \
            --w-ld-chr {input.dir_weights}/weights.hm3_noMHC.
        """

# Functional enrichment with S-LDSC: https://github.com/omerwe/polyfun/wiki/5.-Estimating-functional-enrichment-using-S-LDSC
# Note: example is in the context of polyfun (but should be the same script under the hood)
rule ldsc_enrichment:
    input:
        dir_ref = "data/ldsc/{annot_set}",
        # sumstats defined in finemap.smks
        sumstats = config['dir_res_post'] + '{pheno}_{ancestry}/finemap/sumstats.parquet',
        script = "scripts/polyfun/ldsc.py"
    output:
        config['dir_res_post'] + '{pheno}_{ancestry}/LDSC/{annot_set}_enrichment.results'
    params:
        output_prefix = lambda wc, output: Path(str(output)).with_suffix('').as_posix(),
        annot_prefix = lambda wc,input: input.dir_ref + "/baselineLF2.2.UKB.",
        weights_prefix = lambda wc,input: input.dir_ref + "/weights.UKB."
    conda:
        "envs/polyfun.yml"
    shell:
        """
        python {input.script} \
            --h2 {input.sumstats} \
            --ref-ld-chr {params.annot_prefix} \
            --w-ld-chr {params.weights_prefix} \
            --out {params.output_prefix} \
            --overlap-annot \
            --not-M-5-50
        """
