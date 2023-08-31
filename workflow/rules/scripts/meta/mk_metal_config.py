# script to write METAL config file and run
from snakemake.shell import shell
from pathlib import Path

# write config file
with open(str(snakemake.output), 'a+') as o:
    if snakemake.params.stderr is True:
        o.write('\nSCHEME\tSTDERR\n')
    if snakemake.params.avg_freq is True:
        o.write('\nAVERAGEFREQ\tON\n')
    if snakemake.params.minmax_freq is True:
        o.write('\nMINMAXFREQ\tON\n')

    o.write(f"""
CUSTOMVARIABLE  N_case
CUSTOMVARIABLE  N_total
""")

    for i in snakemake.input:
        file = str(i)
        config=f"""
MARKER  #key
ALLELE  A1 A2
EFFECT  A1_beta
FREQ    A1_freq
PVALUE  pval
STDERR  se
COLUMNCOUNTING LENIENT
LABEL N_case as n_events
LABEL N_total as n_total
PROCESS {file}
"""
        o.write(config)

    analysis_id = snakemake.wildcards.analysis_id
    outfile = f'results/meta/{analysis_id}/METAL-GWAS'

    if snakemake.params.verbose is True:
        o.write('\nVERBOSE\tON\n')

    o.write(f'\nEFFECT_PRINT_PRECISION\t{snakemake.params.beta_precision}\n')
    o.write(f'\nSTDERR_PRINT_PRECISION\t{snakemake.params.se_precision}\n')
    o.write(f'\nOUTFILE {outfile} .tsv\nANALYZE HETEROGENEITY\n')
