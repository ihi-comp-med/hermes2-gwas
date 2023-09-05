import pandas as pd
from pathlib import Path
from math import floor
from snakemake.shell import shell

df = pd.read_table(str(snakemake.input.locus))

max_n = floor(max(df.n))
n_var = snakemake.params.n_var
chr,start,end = df.loc[0,["Chr", "start_bp_locus", "end_bp_locus"]]

ld_cache = (Path(str(snakemake.output)).parent / 'LD_cache').as_posix()
bfile_prefix = Path(snakemake.input.bfiles[0]).parent.as_posix() + f"/chr{chr}"
Path(ld_cache).mkdir(exist_ok = True)

args = ["python",
    "{snakemake.input.script}",
    "--sumstats {snakemake.input.sumstats}",
    "--chr {chr}",
    "--n {max_n}",
    "--start {start}",
    "--end {end}",
    "--method susie",
    "--max-num-causal {n_var}",
    "--allow-missing",
    "--cache-dir {ld_cache}",
    "--out {snakemake.output}"]

if n_var > 1:
    args += ["--geno {bfile_prefix}"]

shell(" ".join(args))