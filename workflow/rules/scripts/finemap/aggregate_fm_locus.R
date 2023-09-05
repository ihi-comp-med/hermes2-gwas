# R script to define credible sets

pacman::p_load(data.table, fs, tidyverse)

MAX_PIP <- as.numeric(snakemake@params[["pip_percent_cs"]]) / 100
FILES_ANNOT <- snakemake@input[["annot"]]
FILES_FM <- snakemake@input[["fm"]]
OUTPUT <- snakemake@output[[1]]

# testing
# MAX_PIP <- 0.95
# FILES_ANNOT <- dir_ls("results/post_meta/Pheno1_EUR/finemap/GWAS-ALL-loci/", glob = "*annot.tsv.gz")
# FILES_FM <- dir_ls("results/post_meta/Pheno1_EUR/finemap/GWAS-ALL-loci/", glob = "*5causal.tsv.gz")

DT_ANNOT <- FILES_ANNOT %>%
  set_names(~basename(.) %>% str_extract("[[0-9]]+")) %>%
  map_df(fread, .id = "locus") %>%
  mutate(locus = as.numeric(locus))

DT_FM <- FILES_FM %>% 
  set_names(~basename(.) %>% str_extract("[[0-9]]+")) %>%
  map_df(fread, .id = "locus") %>%
  mutate(locus = as.numeric(locus))

DT_ANNOT[DT_FM, CREDIBLE_SET := i.CREDIBLE_SET, on = "SNP"]

# get 95 credible sets
DT_ANNOT[, cum_pip := cumsum(PIP), by = locus]

# Write credible sets only
first_cols <- c("locus", "CREDIBLE_SET", "PIP", "cum_pip", "CHR", "SNP", "BP", "A1", "A2")
setcolorder(DT_ANNOT, c(first_cols, setdiff(first_cols, colnames(DT_ANNOT))))
fwrite(DT_ANNOT[CREDIBLE_SET > 0], OUTPUT, sep = "\t")