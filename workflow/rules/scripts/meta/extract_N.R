# snakemake R scripts to extract N from metal info file

pacman::p_load(tidyverse, data.table, glue)

metal_info <- snakemake@input[["metal_info"]]

df_n <- fread(snakemake@input[["study_info"]])

df_metal_info <- fread(
  cmd = glue("grep 'Input File' {metal_info}"),
  sep = "/", header = F,
)[, c(4, 5)] %>%
  separate(V5, c("file", "pheno", "ancestry"), sep = "_") %>%
  select(file, study = V4, pheno, ancestry) %>%
  mutate(ancestry = str_remove_all(ancestry, "\\..*")) %>%
  left_join(df_n, by = c("study", "ancestry", "pheno"))

fwrite(df_metal_info, snakemake@output[[1]], sep = "\t")
