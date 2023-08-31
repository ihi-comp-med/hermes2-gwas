# Aggregate gcta slct results to define loci based on proximity

pacman::p_load(tidyverse, data.table)

input <- snakemake@input
output <- snakemake@output[[1]]
bp_window <- snakemake@params[["bp_window"]]
setDTthreads(snakemake@threads)
# Interactive
# input <- list(gwas = "results/meta/Pheno1_ALL/FORMAT-METAL_Pheno1_ALL.tsv.gz")
# pThresh <- 5e-8
# bp_window <- 500000

DT <- map_df(keep(input, ~ file.size(.x) > 0L)
, fread)[order(Chr, bp)]

DT[,dist_last_locus := c(0, diff(bp)), by = Chr]
DT[,locus := cumsum(dist_last_locus == 0 | dist_last_locus > bp_window)]
DT[,`:=`(start_bp_locus = min(bp) - bp_window,
         end_bp_locus = max(bp) + bp_window),
   by = locus]

fwrite(DT, output, sep = '\t')
