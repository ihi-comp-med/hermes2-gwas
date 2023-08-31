# R script to split and format 1000G
pacman::p_load(data.table, tidyverse)
input <- snakemake@input[[1]]
output <- snakemake@output
wildcards <- snakemake@wildcards

# for interactive use / testing

DT <- fread(input)

# the columns AFR, AMR, EAS, EUR, SAS, ALL refer to allle freq for the respective population
# Allele will be sorted based on alphabetical order
DT[, `:=`(
  A1 = pmin(a0, a1),
  A2 = pmax(a0, a1)
)][
  , `:=`(`#key` = paste0(
    formatC(as.numeric(chr), 2, width = 2, format = "d", flag = "0"),
    ":", formatC(position, 2, width = 9, format = "d", flag = "0"),
    ":", A1, "_", A2
  ))
]
setnames(DT, c("position", "id"), c("pos_b37", "ID"))
setkey(DT, `#key`)

for (o in output$bgz) {
  ancestry <- basename(o) %>%
    str_remove_all("\\..*") %>%
    str_split("_") %>%
    pluck(1, 3)
  cols <- c("#key", "chr", "pos_b37", "a0", "a1", "A1", "A2", "ID", ancestry)
  DT_ancestry <- DT[, ..cols]
  setnames(DT_ancestry, ancestry, "A1_freq")

  DT_ancestry[A1 != a1, A1_freq := 1 - A1_freq]
  DT_ancestry[, MAF := pmin(A1_freq, 1 - A1_freq)]

  col_order <- c("#key", "chr", "pos_b37", "A1", "A2", "ID", "A1_freq", "MAF")
  dropcols <- setdiff(colnames(DT_ancestry), col_order)

  DT_ancestry[, c(dropcols) := NULL]
  setcolorder(DT_ancestry, col_order)

  # remove very rare MAF
  output_prefix <- str_remove(o, "\\.gz$")
  fwrite(DT_ancestry[MAF >= 0.001], output_prefix, sep = "\t")

  # create tabix-indexed file
  system2("bgzip", c("-f", output_prefix))
  system2("tabix", c("-f", "-s2", "-b3", "-e3", o))
}
