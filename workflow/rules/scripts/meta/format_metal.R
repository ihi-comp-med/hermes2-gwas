# Script to format METAL output

pacman::p_load(data.table, tidyverse)

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params
setDTthreads(snakemake@threads)

DT <- fread(input$metal)
DT_N <- fread(input$metal_N)

min_N_eff <- params$min_N_prop * 4 / (1 / DT_N[, sum(n_case)] + 1 / DT_N[, sum(n_control)])

# Calculate MAF and N effective
DT[, `:=`(
       MAF = pmin(Freq1, 1 - Freq1),
       N_eff = 4 / (1 / N_case + 1 / N_total)
)]

# if only one study, reduce the minimum study parameter
if (nrow(DT_N) == 1) params$min_study <- 1

DT <- DT[MAF >= params$min_maf & N_eff >= min_N_eff & HetDf >= params$min_study - 1]

df_new <- DT %>%
       separate(MarkerName, c("chr", "pos_b37", "A1_A2"), sep = ":", remove = F, extra = "merge") %>%
       separate(A1_A2, c("A1", "A2"), sep = "_") %>%
       mutate(across(c(chr, pos_b37), as.integer),
              A1_beta = ifelse(A1 == toupper(Allele1), Effect, -Effect),
              A1_freq = ifelse(A1 == toupper(Allele1), Freq1, 1 - Freq1),
              logP = -log10(`P-value`),
              N_control = N_total - N_case,
       ) %>%
       select(
              `#key` = MarkerName, chr, pos_b37, A1, A2, A1_beta,
              A1_freq, se = StdErr, pval = `P-value`, logP,
              N_case, N_control, N_total, isq_het = HetISq, p_het = HetPVal
       )

DT_map <- fread(input$map_rsid,
       na.strings = ".",
       select = 1:3,
       col.names = c("chr", "pos_b37", "rsID"),
       key = c("chr", "pos_b37")
)

setDT(df_new, key = c("chr", "pos_b37"))

df_new[DT_map, rsID := i.rsID]

setkey(df_new, `#key`)

setcolorder(df_new, c("#key", "rsID", setdiff(c("#key", "rsID"), colnames(df_new))))

fwrite(df_new, params$output_prefix, sep = "\t")

system2("bgzip", c("-f", params$output_prefix))
system2("tabix", c("-f", "-s3", "-b4", "-e4", output$bgz))
system2("md5sum", c(output$bgz), stdout = paste0(output$bgz, ".md5"))

# Add rsID
# local copy of map file  is downloaded from KAVIAR grch37 http://db.systemsbiology.net/kaviar/Kaviar.downloads.html
# file_map <- "data/reference/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19-trim.vcf.gz"
# range_query <- GRanges(seqnames = df_new$chr,
#                        ranges = df_new$pos_b37)
# map_cols <- c("chr", "pos_b37", "rsID", "REF", "ALT", "QUAL", "FILTER", "INFO")
#
# df_map <- scanTabix(TabixFile(file_map), param=range_query) %>%
#  fread(col.names = map_cols)
