# R script to check allele frequency with reference panel

pacman::p_load(data.table, purrr, magrittr, ggplot2, glue, ragg)

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

setDTthreads(snakemake@threads)

col_classes <- list(
  character = c("#key", "A1", "A2", "ID"),
  double = c("A1_freq", "MAF", "MAC", "A1_beta", "se", "pval", "n_eff"),
  integer = c("chr", "pos_b37", "n_events", "n_total", "imputed")
)

DT_gwas <- fread(input$gwas, key = "#key", colClasses = col_classes)

# use tabix to query
# temp <- tempfile()
#
# DT_gwas[order(chr, pos_b37),2:3] %>% fwrite(temp, sep = '\t', col.names = FALSE)
#
# DT_ref <- fread(cmd = glue("tabix -hR {temp} {args$ref}"), key = "#key", select = c("#key", "MAF", "ID"),
#                 colClasses = list(character = c("#key", "ID"), double = "MAF"))

DT_ref <- fread(input$ref,
  key = "#key", select = c("#key", "A1_freq", "ID"),
  colClasses = list(character = c("#key", "ID"), double = "A1_freq")
)

# left join
# variants not available in DT_ref will be assigned A1_freq_ref = NA
# ID will be replaced with rsID from reference
# DT_gwas[DT_ref, `:=`(ID = i.ID, MAF_ref = i.MAF)]

DT_gwas[DT_ref, `:=`(ID = i.ID, A1_freq_ref = i.A1_freq)]

# filter by MAF and exclude based on AF difference
DT_missing <- DT_gwas[is.na(A1_freq_ref)]

DT_gwas <- DT_gwas[!DT_missing]

# DT_outliers <- DT_gwas[abs(MAF - MAF_ref) > args$max_AF_diff]
DT_outliers <- DT_gwas[abs(A1_freq - A1_freq_ref) > params$max_AF_diff]

DT_output <- DT_gwas[!DT_outliers]

DT_excluded <- rbind(DT_missing, DT_outliers)

# (over)write variant summary statistics
DT_var_stats <- fread(input$var_stats_qc1)
DT_var_stats[, `:=`(
  n_qc2 = nrow(DT_output),
  n_qc2_missing = nrow(DT_missing),
  n_qc2_outliers = nrow(DT_outliers)
)]
fwrite(DT_var_stats, output$var_stats, sep = "\t")

# write to output
fwrite(DT_output, params$prefix_bgz, sep = "\t")
fwrite(DT_excluded, output$excluded, sep = "\t")

# make tabix index
system(paste("bgzip -f", params$prefix_bgz))
system(paste("tabix -f -s2 -b3 -e3", output$bgz))

# ------------------
# make AFCHECK plot
# ------------------
study <- snakemake@wildcards[["study"]]
ancestry <- snakemake@wildcards[["ancestry"]]
ref <- snakemake@config[["ancestry_ref"]][[ancestry]]

title <- paste(
  snakemake@wildcards[["study"]],
  snakemake@wildcards[["pheno"]],
  snakemake@wildcards[["ref"]],
  snakemake@wildcards[["ancestry"]]
)

n_unmatched <- nrow(DT_missing)
n_excluded <- nrow(DT_outliers)

# this is to prune data points for plotting
for (j in c("A1_freq", "A1_freq_ref")) set(DT_gwas, j = j, value = round(DT_gwas[[j]], params$decimal))

DT_plot <- DT_gwas[!duplicated(DT_gwas[, .(A1_freq, A1_freq_ref)])]

# add flag for variants with A1_freq_diff > expected
DT_plot[, A1_freq_diff_flag := abs(A1_freq - A1_freq_ref) > params$max_AF_diff]

corr <- tryCatch(
  {
    round(cor(DT_plot$A1_freq, DT_plot$A1_freq_ref, use = "complete.obs"), 3)
  },
  error = function(e) NA
)

plot <- ggplot(DT_plot, aes(x = A1_freq_ref, y = A1_freq, color = A1_freq_diff_flag)) +
  theme_minimal() +
  # annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -5.5,
  #          label = bquote(italic(r) == .(corr))) +
  # annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -3.5,
  #          label = paste("Variants not in ref =", n_unmatched)) +
  # annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -1.5,
  #          label = paste("Variants excluded =", n_excluded)) +
  geom_abline(intercept = params$max_AF_diff, slope = 1, linetype = "dashed", color = params$color$ref_line) +
  geom_abline(intercept = -params$max_AF_diff, slope = 1, linetype = "dashed", color = params$color$ref_line) +
  geom_point(alpha = 0.8, size = 0.3) +
  coord_cartesian(xlim = c(0, 1.02), ylim = c(0, 1.02), expand = F) +
  scale_color_manual(values = as.character(params$color[1:2])) +
  labs(
    title = title,
    y = paste("EAF", study), x = paste("EAF", ref),
    subtitle = bquote(italic(r) == .(corr) ~ ";" ~
      N[missing] == .(n_unmatched) ~ ";" ~
      N[outlier] == .(n_excluded))
  ) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(output$plot, plot,
  width = params$width,
  height = params$height,
  scaling = params$scaling,
  device = agg_png
)
