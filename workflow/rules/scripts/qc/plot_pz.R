# create pzplot

library(data.table)
library(ggplot2)
library(ragg)

gwas_file <- snakemake@input[["gwas"]]
color <- snakemake@params[["color"]]
decimal <- snakemake@params[["decimal"]]

# for testing
# gwas_file <- "data/qc_ref_check/ARIC-test/REFCHECK-GWAS_ARIC-test_Pheno1_HRC_EUR.tsv.gz"
# color <- list(main = "dodgerblue3", ref_line = "gray60")
# decimal <- 3

study <- paste(
  snakemake@wildcards[["study"]],
  snakemake@wildcards[["pheno"]],
  snakemake@wildcards[["ancestry"]]
)

col_classes <- list(
  character = c("#key", "A1", "A2", "ID"),
  double = c("A1_freq", "MAF", "MAC", "A1_beta", "se", "pval", "A1_freq_ref"),
  integer = c("chr", "pos_b37", "n_events", "n_total", "imputed")
)

DT <- fread(gwas_file,
  colClasses = col_classes,
  select = c("A1_beta", "se", "pval")
)

# quit if no variants
if (nrow(DT) == 0) quit("no")

DT[, p_calculated := 2 * pnorm(-abs(A1_beta / se))]

corr <- round(cor(DT$pval, DT$p_calculated, use = "complete.obs"), 3)

# this is to prune data points for plotting
for (j in c("pval", "p_calculated")) {
  set(DT, j = paste0("log10_", j), value = round(-log10(DT[[j]]), decimal))
}


DT_plot <- DT[!duplicated(DT[, .(log10_pval, log10_p_calculated)])]


plot <- ggplot(DT_plot, aes(y = log10_pval, x = log10_p_calculated)) +
  theme_minimal() +
  #   annotate("text", x = Inf, y = -Inf, hjust = 2, vjust = -2,
  #            label = bquote(italic(r) == .(corr))) +
  geom_abline(
    intercept = 0, slope = 1, linetype = "dashed", color = color$ref_line,
    alpha = 0.8
  ) +
  geom_point(color = color$main, alpha = 0.8, size = 0.8) +
  labs(
    title = study,
    subtitle = bquote(italic(r) == .(corr)),
    y = bquote(-log[10] ~ italic(P)[observed]),
    x = bquote(-log[10] ~ italic(P)[Z - stats])
  ) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  )


par <- snakemake@params

ggsave(snakemake@output[[1]], plot,
  width = par$width,
  height = par$height,
  scaling = par$scaling,
  device = agg_png
)
