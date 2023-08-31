# create pzplot

library(data.table)
library(ggplot2)
library(ragg)

gwas_file <- snakemake@input[["gwas"]]
color <- snakemake@params[["color"]]
decimal <- snakemake@params[["decimal"]]
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
  select = c("pval"), col.names = "p_observed"
)

# quit if no variants
if (nrow(DT) == 0) quit("no")

# calculate expected P-value
DT[order(p_observed), p_expected := ppoints(.N)]

# calculate genomic inflation
lambda_gc <- round(median(qchisq(1 - DT$p_observed, 1)) / qchisq(0.5, 1), 3)

# this is to prune data points for plotting
for (j in c("p_expected", "p_observed")) {
  set(DT, j = paste0("log10_", j), value = round(-log10(DT[[j]]), decimal))
}

DT_plot <- DT[!duplicated(DT[, .(log10_p_expected, log10_p_observed)])]

plot <- ggplot(DT_plot, aes(x = log10_p_expected, y = log10_p_observed)) +
  theme_minimal() +
  # annotate("text", x = Inf, y = -Inf, hjust = 2, vjust = -2,
  #          label = bquote(lambda[GC] == .(lambda_gc))) +
  geom_abline(
    intercept = 0, slope = 1, linetype = "dashed", color = color$ref_line,
    alpha = 0.8
  ) +
  geom_point(color = color$main, alpha = 0.8, size = 0.5) +
  labs(
    title = study,
    x = bquote(-log[10] ~ italic(P)[expected]),
    y = bquote(-log[10] ~ italic(P)[observed]),
    subtitle = bquote(lambda[GC] == .(lambda_gc))
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
