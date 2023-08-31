# script for formatting GWAS results
# update 2021-08-25: remove MAF & MAF + INFO parameter as N_eff is already used
pacman::p_load(data.table, purrr, magrittr, stringr)

# Snakemake objects
file_in <- snakemake@input[[1]]
output <- snakemake@output
output_prefix <- snakemake@params[['output_prefix']]
max_beta <- snakemake@params[['max_beta']]
max_se <- snakemake@params[['max_se']]
max_gc <- snakemake@params[['max_gc']]
min_n_eff <- snakemake@params[['min_n_eff']]
min_maf <- snakemake@params[['min_maf']]
min_imp_score <- snakemake@params[['min_imp_score']]
# other_pars <- snakemake@params[['others']]

# for testing / interactive
# file_in <- 'data/ready/Chb/GWAS_Chb_Pheno2_HRC_EUR.tsv.gz'
# output <- list(
#   format = 'data/qc_format/CHS/FORMAT-GWAS_CHS_Pheno1_TOPMED_AFR.tsv.gz',
#   exclude = 'data/qc_format/CHS/EXCLUDE-FORMAT-GWAS_CHS_Pheno1_TOPMED_AFR.tsv.gz',
#   varstats = 'data/qc_format/CHS/VAR-STATS-QC1_CHS_Pheno1_TOPMED_AFR.tsv.gz'
# )
# max_beta <- 10
# max_se <- 10
# min_n_eff <- 50
# min_maf <- 0.01
# min_imp_score <- 0.7
# max_gc <- 1.1
# other_pars <- c("MAF < 0.01 & oevar_imp < 0.9")

key_cols <- c("chr", "pos_b37", "A_coded", "A_noncoded")
DT_in <- fread(file_in, key = key_cols)

# calculate MAF and n_eff
DT_in[, MAF := pmin(AFreq_coded, 1-AFreq_coded)]
DT_in[, n_eff := 2*MAF*(1-MAF)*n_total*oevar_imp]

make_expr <- function(conds) {
  conds %>%
    paste0("(", ., ")") %>%
    paste(collapse = " | ") %>%
    parse(text=.)
}


# Remove missing values
DT <- DT_in[complete.cases(DT_in[,.(chr, pos_b37, A_coded, A_noncoded, beta, SE, pval)])]

# Sanity check
exclude_cond_1 <- c(
  chr = "!chr %in% c(1:22)",
  beta = paste("is.infinite(beta) | abs(beta) >=", max_beta),
  SE =  paste("is.infinite(SE) | SE <= 0 | SE >=", max_se),
  pval = "is.infinite(pval) | pval < 0",
  AFreq_coded = "AFreq_coded <= 0 | AFreq_coded >= 1"
) %>% make_expr()

DT <- DT[!eval(exclude_cond_1)]

# only evaluate if non missing INFO col

if (DT[!is.na(oevar_imp), .N] > 0) {
  # N, MAF, imp_score filters
  exclude_cond_2 <- c(
    "oevar_imp < 0 | oevar_imp > 1",
    paste("n_eff <", min_n_eff),
    # paste("MAF <", min_maf),
    paste("oevar_imp <", min_imp_score)
    # other_pars
  ) %>% make_expr()

  DT <- DT[complete.cases(DT[,.(oevar_imp)])]
} else {
  exclude_cond_2 <- c(
     paste("MAF <", min_maf)
  ) %>% 
  make_expr()
}

DT <- DT[!eval(exclude_cond_2)]

# list excluded variants and write
DT_exclude <- DT_in[!DT]

dir.create(dirname(output$exclude), showWarnings = F, recursive = T)
fwrite(DT_exclude, output$exclude, sep = "\t")

# Sort allele by alphabetical order and create index, format column names
DT[,`:=`(A1 = pmin(A_coded, A_noncoded),
         A2 = pmax(A_coded, A_noncoded),
         A1_freq = AFreq_coded,
         MAC = 2 * MAF * n_total)][
   ,`:=`(`#key` = paste0(formatC(as.numeric(chr), 2, width = 2, format = "d", flag = "0"),
                      ":", formatC(pos_b37, 2, width = 9, format = "d", flag = "0"),
                      ":", A1, "_", A2))]
DT[A1 != A_coded, `:=`(beta = -beta, A1_freq = 1 - AFreq_coded)]

setnames(DT, c('beta', 'SE', 'oevar_imp', 'SNPID'),
         c('A1_beta', 'se', 'imp_score', 'ID'), skip_absent = T)
cols <- c('#key', 'chr', 'pos_b37', 'A1', 'A2', 'ID', 'A1_freq', 'MAF', 'MAC',
          'A1_beta', 'se', 'pval', 'n_events', 'n_total', 'n_eff', 'imputed', 'imp_score')
dropcols <- setdiff(colnames(DT), cols)

DT[, c(dropcols) := NULL]
setcolorder(DT, cols)
setkey(DT, `#key`)

# calculate genomic control and adjust if above specified threshold
lambda_gc <-  median(qchisq(DT$pval, 1, lower.tail = FALSE)) / qchisq(0.5,1)

if (lambda_gc > max_gc){
  DT[, `:=`(se = se * sqrt(lambda_gc),
            pval = pchisq(qchisq(pval, 1, lower.tail = FALSE) / lambda_gc, 1, lower.tail = FALSE))]
}

# write variant summary statistics
DT_var_stats <- basename(file_in) %>% str_remove_all("\\..*") %>%
  str_split("_") %>% .[[1]] %>% .[2:5] %>%
  set_names(c("study", "pheno", "ref", "ancestry")) %>%
  as.list %>%
  as.data.table
DT_var_stats[, `:=`(n_raw = nrow(DT_in),
                    n_qc1 = nrow(DT),
                    n_qc1_exclude = nrow(DT_exclude),
                    lambda_gc = lambda_gc,
                    gc_corrected = lambda_gc > max_gc)]
fwrite(DT_var_stats, output$varstats, sep="\t")

# write file to output
# if (nrow(DT) == 0) quit(save = "no")

DT[,pos_b37 := as.integer(pos_b37)]
fwrite(DT, output_prefix, sep="\t")

# create tabix
system2("bgzip", c("-f", output_prefix))
system2("tabix", c("-f", "-s2", "-b3", "-e3", output$bgz))


# Reduce precision to 4 decimal (necessary?)
# dbl_cols <- colnames(DT)[map_lgl(DT, is.double)]
# for(col in dbl_cols) set(DT, j = col, value = signif(DT[[col]],4))

# DT_bad_col <- imap(list_bad_values,
#                    ~DT_in[eval(.x)][, condition := paste("bad", .y, "value")]) %>%
#   rbindlist %>% setkeyv(key_cols)
