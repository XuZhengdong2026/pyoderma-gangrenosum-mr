# MR: 血细胞计数暴露（GWAS Catalog 显著性位点）-> PG (FinnGen)
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ieugwasr)
})

outcome_id <- "finn-b-L12_PYODERMA"
out_dir <- "F:/坏疽性脓皮病/outputs/mr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

exposures <- list(
  list(name = "Neutrophil_count", id = "GCST90002398", csv = "F:/gwas_data/neutrophil_ivs.csv"),
  list(name = "WBC_count",        id = "GCST90002407", csv = "F:/gwas_data/WBC_ivs.csv"),
  list(name = "Neutrophil_pct",   id = "GCST90002399", csv = "F:/gwas_data/NEUTPCT_ivs.csv")
)

run_one <- function(expo) {
  cat("\n=====", expo$name, "(", expo$id, ") =====\n")
  exp0 <- read.csv(expo$csv, stringsAsFactors = FALSE)
  exp0 <- exp0[!is.na(exp0$beta) & !is.na(exp0$se) & exp0$pval < 5e-8, ]
  cat("GWAS Catalog significant SNPs:", nrow(exp0), "\n")

  varinfo <- ieugwasr::variants_rsid(rsid = exp0$SNP)
  varinfo <- varinfo[, c("name", "chr", "pos")]
  names(varinfo)[1] <- "SNP"
  exp0 <- merge(exp0, varinfo, by = "SNP")
  cat("With positions:", nrow(exp0), "\n")

  exp0 <- exp0[order(exp0$pval), ]
  keep <- logical(nrow(exp0))
  for (ch in unique(exp0$chr)) {
    idx <- which(exp0$chr == ch)
    ord <- idx[order(exp0$pos[idx])]
    last <- -Inf
    for (i in ord) {
      if (exp0$pos[i] - last >= 1e6) { keep[i] <- TRUE; last <- exp0$pos[i] }
    }
  }
  exp0 <- exp0[keep, ]
  cat("After 1Mb clump:", nrow(exp0), "\n")

  exp_dat <- data.frame(
    SNP = exp0$SNP,
    effect_allele.exposure = exp0$effect_allele,
    other_allele.exposure = NA_character_,
    eaf.exposure = exp0$eaf,
    beta.exposure = exp0$beta,
    se.exposure = exp0$se,
    pval.exposure = exp0$pval,
    exposure = expo$name,
    id.exposure = expo$id,
    stringsAsFactors = FALSE
  )

  outc <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id)
  m <- match(exp_dat$SNP, outc$SNP)
  for (i in seq_len(nrow(exp_dat))) {
    j <- m[i]
    if (is.na(j)) next
    oe <- outc$effect_allele.outcome[j]
    oo <- outc$other_allele.outcome[j]
    if (oe == exp_dat$effect_allele.exposure[i]) exp_dat$other_allele.exposure[i] <- oo
    else if (oo == exp_dat$effect_allele.exposure[i]) exp_dat$other_allele.exposure[i] <- oe
  }
  exp_dat <- exp_dat[!is.na(exp_dat$other_allele.exposure), ]
  cat("SNPs with matched alleles:", nrow(exp_dat), "\n")

  dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outc)
  dat <- dat[dat$mr_keep, ]
  cat("Harmonised SNPs:", nrow(dat), "\n")
  if (nrow(dat) < 3) { cat("Too few SNPs.\n"); return(invisible(NULL)) }

  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                 "mr_weighted_median", "mr_weighted_mode"))
  het <- mr_heterogeneity(dat)
  ple <- mr_pleiotropy_test(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)

  write.csv(res, file.path(out_dir, paste0(expo$name, "_gwascatalog_mr_results.csv")), row.names = FALSE)
  write.csv(het, file.path(out_dir, paste0(expo$name, "_gwascatalog_heterogeneity.csv")), row.names = FALSE)
  write.csv(ple, file.path(out_dir, paste0(expo$name, "_gwascatalog_pleiotropy.csv")), row.names = FALSE)
  write.csv(single, file.path(out_dir, paste0(expo$name, "_gwascatalog_singlesnp.csv")), row.names = FALSE)
  write.csv(loo, file.path(out_dir, paste0(expo$name, "_gwascatalog_leaveoneout.csv")), row.names = FALSE)
  write.csv(dat, file.path(out_dir, paste0(expo$name, "_gwascatalog_harmonised.csv")), row.names = FALSE)
  pdf(file.path(out_dir, paste0(expo$name, "_gwascatalog_forest.pdf"))); print(mr_forest_plot(single)); dev.off()
  pdf(file.path(out_dir, paste0(expo$name, "_gwascatalog_loo.pdf"))); print(mr_leaveoneout_plot(loo)); dev.off()

  print(res[, c("method", "nsnp", "b", "se", "pval")])
  invisible(res)
}

for (e in exposures) run_one(e)
cat("\nAll done. Results in", out_dir, "\n")
