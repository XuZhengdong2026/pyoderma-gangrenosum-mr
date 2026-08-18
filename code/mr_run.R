# MR: 免疫/炎症暴露 -> Pyoderma gangrenosum (FinnGen finn-b-L12_PYODERMA)
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

outcome_id <- "finn-b-L12_PYODERMA"
exposures <- c(
  IBD = "ieu-a-294",
  Neutrophil_count = "ebi-a-GCST90002398",
  Neutrophil_pct = "ebi-a-GCST90002399",
  WBC_count = "ebi-a-GCST90002407"
)

out_dir <- "F:/坏疽性脓皮病/outputs/mr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

run_one <- function(expo_name, expo_id) {
  cat("\n=====", expo_name, "(", expo_id, ") =====\n")
  exp <- tryCatch(
    extract_instruments(outcomes = expo_id, p1 = 5e-8, clump = TRUE),
    error = function(e) { cat("extract_instruments error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(exp) || nrow(exp) == 0) {
    cat("No instruments.\n")
    return(invisible(NULL))
  }
  cat("Instruments:", nrow(exp), "\n")

  outc <- tryCatch(
    extract_outcome_data(snps = exp$SNP, outcomes = outcome_id),
    error = function(e) { cat("extract_outcome error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(outc) || nrow(outc) == 0) {
    cat("No outcome data.\n")
    return(invisible(NULL))
  }

  dat <- harmonise_data(exposure_dat = exp, outcome_dat = outc)
  dat <- dat[dat$mr_keep, ]
  cat("Harmonised SNPs:", nrow(dat), "\n")
  if (nrow(dat) < 3) {
    cat("Too few SNPs for MR.\n")
    return(invisible(NULL))
  }

  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                 "mr_weighted_median", "mr_weighted_mode"))
  res$exposure_name <- expo_name
  het <- mr_heterogeneity(dat)
  ple <- mr_pleiotropy_test(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)

  write.csv(res, file.path(out_dir, paste0(expo_name, "_mr_results.csv")), row.names = FALSE)
  write.csv(het, file.path(out_dir, paste0(expo_name, "_heterogeneity.csv")), row.names = FALSE)
  write.csv(ple, file.path(out_dir, paste0(expo_name, "_pleiotropy.csv")), row.names = FALSE)
  write.csv(single, file.path(out_dir, paste0(expo_name, "_singlesnp.csv")), row.names = FALSE)
  write.csv(loo, file.path(out_dir, paste0(expo_name, "_leaveoneout.csv")), row.names = FALSE)
  write.csv(dat, file.path(out_dir, paste0(expo_name, "_harmonised.csv")), row.names = FALSE)

  pdf(file.path(out_dir, paste0(expo_name, "_forest.pdf")))
  print(mr_forest_plot(single))
  dev.off()
  pdf(file.path(out_dir, paste0(expo_name, "_loo.pdf")))
  print(mr_leaveoneout_plot(loo))
  dev.off()

  print(res[, c("method", "nsnp", "b", "se", "pval")])
  invisible(list(res = res, het = het, ple = ple, dat = dat))
}

all_res <- lapply(seq_along(exposures), function(i) run_one(names(exposures)[i], exposures[i]))
cat("\nAll done. Results in", out_dir, "\n")
