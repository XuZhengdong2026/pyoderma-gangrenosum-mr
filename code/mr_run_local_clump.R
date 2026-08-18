# MR 主流程（本地距离聚类版，避免 API LD clump 超时）
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

outcome_id <- "finn-b-L12_PYODERMA"
out_dir <- "F:/坏疽性脓皮病/outputs/mr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

local_clump <- function(exp, kb = 1000) {
  # 按染色体、按 p 值排序，取每个 1Mb 窗口内最显著 SNP（近似聚类，未考虑 LD r2）
  exp <- exp[order(exp$pval.exposure), ]
  keep <- logical(nrow(exp))
  for (ch in unique(exp$chr.exposure)) {
    idx <- which(exp$chr.exposure == ch)
    pos_sorted <- exp$pos.exposure[idx]
    last_pos <- -Inf
    for (i in idx[order(pos_sorted)]) {
      if (exp$pos.exposure[i] - last_pos >= kb * 1000) {
        keep[i] <- TRUE
        last_pos <- exp$pos.exposure[i]
      }
    }
  }
  exp[keep, ]
}

run_one <- function(expo_name, expo_id) {
  cat("\n=====", expo_name, "(", expo_id, ") =====\n")
  exp <- tryCatch(
    extract_instruments(outcomes = expo_id, p1 = 5e-8, clump = FALSE),
    error = function(e) { cat("extract error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(exp) || nrow(exp) == 0) { cat("No instruments.\n"); return(invisible(NULL)) }
  cat("Raw instruments:", nrow(exp), "\n")
  exp <- local_clump(exp, kb = 1000)
  cat("After local 1Mb clump:", nrow(exp), "\n")
  if (nrow(exp) > 500) exp <- exp[seq_len(500), ]  # 保护 API 配额

  outc <- tryCatch(
    extract_outcome_data(snps = exp$SNP, outcomes = outcome_id),
    error = function(e) { cat("outcome error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(outc) || nrow(outc) == 0) { cat("No outcome data.\n"); return(invisible(NULL)) }

  dat <- harmonise_data(exposure_dat = exp, outcome_dat = outc)
  dat <- dat[dat$mr_keep, ]
  cat("Harmonised SNPs:", nrow(dat), "\n")
  if (nrow(dat) < 3) { cat("Too few SNPs.\n"); return(invisible(NULL)) }

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
  invisible(list(res = res, dat = dat))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  run_one(args[1], args[2])
} else {
  run_one("IBD", "ieu-a-294")
}
cat("\nDone.\n")
