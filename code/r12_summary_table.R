suppressPackageStartupMessages({library(dplyr)})
files <- c("IBD", "CD_deLange", "UC_deLange")
rows <- list()
for (f in files) {
  res <- read.csv(paste0("F:/坏疽性脓皮病/outputs/mr/", f, "_R12_mr_results.csv"),
                  stringsAsFactors = FALSE)
  het <- read.csv(paste0("F:/坏疽性脓皮病/outputs/mr/", f, "_R12_heterogeneity.csv"),
                  stringsAsFactors = FALSE)
  ple <- read.csv(paste0("F:/坏疽性脓皮病/outputs/mr/", f, "_R12_pleiotropy.csv"),
                  stringsAsFactors = FALSE)
  ivw <- res[res$method == "Inverse variance weighted", ]
  em  <- res[res$method == "MR Egger", ]
  wm  <- res[res$method == "Weighted median", ]
  q <- het[het$Q_df == max(het$Q_df), ]
  rows[[f]] <- data.frame(
    exposure = f, nsnp = ivw$nsnp,
    OR = round(exp(ivw$b), 3),
    lo = round(exp(ivw$b - 1.96 * ivw$se), 3),
    hi = round(exp(ivw$b + 1.96 * ivw$se), 3),
    ivw_p = signif(ivw$pval, 3),
    egger_p = signif(em$pval, 3),
    wm_p = signif(wm$pval, 3),
    Q = round(q$Q, 2), Q_p = signif(q$Q_pval, 3),
    egger_int_p = signif(ple$pval, 3)
  )
}
out <- do.call(rbind, rows)
print(out)
write.csv(out, "F:/坏疽性脓皮病/outputs/mr/R12_sensitivity_summary.csv", row.names = FALSE)
cat("\n--- MR-PRESSO ---\n")
for (f in files) {
  p <- tryCatch(readLines(paste0("F:/坏疽性脓皮病/outputs/mr/", f, "_R12_mrpresso.txt"),
                          warn = FALSE), error = function(e) NA)
  cat(f, ":", paste(grep("Raw|Global|P-value|Pvalue", p, value = TRUE), collapse = " | "), "\n")
}
