# Add Steiger directionality filtering to the relaxed-threshold reverse MR.
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

out_dir <- "F:/坏疽性脓皮病/outputs/mr/reverse_relaxed"
thr_list <- c(5e-7, 5e-6, 5e-5)
outcomes <- c("IBD", "CD", "UC")
pg_n <- 280 + 208449

for (thr in thr_list) {
  for (nm in outcomes) {
    base <- file.path(out_dir, paste0("PG_to_", nm, "_p", thr))
    hpath <- paste0(base, "_harmonised.csv")
    if (!file.exists(hpath)) next
    d <- read.csv(hpath, stringsAsFactors = FALSE)
    d <- d[d$mr_keep, ]
    d$samplesize.exposure <- pg_n
    d$F.exposure <- d$beta.exposure^2 / d$se.exposure^2
    st <- tryCatch(steiger_filtering(d), error = function(e) NULL)
    if (is.null(st)) {
      cat(thr, nm, ": steiger failed, rows =", nrow(d), "\n")
      next
    }
    n_dir <- sum(st$steiger_dir)
    cat(thr, nm, ": harmonised", nrow(st), "| steiger pass", n_dir, "\n")
    write.csv(st, paste0(base, "_steiger.csv"), row.names = FALSE)
    d_f <- st[st$steiger_dir, ]
    if (nrow(d_f) >= 3) {
      res <- mr(d_f, method_list = c("mr_ivw", "mr_egger_regression",
                                     "mr_weighted_median", "mr_weighted_mode"))
      res$threshold <- thr
      write.csv(res, paste0(base, "_steiger_mr.csv"), row.names = FALSE)
      cat("  IVW:", round(res$b[res$method == "Inverse variance weighted"], 4),
          "P =", signif(res$pval[res$method == "Inverse variance weighted"], 3),
          "nsnp =", res$nsnp[res$method == "Inverse variance weighted"], "\n")
    } else {
      cat("  after Steiger <3 SNPs, no MR\n")
    }
  }
}
cat("\nSteiger filtering done.\n")
