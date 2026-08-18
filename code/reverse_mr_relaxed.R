# Reverse-direction MR: PG -> IBD / CD / UC
# PG instruments selected with relaxed p thresholds (5e-7 / 5e-6 / 5e-5),
# strict LD clumping (r2 < 0.001, 10 Mb, EUR), Steiger filtering.
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ieugwasr)
})

out_dir <- "F:/坏疽性脓皮病/outputs/mr/reverse_relaxed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pg_id <- "finn-b-L12_PYODERMA"
outcomes <- c(
  IBD = "ieu-a-294",
  CD  = "ebi-a-GCST004132",
  UC  = "ebi-a-GCST004133"
)
thresholds <- c(5e-7, 5e-6, 5e-5)

extract_outcome_chunked <- function(snps, id, size = 40) {
  ss <- unique(snps)
  chunks <- split(ss, ceiling(seq_along(ss) / size))
  out <- lapply(chunks, function(ch) {
    tryCatch(extract_outcome_data(snps = ch, outcomes = id), error = function(e) NULL)
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(NULL)
  dplyr::bind_rows(out)
}

for (thr in thresholds) {
  cat("\n################################################\n")
  cat("PG exposure threshold p <", thr, "\n")
  exp <- tryCatch(
    extract_instruments(outcomes = pg_id, p1 = thr, clump = FALSE),
    error = function(e) {
      cat("extract_instruments error:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(exp) || nrow(exp) == 0) {
    cat("No PG instruments.\n")
    next
  }
  cat("Raw PG instruments:", nrow(exp), "\n")

  exp_cl <- tryCatch(
    clump_data(exp, pop = "EUR"),
    error = function(e) {
      cat("clump_data error:", conditionMessage(e), "\n")
      exp
    }
  )
  cat("After strict LD clump (r2<0.001, 10Mb):", nrow(exp_cl), "\n")
  if (nrow(exp_cl) < 3) {
    cat("Too few instruments after clumping.\n")
    next
  }

  fstat <- exp_cl$beta.exposure^2 / exp_cl$se.exposure^2
  cat("F statistics: mean =", round(mean(fstat), 1),
      "min =", round(min(fstat), 1),
      "| n F<10 =", sum(fstat < 10), "\n")

  for (nm in names(outcomes)) {
    outc <- extract_outcome_chunked(exp_cl$SNP, outcomes[nm])
    if (is.null(outc)) {
      cat(nm, ": no outcome data\n")
      next
    }
    dat <- harmonise_data(exp_cl, outc)
    dat <- dat[dat$mr_keep, ]
    cat("\n---", nm, "| harmonised:", nrow(dat), "---\n")
    if (nrow(dat) < 3) {
      cat("Too few harmonised SNPs.\n")
      next
    }
    steiger <- tryCatch(steiger_filtering(dat), error = function(e) NULL)
    if (!is.null(steiger)) {
      dat_s <- steiger[steiger$steiger_dir, ]
      cat("Steiger pass:", nrow(dat_s), "of", nrow(dat), "\n")
      if (nrow(dat_s) >= 3) dat <- dat_s
    }
    res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                   "mr_weighted_median", "mr_weighted_mode"))
    res$threshold <- thr
    res$exposure <- nm
    write.csv(res, file.path(out_dir,
              paste0("PG_to_", nm, "_p", thr, "_mr.csv")), row.names = FALSE)
    write.csv(dat, file.path(out_dir,
              paste0("PG_to_", nm, "_p", thr, "_harmonised.csv")), row.names = FALSE)
    print(res[, c("method", "nsnp", "b", "se", "pval")])
  }
}
cat("\nReverse MR (relaxed) done.\n")
