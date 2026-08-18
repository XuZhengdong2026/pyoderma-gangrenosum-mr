# 敏感性分析：CD/UC -> PG；反向 PG -> IBD；Steiger 方向性检验
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

outcome_id <- "finn-b-L12_PYODERMA"
out_dir <- "F:/坏疽性脓皮病/outputs/mr"

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

run_mr_binary <- function(expo_id, name, outcome_id) {
  cat("\n=====", name, "(", expo_id, ") =====\n")
  exp <- tryCatch(extract_instruments(outcomes = expo_id, p1 = 5e-8, clump = TRUE),
                  error = function(e) paste("ERR", conditionMessage(e)))
  if (is.character(exp)) { cat(exp, "\n"); return(invisible(NULL)) }
  cat("Instruments:", nrow(exp), "\n")
  if (nrow(exp) < 3) return(invisible(NULL))
  outc <- extract_outcome_chunked(exp$SNP, outcome_id)
  if (is.null(outc)) { cat("No outcome data\n"); return(invisible(NULL)) }
  dat <- harmonise_data(exp, outc)
  dat <- dat[dat$mr_keep, ]
  cat("Harmonised:", nrow(dat), "\n")
  if (nrow(dat) < 3) return(invisible(NULL))
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                 "mr_weighted_median", "mr_weighted_mode"))
  write.csv(res, file.path(out_dir, paste0(name, "_mr_results.csv")), row.names = FALSE)
  write.csv(dat, file.path(out_dir, paste0(name, "_harmonised.csv")), row.names = FALSE)
  print(res[, c("method", "nsnp", "b", "se", "pval")])
  invisible(res)
}

# 1. Steiger 方向性检验（主分析）
ibd_dat <- read.csv(file.path(out_dir, "IBD_ldclump_harmonised.csv"), stringsAsFactors = FALSE)
ibd_dat <- ibd_dat[ibd_dat$mr_keep, ]
st <- tryCatch(directionality_test(ibd_dat), error = function(e) NULL)
if (!is.null(st)) {
  write.csv(st, file.path(out_dir, "IBD_ldclump_steiger.csv"), row.names = FALSE)
  print(st)
} else {
  cat("Steiger failed\n")
}

# 2. CD / UC -> PG
run_mr_binary("ebi-a-GCST004132", "CD_deLange_to_PG", outcome_id)
run_mr_binary("ebi-a-GCST004133", "UC_deLange_to_PG", outcome_id)

# 3. 反向 MR：PG -> IBD（若 PG 有 >=3 个 IV）
pg_exp <- tryCatch(extract_instruments(outcomes = outcome_id, p1 = 5e-8, clump = TRUE),
                   error = function(e) paste("ERR", conditionMessage(e)))
if (is.character(pg_exp)) {
  cat("PG instruments error:", pg_exp, "\n")
} else {
  cat("PG instruments:", nrow(pg_exp), "\n")
  if (!is.null(pg_exp) && nrow(pg_exp) >= 3) {
    ibd_outc <- extract_outcome_chunked(pg_exp$SNP, "ieu-a-294")
    if (!is.null(ibd_outc)) {
      dat <- harmonise_data(pg_exp, ibd_outc)
      dat <- dat[dat$mr_keep, ]
      if (nrow(dat) >= 3) {
        res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
        write.csv(res, file.path(out_dir, "PG_to_IBD_mr_results.csv"), row.names = FALSE)
        print(res[, c("method", "nsnp", "b", "se", "pval")])
      }
    }
  }
}
cat("\nSensitivity done.\n")
