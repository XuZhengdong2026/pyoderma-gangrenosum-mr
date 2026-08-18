# FinnGen R12 (L12_PYODERMA, 703 cases / 470,507 controls) sensitivity MR
# Same exposure instruments as the main analysis; outcome switched to R12.
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

out_dir <- "F:/坏疽性脓皮病/outputs/mr"
r12_csv <- "F:/gwas_data/r12_pg_outcome.csv"

read_r12_outcome <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$effect_allele.outcome <- d$alt38
  d$other_allele.outcome  <- d$ref38
  d$eaf.outcome           <- as.numeric(d$af_alt)
  d$beta.outcome          <- as.numeric(d$beta)
  d$se.outcome            <- as.numeric(d$sebeta)
  d$pval.outcome          <- as.numeric(d$pval)
  d$samplesize.outcome    <- 703 + 470507
  d$outcome               <- "Pyoderma (FinnGen R12)"
  d$id.outcome            <- "finngen_R12_L12_PYODERMA"
  # keep all candidate rows; per-exposure allele matching happens later
  d[, c("SNP", "effect_allele.outcome", "other_allele.outcome",
        "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome",
        "samplesize.outcome", "outcome", "id.outcome")]
}

compl <- function(a) {
  a <- toupper(a)
  ifelse(a == "A", "T", ifelse(a == "T", "A",
         ifelse(a == "C", "G", ifelse(a == "G", "C", NA))))
}

build_exposure <- function(harm_path, name) {
  h <- read.csv(harm_path, stringsAsFactors = FALSE)
  h <- h[h$mr_keep, ]
  h[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
        "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure",
        "samplesize.exposure", "chr.exposure", "pos.exposure",
        "exposure", "id.exposure")]
}

run_r12 <- function(label, harm_path, out_prefix) {
  exp <- build_exposure(harm_path, label)
  outc <- read_r12_outcome(r12_csv)
  # For multi-allelic sites, prefer the outcome row whose alt allele matches
  # the exposure effect allele (or its complement).
  pick <- function(snp, ea, oa) {
    cand <- outc[outc$SNP == snp, ]
    if (nrow(cand) <= 1) return(cand)
    idx <- which(cand$effect_allele.outcome %in% c(ea, compl(oa)) |
                   cand$other_allele.outcome %in% c(ea, compl(oa)))
    if (length(idx) > 0) cand[idx[1], , drop = FALSE] else cand[1, , drop = FALSE]
  }
  outc <- do.call(rbind, lapply(seq_len(nrow(exp)), function(i) {
    pick(exp$SNP[i], exp$effect_allele.exposure[i], exp$other_allele.exposure[i])
  }))
  outc <- outc[!duplicated(outc$SNP), ]
  dat <- harmonise_data(exposure_dat = exp, outcome_dat = outc)
  dat <- dat[dat$mr_keep, ]
  cat("\n===== R12 sensitivity:", label, "=====\n")
  cat("Instruments (main):", nrow(exp),
      "| harmonised with R12:", nrow(dat), "\n")
  if (nrow(dat) < 3) {
    cat("Too few SNPs.\n")
    return(invisible(NULL))
  }
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                 "mr_weighted_median", "mr_weighted_mode"))
  res$exposure_name <- label
  het <- mr_heterogeneity(dat)
  ple <- mr_pleiotropy_test(dat)
  write.csv(res, file.path(out_dir, paste0(out_prefix, "_R12_mr_results.csv")),
            row.names = FALSE)
  write.csv(het, file.path(out_dir, paste0(out_prefix, "_R12_heterogeneity.csv")),
            row.names = FALSE)
  write.csv(ple, file.path(out_dir, paste0(out_prefix, "_R12_pleiotropy.csv")),
            row.names = FALSE)
  write.csv(dat, file.path(out_dir, paste0(out_prefix, "_R12_harmonised.csv")),
            row.names = FALSE)
  print(res[, c("method", "nsnp", "b", "se", "pval")])

  presso <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      data = dat, NbDistribution = 1000, SignifThreshold = 0.05
    ),
    error = function(e) paste("ERR", conditionMessage(e))
  )
  if (is.character(presso)) {
    cat("MR-PRESSO:", presso, "\n")
  } else {
    out <- capture.output(print(presso))
    writeLines(out, file.path(out_dir, paste0(out_prefix, "_R12_mrpresso.txt")))
    cat(paste(out, collapse = "\n"), "\n")
  }
  invisible(list(res = res, het = het, ple = ple, dat = dat))
}

run_r12("IBD -> PG (R12)", 
        "F:/坏疽性脓皮病/outputs/mr/IBD_ldclump_harmonised.csv", "IBD")
run_r12("CD -> PG (R12)",
        "F:/坏疽性脓皮病/outputs/mr/CD_deLange_to_PG_harmonised.csv", "CD_deLange")
run_r12("UC -> PG (R12)",
        "F:/坏疽性脓皮病/outputs/mr/UC_deLange_to_PG_harmonised.csv", "UC_deLange")
cat("\nR12 sensitivity done.\n")
