# Save MR-PRESSO output for the R5 (FinnGen R5) CD and UC -> PG analyses
# as reproducible text and RDS files.
suppressMessages(library(MRPRESSO))

cd_harm <- Sys.getenv("CD_HARM")
uc_harm <- Sys.getenv("UC_HARM")
out_txt <- Sys.getenv("OUT_TXT")
out_rds <- Sys.getenv("OUT_RDS")
out_csv <- Sys.getenv("OUT_CSV")

set.seed(20260814)

cat("CD harmonised:", cd_harm, "\n")
cat("UC harmonised:", uc_harm, "\n")

run_presso <- function(label, path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$mr_keep == TRUE, ]
  cat(label, "n =", nrow(d), "\n")
  res <- tryCatch(
    mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      data = d, NbDistribution = 1000, SignifThreshold = 0.05
    ),
    error = function(e) e
  )
  list(label = label, n = nrow(d), res = res)
}

cd_res <- run_presso("CD -> PG (FinnGen R5)", cd_harm)
uc_res <- run_presso("UC -> PG (FinnGen R5)", uc_harm)

extract_row <- function(run) {
  res <- run$res
  if (inherits(res, "error")) {
    return(data.frame(exposure = run$label, n_snp = run$n,
                      raw_b = NA, raw_se = NA, raw_t = NA, raw_p = NA,
                      global_RSSobs = NA, global_P = NA,
                      outliers = "ERROR", stringsAsFactors = FALSE))
  }
  main <- res$`Main MR results`
  raw <- main[main$`MR Analysis` == "Raw", ]
  gt <- res$`MR-PRESSO results`$`Global Test`
  data.frame(
    exposure = sub(" -> PG \\(FinnGen R5\\).*", "", run$label),
    n_snp = run$n,
    raw_b = raw$`Causal Estimate`[1],
    raw_se = raw$Sd[1],
    raw_t = raw$`T-stat`[1],
    raw_p = raw$`P-value`[1],
    global_RSSobs = gt$RSSobs,
    global_P = gt$Pvalue,
    outliers = "none detected",
    stringsAsFactors = FALSE
  )
}

out_df <- rbind(extract_row(cd_res), extract_row(uc_res))
write.csv(out_df, out_csv, row.names = FALSE)

fmt <- function(x) {
  if (inherits(x, "error")) return(paste("ERROR:", conditionMessage(x)))
  paste(capture.output(print(x)), collapse = "\n")
}

txt <- c(
  "MR-PRESSO (1,000 permutations; 0.05 outlier threshold)",
  paste0("CD -> PG (FinnGen R5), harmonised SNPs: ", cd_res$n),
  "------------------------------------------------------------",
  fmt(cd_res$res),
  "",
  paste0("UC -> PG (FinnGen R5), harmonised SNPs: ", uc_res$n),
  "------------------------------------------------------------",
  fmt(uc_res$res),
  ""
)
writeLines(txt, out_txt, useBytes = TRUE)
saveRDS(list(cd = cd_res$res, uc = uc_res$res), out_rds)
cat("WRITTEN", out_txt, "\n")
