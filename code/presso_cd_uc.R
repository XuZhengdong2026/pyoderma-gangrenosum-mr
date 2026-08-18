suppressMessages(library(MRPRESSO))

run_presso <- function(label, path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$mr_keep == TRUE, ]
  cat("=== ", label, " n =", nrow(d), " ===\n", sep = "")
  res <- tryCatch(
    mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      data = d, NbDistribution = 1000, SignifThreshold = 0.05
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    cat("ERROR:", conditionMessage(res), "\n")
    return()
  }
  cat("res class:", paste(class(res), collapse = ","), "length:", length(res), "\n")
  cat("res names:", paste(names(res), collapse = " | "), "\n")
  for (i in seq_along(res)) {
    cat("--- element", i, "class:", paste(class(res[[i]]), collapse = ","), "\n")
    if (is.data.frame(res[[i]])) print(res[[i]])
    else if (is.list(res[[i]])) {
      cat("names:", paste(names(res[[i]]), collapse = " | "), "\n")
      print(res[[i]])
    }
  }
}

run_presso("CD", "F:/坏疽性脓皮病/outputs/mr/CD_deLange_to_PG_harmonised.csv")
run_presso("UC", "F:/坏疽性脓皮病/outputs/mr/UC_deLange_to_PG_harmonised.csv")
