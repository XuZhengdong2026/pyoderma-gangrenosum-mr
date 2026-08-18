# SuSiE-based multi-signal colocalisation (coloc.susie) for CXCL8 and TNF
suppressMessages(library(coloc))
suppressMessages(library(jsonlite))

out_dir <- "F:/gwas_data/coloc_susie_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

loci <- list(
  CXCL8 = list(
    json = "F:/gwas_data/coloc_susie_CXCL8_Folkersen2020_vs_CD.json",
    Np = 21758, Ng = 40266, sg = 12194 / 40266, label = "CXCL8"
  ),
  TNF = list(
    json = "F:/gwas_data/coloc_susie_TNF_Sun2018_vs_CD.json",
    Np = 3301, Ng = 40266, sg = 12194 / 40266, label = "TNF"
  )
)

for (nm in names(loci)) {
  cat("=====", nm, "=====\n")
  obj <- fromJSON(loci[[nm]]$json, simplifyVector = FALSE)
  rows <- obj$rows
  ld_chunks <- obj$ld
  alleles <- obj$alleles
  # single LD chunk expected
  snplist <- unlist(ld_chunks[[1]]$snplist)
  mat_raw <- ld_chunks[[1]]$matrix
  parsed <- strsplit(snplist, "_")
  rsids <- vapply(parsed, function(x) x[[1]], character(1))
  a1 <- vapply(parsed, function(x) x[[2]], character(1))
  a2 <- vapply(parsed, function(x) x[[3]], character(1))
  LD <- matrix(as.numeric(unlist(mat_raw)), nrow = length(rsids),
               dimnames = list(rsids, rsids))
  cat("LD SNPs:", length(rsids), "\n")

  rdf <- do.call(rbind, lapply(rows, function(r) data.frame(
    snp = r$snp, pos = as.numeric(r$pos),
    beta_p = as.numeric(r$beta_p), varbeta_p = as.numeric(r$varbeta_p),
    maf_p = as.numeric(r$maf_p),
    beta_g = as.numeric(r$beta_g), varbeta_g = as.numeric(r$varbeta_g),
    stringsAsFactors = FALSE
  )))
  keep <- intersect(rsids, rdf$snp)
  rdf <- rdf[rdf$snp %in% keep, ]
  rdf <- rdf[order(rdf$pos), ]
  # align LD reference allele to the pQTL effect allele
  sign <- rep(NA_real_, nrow(rdf))
  for (i in seq_len(nrow(rdf))) {
    rs <- rdf$snp[i]
    al <- alleles[[rs]]
    if (is.null(al)) next
    ea <- al$ea
    if (ea == a1[match(rs, rsids)]) sign[i] <- 1
    else if (ea == a2[match(rs, rsids)]) sign[i] <- -1
  }
  rdf <- rdf[!is.na(sign), ]
  sign <- sign[!is.na(sign)]
  cat("usable SNPs after allele alignment:", nrow(rdf), "\n")
  if (nrow(rdf) < 50) {
    cat("too few SNPs\n")
    next
  }
  D <- diag(sign)
  LD <- LD[rdf$snp, rdf$snp]
  LD <- D %*% LD %*% D
  dimnames(LD) <- list(rdf$snp, rdf$snp)
  rdf$beta_p <- rdf$beta_p * sign
  rdf$beta_g <- rdf$beta_g * sign
  cat("usable SNPs:", nrow(rdf), "\n")

  d1 <- list(snp = rdf$snp, beta = rdf$beta_p, varbeta = rdf$varbeta_p,
             MAF = rdf$maf_p, LD = LD, N = loci[[nm]]$Np, type = "quant")
  d2 <- list(snp = rdf$snp, beta = rdf$beta_g, varbeta = rdf$varbeta_g,
             LD = LD, N = loci[[nm]]$Ng, type = "cc", s = loci[[nm]]$sg)
  res <- tryCatch(
    coloc.susie(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5,
                susie.args = list(L = 10, estimate_residual_variance = FALSE,
                                  estimate_prior_method = "simple")),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    cat("ERROR:", conditionMessage(res), "\n")
    next
  }
  print(res$summary)
  write.csv(as.data.frame(res$summary),
            file.path(out_dir, paste0("susie_summary_", nm, ".csv")),
            row.names = FALSE)
  write.csv(as.data.frame(res$results),
            file.path(out_dir, paste0("susie_results_", nm, ".csv")),
            row.names = FALSE)
  cat("saved", nm, "\n")
}
