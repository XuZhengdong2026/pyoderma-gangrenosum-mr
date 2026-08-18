# coloc.abf for IBD/CD GWAS vs plasma cis-pQTL
suppressMessages(library(coloc))

input_dir <- "F:/gwas_data/coloc_input"
out_dir <- "F:/gwas_data/coloc_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pairs <- jsonlite::fromJSON(file.path(input_dir, "pairs.json"))

N_GWAS <- c(CD = 40266, IBD = 65642)
S_GWAS <- c(CD = 12194 / 40266, IBD = 31665 / 65642)
N_PQTL <- c(
  Folkersen2020 = 21758, Sun2018 = 3301, Folkersen2018 = 3394,
  Pietzner2020 = 10708
)

summary_rows <- list()

for (i in seq_len(nrow(pairs))) {
  gene <- pairs$gene[i]
  pname <- pairs$pqtl[i]
  gname <- pairs$gwas[i]
  dat <- read.csv(pairs$file[i], stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))
  dat <- dat[!is.na(dat$snp) & dat$snp != "" &
               !is.na(dat$beta_p) & !is.na(dat$beta_g) &
               dat$varbeta_p > 0 & dat$varbeta_g > 0, ]
  if (nrow(dat) < 2) {
    cat(sprintf("%s %s vs %s: too few SNPs (%d)\n",
                gene, pname, gname, nrow(dat)))
    next
  }
  d1 <- list(
    snp = dat$snp, beta = dat$beta_p, varbeta = dat$varbeta_p,
    MAF = dat$maf_p, type = "quant", N = N_PQTL[[pname]]
  )
  d2 <- list(
    snp = dat$snp, beta = dat$beta_g, varbeta = dat$varbeta_g,
    type = "cc", N = N_GWAS[[gname]], s = S_GWAS[[gname]]
  )
  res <- tryCatch(
    coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    cat(sprintf("%s %s vs %s ERROR: %s\n", gene, pname, gname,
                conditionMessage(res)))
    next
  }
  sm <- res$summary
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    gene = gene, pqtl = pname, gwas = gname,
    n_snps = sm[["nsnps"]],
    PP_H0 = sm[["PP.H0.abf"]], PP_H1 = sm[["PP.H1.abf"]],
    PP_H2 = sm[["PP.H2.abf"]], PP_H3 = sm[["PP.H3.abf"]],
    PP_H4 = sm[["PP.H4.abf"]],
    stringsAsFactors = FALSE
  )
  pp <- res$results[, c("snp", "SNP.PP.H4")]
  write.csv(pp, file.path(out_dir, sprintf(
    "pp4_%s_%s_vs_%s.csv", gene, pname, gname)), row.names = FALSE)
  cat(sprintf("%s %s vs %s: n=%d PP.H4=%.4f\n",
              gene, pname, gname, sm[["nsnps"]], sm[["PP.H4.abf"]]))
}

if (length(summary_rows)) {
  sum_df <- do.call(rbind, summary_rows)
  write.csv(sum_df, file.path(out_dir, "coloc_summary.csv"),
            row.names = FALSE)
  cat("saved coloc_summary.csv\n")
}
