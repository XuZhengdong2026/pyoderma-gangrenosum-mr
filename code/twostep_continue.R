# 两步 MR 续跑（断点续传式）：先算 step1，再算/复用 step2，最后汇总中介比例
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ieugwasr)
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

local_clump <- function(exp, kb = 1000) {
  exp <- exp[order(exp$pval), ]
  keep <- logical(nrow(exp))
  for (ch in unique(exp$chr)) {
    idx <- which(exp$chr == ch)
    ord <- idx[order(exp$pos[idx])]
    last <- -Inf
    for (i in ord) {
      if (exp$pos[i] - last >= kb * 1000) { keep[i] <- TRUE; last <- exp$pos[i] }
    }
  }
  exp[keep, ]
}

build_mediator_exposure <- function(csv_path, name, gwas_id) {
  exp0 <- read.csv(csv_path, stringsAsFactors = FALSE)
  exp0 <- exp0[!is.na(exp0$beta) & !is.na(exp0$se) & exp0$pval < 5e-8, ]
  if (nrow(exp0) == 0) return(NULL)
  exp0 <- exp0[order(exp0$pval), ][seq_len(min(nrow(exp0), 400)), ]
  varinfo <- variants_rsid(rsid = exp0$SNP)
  varinfo <- varinfo[, c("name", "chr", "pos")]
  names(varinfo)[1] <- "SNP"
  exp0 <- merge(exp0, varinfo, by = "SNP")
  exp0 <- local_clump(exp0)
  if (nrow(exp0) == 0) return(NULL)
  data.frame(
    SNP = exp0$SNP,
    effect_allele.exposure = exp0$effect_allele,
    other_allele.exposure = NA_character_,
    eaf.exposure = exp0$eaf,
    beta.exposure = exp0$beta,
    se.exposure = exp0$se,
    pval.exposure = exp0$pval,
    exposure = name,
    id.exposure = gwas_id,
    stringsAsFactors = FALSE
  )
}

fill_other_allele <- function(exp_dat, outc) {
  m <- match(exp_dat$SNP, outc$SNP)
  for (i in seq_len(nrow(exp_dat))) {
    j <- m[i]
    if (is.na(j)) next
    oe <- outc$effect_allele.outcome[j]
    oo <- outc$other_allele.outcome[j]
    if (oe == exp_dat$effect_allele.exposure[i]) exp_dat$other_allele.exposure[i] <- oo
    else if (oo == exp_dat$effect_allele.exposure[i]) exp_dat$other_allele.exposure[i] <- oe
  }
  exp_dat[!is.na(exp_dat$other_allele.exposure), ]
}

ibd_exp <- read.csv("F:/gwas_data/ibd_ivs_ldclump.csv", stringsAsFactors = FALSE)
ibd_dat <- read.csv(file.path(out_dir, "IBD_ldclump_harmonised.csv"), stringsAsFactors = FALSE)
ibd_dat <- ibd_dat[ibd_dat$mr_keep, ]
tr <- mr(ibd_dat, method_list = "mr_ivw")
b_total <- tr$b; se_total <- tr$se
cat(sprintf("Total IBD->PG: b=%.4f se=%.4f p=%.3g OR=%.2f\n", b_total, se_total, tr$pval, exp(b_total)))

mediators <- list(
  list(name = "Neutrophil_count", id = "ebi-a-GCST90002398", csv = "F:/gwas_data/neutrophil_ivs.csv"),
  list(name = "Neutrophil_pct",   id = "ebi-a-GCST90002399", csv = "F:/gwas_data/NEUTPCT_ivs.csv"),
  list(name = "Monocyte_pct",     id = "ebi-a-GCST90002394", csv = "F:/gwas_data/MONO_PCT_ivs.csv"),
  list(name = "WBC_count",        id = "ebi-a-GCST90002407", csv = "F:/gwas_data/WBC_ivs.csv"),
  list(name = "Platelet_count",   id = "ebi-a-GCST90002402", csv = "F:/gwas_data/PLATELET_ivs.csv")
)

for (m in mediators) {
  cat("\n=====", m$name, "=====\n")
  s1_path <- file.path(out_dir, paste0("mediator_", m$name, "_step1.csv"))
  if (file.exists(s1_path)) {
    s1 <- read.csv(s1_path)
  } else {
    m_outc <- tryCatch(extract_outcome_chunked(ibd_exp$SNP, m$id), error = function(e) NULL)
    s1 <- data.frame(nsnp1 = 0, b1 = NA, se1 = NA, p1 = NA)
    if (!is.null(m_outc) && nrow(m_outc) > 0) {
      d1 <- harmonise_data(ibd_exp, m_outc)
      d1 <- d1[d1$mr_keep, ]
      if (nrow(d1) >= 3) {
        r1 <- mr(d1, method_list = "mr_ivw")
        s1 <- data.frame(nsnp1 = r1$nsnp, b1 = r1$b, se1 = r1$se, p1 = r1$pval)
      }
    }
    write.csv(s1, s1_path, row.names = FALSE)
  }
  cat("Step1:", s1$nsnp1, "SNPs b=", round(s1$b1, 4), "p=", signif(s1$p1, 3), "\n")

  s2_path <- file.path(out_dir, paste0("mediator_", m$name, "_harmonised.csv"))
  if (file.exists(s2_path)) {
    d2 <- read.csv(s2_path, stringsAsFactors = FALSE)
    d2 <- d2[d2$mr_keep, ]
  } else {
    m_exp <- build_mediator_exposure(m$csv, m$name, m$id)
    d2 <- NULL
    if (!is.null(m_exp) && nrow(m_exp) >= 3) {
      pg_outc <- extract_outcome_chunked(m_exp$SNP, outcome_id)
      if (!is.null(pg_outc)) {
        m_exp <- fill_other_allele(m_exp, pg_outc)
        if (nrow(m_exp) >= 3) {
          d2 <- harmonise_data(m_exp, pg_outc)
          d2 <- d2[d2$mr_keep, ]
          write.csv(d2, s2_path, row.names = FALSE)
        }
      }
    }
  }
  s2 <- data.frame(nsnp2 = 0, b2 = NA, se2 = NA, p2 = NA)
  if (!is.null(d2) && nrow(d2) >= 3) {
    r2 <- mr(d2, method_list = "mr_ivw")
    s2 <- data.frame(nsnp2 = r2$nsnp, b2 = r2$b, se2 = r2$se, p2 = r2$pval)
  }
  cat("Step2:", s2$nsnp2, "SNPs b=", round(s2$b2, 4), "p=", signif(s2$p2, 3), "\n")

  prop <- NA; prop_se <- NA; prop_p <- NA
  if (!is.na(s1$b1) && !is.na(s2$b2) && !is.na(b_total) && b_total != 0) {
    prop <- s1$b1 * s2$b2 / b_total
    var_prop <- (s2$b2^2 * s1$se1^2 + s1$b1^2 * s2$se2^2) / b_total^2
    prop_se <- sqrt(var_prop)
    prop_p <- 2 * pnorm(-abs(prop / prop_se))
  }
  row <- data.frame(
    mediator = m$name, mediator_id = m$id,
    nsnp1 = s1$nsnp1, b1 = s1$b1, se1 = s1$se1, p1 = s1$p1,
    nsnp2 = s2$nsnp2, b2 = s2$b2, se2 = s2$se2, p2 = s2$p2,
    b_total = b_total, se_total = se_total,
    prop_mediated = prop, prop_se = prop_se, prop_p = prop_p
  )
  write.csv(row, file.path(out_dir, paste0("mediator_", m$name, "_mediation.csv")), row.names = FALSE)
}

# 汇总
rows <- lapply(mediators, function(m) {
  read.csv(file.path(out_dir, paste0("mediator_", m$name, "_mediation.csv")))
})
med_df <- do.call(rbind, rows)
write.csv(med_df, file.path(out_dir, "twostep_mr_mediation.csv"), row.names = FALSE)
print(med_df)
cat("\nMediation summary saved.\n")
