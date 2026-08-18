# 两步 MR：IBD -> 免疫/血细胞特征 -> PG（FinnGen）
# 同时输出 IBD->PG 的稳健性检验（MR-PRESSO、Steiger、诊断图）
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ieugwasr)
  library(MRPRESSO)
})

outcome_id <- "finn-b-L12_PYODERMA"
out_dir <- "F:/坏疽性脓皮病/outputs/mr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

# ---------- 1. IBD -> PG 主分析（132 个 LD-clumped IV） ----------
ibd_exp <- read.csv("F:/gwas_data/ibd_ivs_ldclump.csv", stringsAsFactors = FALSE)
cat("IBD LD-clumped IVs:", nrow(ibd_exp), "\n")

ibd_outc <- extract_outcome_chunked(ibd_exp$SNP, outcome_id)
ibd_dat <- harmonise_data(ibd_exp, ibd_outc)
ibd_dat <- ibd_dat[ibd_dat$mr_keep, ]
cat("IBD harmonised SNPs:", nrow(ibd_dat), "\n")

total_res <- mr(ibd_dat, method_list = c("mr_ivw", "mr_egger_regression",
                                         "mr_weighted_median", "mr_weighted_mode"))
total_ivw <- total_res[total_res$method == "Inverse variance weighted", ]
b_total <- total_ivw$b
se_total <- total_ivw$se
write.csv(total_res, file.path(out_dir, "IBD_ldclump_mr_results.csv"), row.names = FALSE)
write.csv(ibd_dat, file.path(out_dir, "IBD_ldclump_harmonised.csv"), row.names = FALSE)

# MR-PRESSO
presso <- tryCatch(
  mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
            SdOutcome = "se.outcome", SdExposure = "se.exposure",
            OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
            data = ibd_dat, NbDistribution = 1000),
  error = function(e) paste("MR-PRESSO error:", conditionMessage(e))
)
if (is.character(presso)) {
  cat(presso, "\n")
} else {
  capture.output(print(presso), file = file.path(out_dir, "IBD_ldclump_mrpresso.txt"))
  saveRDS(presso, file.path(out_dir, "IBD_ldclump_mrpresso.rds"))
  cat("MR-PRESSO done\n")
}

# Steiger 方向性检验
st <- tryCatch(directionality_test(ibd_dat), error = function(e) NULL)
if (!is.null(st)) write.csv(st, file.path(out_dir, "IBD_ldclump_steiger.csv"), row.names = FALSE)

# 诊断图
pdf(file.path(out_dir, "IBD_ldclump_scatter.pdf")); print(mr_scatter_plot(total_res, ibd_dat)); dev.off()
pdf(file.path(out_dir, "IBD_ldclump_forest.pdf")); print(mr_forest_plot(mr_singlesnp(ibd_dat))); dev.off()
pdf(file.path(out_dir, "IBD_ldclump_funnel.pdf")); print(mr_funnel_plot(mr_singlesnp(ibd_dat))); dev.off()
pdf(file.path(out_dir, "IBD_ldclump_loo.pdf")); print(mr_leaveoneout_plot(mr_leaveoneout(ibd_dat))); dev.off()

cat(sprintf("Total effect IBD->PG: b=%.4f se=%.4f p=%.3g OR=%.2f (%.2f-%.2f)\n",
            b_total, se_total, total_ivw$pval, exp(b_total),
            exp(b_total - 1.96 * se_total), exp(b_total + 1.96 * se_total)))

# ---------- 2. 两步 MR ----------
mediators <- list(
  list(name = "Neutrophil_count", id = "ebi-a-GCST90002398", csv = "F:/gwas_data/neutrophil_ivs.csv"),
  list(name = "Neutrophil_pct",   id = "ebi-a-GCST90002399", csv = "F:/gwas_data/NEUTPCT_ivs.csv"),
  list(name = "Monocyte_pct",     id = "ebi-a-GCST90002394", csv = "F:/gwas_data/MONO_PCT_ivs.csv"),
  list(name = "WBC_count",        id = "ebi-a-GCST90002407", csv = "F:/gwas_data/WBC_ivs.csv"),
  list(name = "Platelet_count",   id = "ebi-a-GCST90002402", csv = "F:/gwas_data/PLATELET_ivs.csv")
)

med_rows <- list()
for (m in mediators) {
  cat("\n===== Mediator:", m$name, "=====\n")
  # Step 1: IBD -> mediator
  m_outc <- tryCatch(extract_outcome_chunked(ibd_exp$SNP, m$id), error = function(e) NULL)
  s1 <- NULL
  if (!is.null(m_outc) && nrow(m_outc) > 0) {
    m_dat1 <- harmonise_data(ibd_exp, m_outc)
    m_dat1 <- m_dat1[m_dat1$mr_keep, ]
    if (nrow(m_dat1) >= 3) {
      r1 <- mr(m_dat1, method_list = "mr_ivw")
      s1 <- data.frame(nsnp1 = r1$nsnp, b1 = r1$b, se1 = r1$se, p1 = r1$pval)
    }
  }
  if (is.null(s1)) { cat("Step1 failed/too few SNPs\n"); s1 <- data.frame(nsnp1 = 0, b1 = NA, se1 = NA, p1 = NA) }

  # Step 2: mediator -> PG
  m_exp <- build_mediator_exposure(m$csv, m$name, m$id)
  s2 <- NULL
  if (!is.null(m_exp) && nrow(m_exp) >= 3) {
    pg_outc <- extract_outcome_chunked(m_exp$SNP, outcome_id)
    if (!is.null(pg_outc)) {
      m_exp <- fill_other_allele(m_exp, pg_outc)
      if (nrow(m_exp) >= 3) {
        m_dat2 <- harmonise_data(m_exp, pg_outc)
        m_dat2 <- m_dat2[m_dat2$mr_keep, ]
        if (nrow(m_dat2) >= 3) {
          r2 <- mr(m_dat2, method_list = "mr_ivw")
          s2 <- data.frame(nsnp2 = r2$nsnp, b2 = r2$b, se2 = r2$se, p2 = r2$pval)
          write.csv(m_dat2, file.path(out_dir, paste0("mediator_", m$name, "_harmonised.csv")), row.names = FALSE)
        }
      }
    }
  }
  if (is.null(s2)) { cat("Step2 failed/too few SNPs\n"); s2 <- data.frame(nsnp2 = 0, b2 = NA, se2 = NA, p2 = NA) }

  # 中介比例（delta 法）
  prop <- NA; prop_se <- NA; prop_p <- NA
  if (!is.na(s1$b1) && !is.na(s2$b2) && !is.na(b_total) && b_total != 0) {
    prop <- s1$b1 * s2$b2 / b_total
    var_prop <- (s2$b2^2 * s1$se1^2 + s1$b1^2 * s2$se2^2) / b_total^2
    prop_se <- sqrt(var_prop)
    prop_p <- 2 * pnorm(-abs(prop / prop_se))
  }
  med_rows[[m$name]] <- data.frame(
    mediator = m$name, mediator_id = m$id,
    nsnp1 = s1$nsnp1, b1 = s1$b1, se1 = s1$se1, p1 = s1$p1,
    nsnp2 = s2$nsnp2, b2 = s2$b2, se2 = s2$se2, p2 = s2$p2,
    b_total = b_total, se_total = se_total,
    prop_mediated = prop, prop_se = prop_se, prop_p = prop_p
  )
  cat(sprintf("Step1 b=%.4g p=%.3g | Step2 b=%.4g p=%.3g | prop=%.3f p=%.3g\n",
              s1$b1, s1$p1, s2$b2, s2$p2, prop, prop_p))
}

med_df <- do.call(rbind, med_rows)
write.csv(med_df, file.path(out_dir, "twostep_mr_mediation.csv"), row.names = FALSE)
cat("\nAll done. Results in", out_dir, "\n")
