# CD-specific two-step mediation + skin-specificity MR controls
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ieugwasr)
})

jwt_lines <- readLines("F:/gwas_data/opengwas_jwt.txt", warn = FALSE)
jwt_lines <- jwt_lines[nzchar(trimws(jwt_lines))]
if (length(jwt_lines) > 0) {
  Sys.setenv(OPENGWAS_JWT = trimws(jwt_lines[1]))
  cat("OpenGWAS JWT loaded\n")
} else {
  stop("OpenGWAS JWT token file is empty")
}

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

load_exposure_from_harmonised <- function(csv_path, label, id) {
  d <- read.csv(csv_path, stringsAsFactors = FALSE)
  d <- d[!duplicated(d$SNP), ]
  data.frame(
    SNP = d$SNP,
    effect_allele.exposure = d$effect_allele.exposure,
    other_allele.exposure = d$other_allele.exposure,
    eaf.exposure = d$eaf.exposure,
    beta.exposure = d$beta.exposure,
    se.exposure = d$se.exposure,
    pval.exposure = d$pval.exposure,
    exposure = label,
    id.exposure = id,
    stringsAsFactors = FALSE
  )
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

# ---------------- CD / UC instruments ----------------
cd_exp <- load_exposure_from_harmonised(
  "F:/坏疽性脓皮病/outputs/mr/CD_deLange_to_PG_harmonised.csv",
  "Crohn's disease", "ebi-a-GCST004132")
uc_exp <- load_exposure_from_harmonised(
  "F:/坏疽性脓皮病/outputs/mr/UC_deLange_to_PG_harmonised.csv",
  "Ulcerative colitis", "ebi-a-GCST004133")
cat("CD instruments:", nrow(cd_exp), "| UC instruments:", nrow(uc_exp), "\n")

# ---------------- 1. CD -> blood-cell traits -> PG mediation ----------------
cd_total <- read.csv("F:/坏疽性脓皮病/outputs/mr/CD_deLange_to_PG_mr_results.csv",
                     stringsAsFactors = FALSE)
cd_ivw <- cd_total[cd_total$method == "Inverse variance weighted", ]
b_total <- cd_ivw$b[1]
se_total <- cd_ivw$se[1]
cat(sprintf("CD->PG total: b=%.4f se=%.4f p=%.4g\n", b_total, se_total, cd_ivw$pval[1]))

mediators <- list(
  list(name = "Neutrophil_count", id = "ebi-a-GCST90002398", csv = "F:/gwas_data/neutrophil_ivs.csv"),
  list(name = "Neutrophil_pct",   id = "ebi-a-GCST90002399", csv = "F:/gwas_data/NEUTPCT_ivs.csv"),
  list(name = "Monocyte_pct",     id = "ebi-a-GCST90002394", csv = "F:/gwas_data/MONO_PCT_ivs.csv"),
  list(name = "WBC_count",        id = "ebi-a-GCST90002407", csv = "F:/gwas_data/WBC_ivs.csv"),
  list(name = "Platelet_count",   id = "ebi-a-GCST90002402", csv = "F:/gwas_data/PLATELET_ivs.csv")
)

cd_med_rows <- list()
for (m in mediators) {
  cat("\n===== CD mediation:", m$name, "=====\n")
  m_outc <- tryCatch(extract_outcome_chunked(cd_exp$SNP, m$id), error = function(e) NULL)
  s1 <- data.frame(nsnp1 = 0, b1 = NA, se1 = NA, p1 = NA)
  if (!is.null(m_outc) && nrow(m_outc) > 0) {
    m_dat1 <- harmonise_data(cd_exp, m_outc)
    m_dat1 <- m_dat1[m_dat1$mr_keep, ]
    if (nrow(m_dat1) >= 3) {
      r1 <- mr(m_dat1, method_list = "mr_ivw")
      s1 <- data.frame(nsnp1 = r1$nsnp, b1 = r1$b, se1 = r1$se, p1 = r1$pval)
      write.csv(m_dat1, file.path(out_dir, paste0("CD_mediator_", m$name, "_step1_harmonised.csv")),
                row.names = FALSE)
    }
  }
  cat(sprintf("Step1 (CD->trait): nsnp=%d b=%.4g p=%.3g\n", s1$nsnp1, s1$b1, s1$p1))

  m_exp <- build_mediator_exposure(m$csv, m$name, m$id)
  s2 <- data.frame(nsnp2 = 0, b2 = NA, se2 = NA, p2 = NA)
  if (!is.null(m_exp) && nrow(m_exp) >= 3) {
    pg_outc <- extract_outcome_chunked(m_exp$SNP, "finn-b-L12_PYODERMA")
    if (!is.null(pg_outc)) {
      m_exp <- fill_other_allele(m_exp, pg_outc)
      if (nrow(m_exp) >= 3) {
        m_dat2 <- harmonise_data(m_exp, pg_outc)
        m_dat2 <- m_dat2[m_dat2$mr_keep, ]
        if (nrow(m_dat2) >= 3) {
          r2 <- mr(m_dat2, method_list = "mr_ivw")
          s2 <- data.frame(nsnp2 = r2$nsnp, b2 = r2$b, se2 = r2$se, p2 = r2$pval)
          write.csv(m_dat2, file.path(out_dir, paste0("CD_mediator_", m$name, "_step2_harmonised.csv")),
                    row.names = FALSE)
        }
      }
    }
  }
  cat(sprintf("Step2 (trait->PG): nsnp=%d b=%.4g p=%.3g\n", s2$nsnp2, s2$b2, s2$p2))

  prop <- NA; prop_se <- NA; prop_p <- NA
  if (!is.na(s1$b1) && !is.na(s2$b2) && !is.na(b_total) && b_total != 0) {
    prop <- s1$b1 * s2$b2 / b_total
    var_prop <- (s2$b2^2 * s1$se1^2 + s1$b1^2 * s2$se2^2) / b_total^2
    prop_se <- sqrt(var_prop)
    prop_p <- 2 * pnorm(-abs(prop / prop_se))
  }
  cd_med_rows[[m$name]] <- data.frame(
    exposure = "CD", mediator = m$name, mediator_id = m$id,
    nsnp1 = s1$nsnp1, b1 = s1$b1, se1 = s1$se1, p1 = s1$p1,
    nsnp2 = s2$nsnp2, b2 = s2$b2, se2 = s2$se2, p2 = s2$p2,
    b_total = b_total, se_total = se_total,
    prop_mediated = prop, prop_se = prop_se, prop_p = prop_p
  )
}
cd_med_df <- do.call(rbind, cd_med_rows)
cd_med_df$prop_fdr <- p.adjust(cd_med_df$prop_p, method = "BH")
cd_med_df$step2_fdr <- p.adjust(cd_med_df$p2, method = "BH")
write.csv(cd_med_df, file.path(out_dir, "CD_twostep_mr_mediation.csv"), row.names = FALSE)
cat("\nCD mediation done\n")

# ---------------- 2. Skin-specificity MR controls ----------------
skin_outcomes <- list(
  list(id = "ebi-a-GCST90019017", label = "Psoriasis (European GWAS)"),
  list(id = "ebi-a-GCST90027161", label = "Atopic dermatitis (European GWAS)"),
  list(id = "finn-b-L12_VITILIGO", label = "Vitiligo (FinnGen R5)"),
  list(id = "finn-b-L12_HIDRADENITISSUP", label = "Hidradenitis suppurativa (FinnGen R5)"),
  list(id = "finn-b-L12_ALOPECAREATA", label = "Alopecia areata (FinnGen R5)"),
  list(id = "finn-b-L12_ERYTHEMANODOSUM", label = "Erythema nodosum (FinnGen R5)")
)

spec_rows <- list()
for (expo in list(list(dat = cd_exp, lab = "CD"), list(dat = uc_exp, lab = "UC"))) {
  for (oc in skin_outcomes) {
    cat("\n====", expo$lab, "->", oc$label, "====\n")
    outc <- tryCatch(extract_outcome_chunked(expo$dat$SNP, oc$id), error = function(e) NULL)
    if (is.null(outc) || nrow(outc) == 0) {
      spec_rows[[length(spec_rows) + 1]] <- data.frame(
        exposure = expo$lab, outcome = oc$label, outcome_id = oc$id,
        nsnp = 0, method = "none", b = NA, se = NA, pval = NA, or = NA, or_lci = NA, or_uci = NA)
      next
    }
    dat <- harmonise_data(expo$dat, outc)
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) < 3) {
      spec_rows[[length(spec_rows) + 1]] <- data.frame(
        exposure = expo$lab, outcome = oc$label, outcome_id = oc$id,
        nsnp = nrow(dat), method = "too_few_snps", b = NA, se = NA, pval = NA, or = NA, or_lci = NA, or_uci = NA)
      next
    }
    res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
    for (i in seq_len(nrow(res))) {
      spec_rows[[length(spec_rows) + 1]] <- data.frame(
        exposure = expo$lab, outcome = oc$label, outcome_id = oc$id,
        nsnp = res$nsnp[i], method = res$method[i],
        b = res$b[i], se = res$se[i], pval = res$pval[i],
        or = exp(res$b[i]), or_lci = exp(res$b[i] - 1.96 * res$se[i]),
        or_uci = exp(res$b[i] + 1.96 * res$se[i]))
    }
    write.csv(dat, file.path(out_dir, paste0("spec_", expo$lab, "_", gsub("[^A-Za-z0-9]", "_", oc$label), "_harmonised.csv")),
              row.names = FALSE)
  }
}
spec_df <- do.call(rbind, spec_rows)
ivw_idx <- spec_df$method == "Inverse variance weighted"
spec_df$fdr_ivw <- NA_real_
spec_df$fdr_ivw[ivw_idx] <- p.adjust(spec_df$pval[ivw_idx], method = "BH")
write.csv(spec_df, file.path(out_dir, "skin_specificity_mr.csv"), row.names = FALSE)

# Forest plot (IVW only)
ivw <- spec_df[spec_df$method == "Inverse variance weighted", ]
if (nrow(ivw) > 0) {
  ivw <- ivw[order(ivw$exposure, ivw$or), ]
  png(file.path(out_dir, "skin_specificity_forest.png"), width = 2200, height = 1600, res = 300)
  par(mar = c(5, 14, 3, 3))
  n <- nrow(ivw)
  plot(NA, xlim = c(0, 3.5), ylim = c(0.5, n + 1.2), xlab = "Odds ratio (95% CI) per SD increase in genetic liability",
       ylab = "", yaxt = "n", main = "CD/UC genetic instruments vs dermatological outcomes")
  abline(v = 1, lty = 2, col = "grey40")
  for (i in seq_len(n)) {
    y <- n - i + 1
    col <- ifelse(ivw$pval[i] < 0.05, "firebrick", "grey30")
    segments(log(ivw$or_lci[i]), y, log(ivw$or_uci[i]), y, col = col, lwd = 2)
    points(log(ivw$or[i]), y, pch = 18, cex = 1.4, col = col)
    text(par("usr")[1], y, adj = 1, labels = sprintf("%s -> %s", ivw$exposure[i], ivw$outcome[i]),
         cex = 0.62, col = "black", xpd = TRUE)
  }
  axis(2, at = seq_len(n), labels = FALSE)
  dev.off()
}
cat("\nSkin specificity MR done. Rows:", nrow(spec_df), "\n")
