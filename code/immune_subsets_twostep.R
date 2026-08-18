# Immune-cell subset two-step MR (CD -> immune subset -> PG)
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ieugwasr)
})

jwt_lines <- readLines("F:/gwas_data/opengwas_jwt.txt", warn = FALSE)
jwt_lines <- jwt_lines[nzchar(trimws(jwt_lines))]
if (length(jwt_lines) > 0) Sys.setenv(OPENGWAS_JWT = trimws(jwt_lines[1]))

out_dir <- "F:/坏疽性脓皮病/outputs/mr/immune_subsets"
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
  pcol <- if ("pval" %in% names(exp)) "pval" else "pval.exposure"
  chrcol <- if ("chr" %in% names(exp)) "chr" else "chr.exposure"
  poscol <- if ("pos" %in% names(exp)) "pos" else "pos.exposure"
  exp <- exp[order(exp[[pcol]]), ]
  keep <- logical(nrow(exp))
  for (ch in unique(exp[[chrcol]])) {
    idx <- which(exp[[chrcol]] == ch)
    ord <- idx[order(exp[[poscol]][idx])]
    last <- -Inf
    for (i in ord) {
      if (exp[[poscol]][i] - last >= kb * 1000) { keep[i] <- TRUE; last <- exp[[poscol]][i] }
    }
  }
  exp[keep, ]
}

cd_exp <- read.csv("F:/坏疽性脓皮病/outputs/mr/CD_deLange_to_PG_harmonised.csv",
                   stringsAsFactors = FALSE)
cd_exp <- cd_exp[!duplicated(cd_exp$SNP), ]
cd_exp <- data.frame(
  SNP = cd_exp$SNP,
  effect_allele.exposure = cd_exp$effect_allele.exposure,
  other_allele.exposure = cd_exp$other_allele.exposure,
  eaf.exposure = cd_exp$eaf.exposure,
  beta.exposure = cd_exp$beta.exposure,
  se.exposure = cd_exp$se.exposure,
  pval.exposure = cd_exp$pval.exposure,
  exposure = "Crohn's disease",
  id.exposure = "ebi-a-GCST004132",
  stringsAsFactors = FALSE
)

cd_ivw <- read.csv("F:/坏疽性脓皮病/outputs/mr/CD_deLange_to_PG_mr_results.csv",
                   stringsAsFactors = FALSE)
cd_ivw <- cd_ivw[cd_ivw$method == "Inverse variance weighted", ]
b_total <- cd_ivw$b[1]
se_total <- cd_ivw$se[1]

traits <- data.frame(
  id = c(
    "ebi-a-GCST90001448", "ebi-a-GCST90001450", "ebi-a-GCST90001477", "ebi-a-GCST90001456",
    "ebi-a-GCST90001513", "ebi-a-GCST90001490", "ebi-a-GCST90001496", "ebi-a-GCST90001502",
    "ebi-a-GCST90001511", "ebi-a-GCST90001504", "ebi-a-GCST90001540", "ebi-a-GCST90001547",
    "ebi-a-GCST90001552", "ebi-a-GCST90001555", "ebi-a-GCST90001544", "ebi-a-GCST90001551",
    "ebi-a-GCST90001548", "ebi-a-GCST90001554", "ebi-a-GCST90001557", "ebi-a-GCST90001560",
    "ebi-a-GCST90001407", "ebi-a-GCST90001409", "ebi-a-GCST90001403", "ebi-a-GCST90001405",
    "ebi-a-GCST90001458", "ebi-a-GCST90001460", "ebi-a-GCST90001461", "ebi-a-GCST90001515",
    "ebi-a-GCST90001524", "ebi-a-GCST90001530", "ebi-a-GCST90001529", "ebi-a-GCST90001532"
  ),
  trait = c(
    "CD11c+ monocyte AC", "CD62L- monocyte AC", "HLA DR++ monocyte AC",
    "CD11c+ HLA DR++ monocyte AC", "CD4 regulatory T cell AC",
    "Resting CD4 regulatory T cell AC", "Activated CD4 regulatory T cell AC",
    "Secreting CD4 regulatory T cell AC", "Activated & secreting CD4 regulatory T cell AC",
    "CD25++ CD4+ T cell AC", "Naive CD4+ T cell AC", "Central memory CD4+ T cell AC",
    "Effector memory CD4+ T cell AC", "Terminally differentiated CD4+ T cell AC",
    "CD45RA- CD4+ T cell AC", "Naive CD8+ T cell AC", "Central memory CD8+ T cell AC",
    "Effector memory CD8+ T cell AC", "Terminally differentiated CD8+ T cell AC",
    "CD45RA+ CD8+ T cell AC", "Memory B cell AC", "Naive-mature B cell AC",
    "Switched memory B cell AC", "Plasma blast/plasma cell AC",
    "Myeloid dendritic cell AC", "Plasmacytoid dendritic cell AC", "Dendritic cell AC",
    "Immature MDSC AC", "Granulocytic MDSC AC", "Monocytic MDSC AC",
    "CD66b++ myeloid cell AC", "Basophil AC"
  ),
  stringsAsFactors = FALSE
)

# remove accidental duplicates and enforce unique ids
traits <- traits[!duplicated(traits$id), ]

rows <- list()
for (k in seq_len(nrow(traits))) {
  tid <- traits$id[k]
  tname <- traits$trait[k]
  cat("\n=====", k, "/", nrow(traits), tname, tid, "=====\n")

  # Step 2: trait instruments -> PG
  inst <- tryCatch(extract_instruments(outcomes = tid, p1 = 5e-8, clump = FALSE),
                   error = function(e) NULL)
  s2 <- data.frame(nsnp2 = 0, b2 = NA, se2 = NA, p2 = NA)
  if (is.data.frame(inst) && nrow(inst) > 0) {
    inst2 <- data.frame(
      SNP = inst$SNP,
      effect_allele.exposure = inst$effect_allele.exposure,
      other_allele.exposure = inst$other_allele.exposure,
      eaf.exposure = inst$eaf.exposure,
      beta.exposure = inst$beta.exposure,
      se.exposure = inst$se.exposure,
      pval.exposure = inst$pval.exposure,
      exposure = inst$exposure,
      id.exposure = inst$id.exposure,
      chr = inst$chr.exposure,
      pos = inst$pos.exposure,
      stringsAsFactors = FALSE
    )
    inst2 <- local_clump(inst2)
    cat("  instruments after clump:", nrow(inst2), "\n")
    if (nrow(inst2) >= 3) {
      pg <- extract_outcome_chunked(inst2$SNP, "finn-b-L12_PYODERMA")
      if (!is.null(pg) && nrow(pg) > 0) {
        dat2 <- harmonise_data(inst2, pg)
        dat2 <- dat2[dat2$mr_keep, ]
        if (nrow(dat2) >= 3) {
          r2 <- mr(dat2, method_list = "mr_ivw")
          s2 <- data.frame(nsnp2 = r2$nsnp, b2 = r2$b, se2 = r2$se, p2 = r2$pval)
          write.csv(dat2, file.path(out_dir, paste0(tid, "_step2_harmonised.csv")),
                    row.names = FALSE)
        }
      }
    }
  }
  cat(sprintf("  Step2 nsnp=%d b=%.4g p=%.3g\n", s2$nsnp2, s2$b2, s2$p2))

  # Step 1: CD -> trait
  s1 <- data.frame(nsnp1 = 0, b1 = NA, se1 = NA, p1 = NA)
  tout <- tryCatch(extract_outcome_chunked(cd_exp$SNP, tid), error = function(e) NULL)
  if (!is.null(tout) && nrow(tout) > 0) {
    dat1 <- harmonise_data(cd_exp, tout)
    dat1 <- dat1[dat1$mr_keep, ]
    if (nrow(dat1) >= 3) {
      r1 <- mr(dat1, method_list = "mr_ivw")
      s1 <- data.frame(nsnp1 = r1$nsnp, b1 = r1$b, se1 = r1$se, p1 = r1$pval)
      write.csv(dat1, file.path(out_dir, paste0(tid, "_step1_harmonised.csv")),
                row.names = FALSE)
    }
  }
  cat(sprintf("  Step1 nsnp=%d b=%.4g p=%.3g\n", s1$nsnp1, s1$b1, s1$p1))

  prop <- NA; prop_se <- NA; prop_p <- NA
  if (!is.na(s1$b1) && !is.na(s2$b2) && !is.na(b_total) && b_total != 0) {
    prop <- s1$b1 * s2$b2 / b_total
    var_prop <- (s2$b2^2 * s1$se1^2 + s1$b1^2 * s2$se2^2) / b_total^2
    prop_se <- sqrt(var_prop)
    prop_p <- 2 * pnorm(-abs(prop / prop_se))
  }
  rows[[k]] <- data.frame(
    trait = tname, trait_id = tid,
    nsnp1 = s1$nsnp1, b1 = s1$b1, se1 = s1$se1, p1 = s1$p1,
    nsnp2 = s2$nsnp2, b2 = s2$b2, se2 = s2$se2, p2 = s2$p2,
    b_total = b_total, se_total = se_total,
    prop_mediated = prop, prop_se = prop_se, prop_p = prop_p
  )
  write.csv(do.call(rbind, rows), file.path(out_dir, "immune_subsets_twostep_results.csv"),
            row.names = FALSE)
}

res <- do.call(rbind, rows)
res$step2_fdr <- p.adjust(res$p2, method = "BH")
res$step1_fdr <- p.adjust(res$p1, method = "BH")
res$prop_fdr <- p.adjust(res$prop_p, method = "BH")
write.csv(res, file.path(out_dir, "immune_subsets_twostep_results.csv"), row.names = FALSE)
cat("\nALL DONE. Rows:", nrow(res), "\n")
