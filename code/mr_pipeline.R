# MR pipeline: UKB 来源暴露 -> FinnGen PG 结局
# 依赖 R 包：TwoSampleMR, MendelianRandomization, ieugwasr, dplyr, readr
#
# 用法（数据就位后）：
#   Rscript mr_pipeline.R --exposure ibd_deLange2017_build37.txt.gz --exposure_type ibd
#   Rscript mr_pipeline.R --exposure 30080_irnt.tsv.bgz --exposure_type neut
# 输出：outputs/mr/<exposure>_mr_results.csv / 森林图 / 漏斗图等

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1]
}

exposure_path <- get_arg("--exposure")
exposure_type <- get_arg("--exposure_type", "ibd")
outcome_path  <- get_arg("--outcome",
  "F:/坏疽性脓皮病/outputs/literature/gwas/finngen_R10_L12_PYODERMA.gz")
out_dir <- "F:/坏疽性脓皮病/outputs/mr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (is.null(exposure_path)) stop("请提供 --exposure 文件路径")

# ---------- 读取暴露 ----------
# IBD (de Lange 2017)：SNP CHR BP A1 A2 INFO N BETA SE P
# Neutrophil (UKB Neale round2)：variant minor_allele minor_AF pval beta se ...
read_exposure <- function(path, type) {
  if (type == "ibd") {
    read_exposure_data(
      filename = path,
      sep = " ",
      snp_col = "SNP", beta_col = "BETA", se_col = "SE",
      effect_allele_col = "A1", other_allele_col = "A2",
      eaf_col = "FRQ", pval_col = "P"
    )
  } else if (type == "neut") {
    # Neale lab UKB: variant=chr:pos_ref_alt, 需要拆出 effect/other allele
    raw <- read.table(gzfile(path), header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, nrows = 1000)
    cat("Neutrophil 表头（前 8 列）：\n")
    print(head(names(raw), 8))
    stop("Neutrophil 读取需按实际表头调整列名（见本文件注释）")
  } else {
    stop("未知 exposure_type")
  }
}

# ---------- 读取结局 (FinnGen R10) ----------
read_finngen_outcome <- function(path) {
  # FinnGen R10: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...
  d <- read.table(gzfile(path), header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, comment.char = "")
  names(d) <- sub("^X\\.chrom$", "chrom", names(d))
  d$SNP <- d$rsids
  d$effect_allele <- d$alt
  d$other_allele <- d$ref
  d$beta <- d$beta
  d$se <- d$sebeta
  d$pval <- d$pval
  d$eaf <- d$af_alt
  d$samplesize <- 429209
  d$outcome <- "Pyoderma gangrenosum (FinnGen R10)"
  d <- d[!is.na(d$SNP) & d$SNP != "", ]
  d[, c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "eaf", "outcome")]
}

# ---------- 主流程 ----------
expo <- read_exposure(exposure_path, exposure_type)
outc <- read_finngen_outcome(outcome_path)

cat("暴露 IV（clump 前）:", nrow(expo), "\n")
expo <- clump_data(expo, pop = "EUR")   # 默认经 OpenGWAS API；无 token 时改用本地 plink
cat("暴露 IV（clump 后）:", nrow(expo), "\n")

dat <- harmonise_data(exposure_dat = expo, outcome_dat = outc)
dat <- dat[dat$mr_keep, ]
cat("成功协调的 SNP 数:", nrow(dat), "\n")

res <- mr(dat, method_list = c(
  "mr_ivw", "mr_egger_regression",
  "mr_weighted_median", "mr_weighted_mode"
))
het <- mr_heterogeneity(dat)
ple <- mr_pleiotropy_test(dat)
single <- mr_singlesnp(dat)

write.csv(res, file.path(out_dir, paste0(exposure_type, "_mr_results.csv")), row.names = FALSE)
write.csv(het, file.path(out_dir, paste0(exposure_type, "_mr_heterogeneity.csv")), row.names = FALSE)
write.csv(ple, file.path(out_dir, paste0(exposure_type, "_mr_pleiotropy.csv")), row.names = FALSE)
write.csv(single, file.path(out_dir, paste0(exposure_type, "_mr_singlesnp.csv")), row.names = FALSE)

pdf(file.path(out_dir, paste0(exposure_type, "_mr_forest.pdf")))
mr_forest_plot(single)
dev.off()

print(res)
cat("MR 结果与图表已输出到:", out_dir, "\n")
