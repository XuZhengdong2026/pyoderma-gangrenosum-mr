# Map GPL6244 probe IDs to gene symbols using hugene10sttranscriptcluster.db
suppressPackageStartupMessages({
  library(hugene10sttranscriptcluster.db)
  library(AnnotationDbi)
})

probes <- readLines("F:/坏疽性脓皮病/outputs/mr/gut_skin/probe_ids.txt", warn = FALSE)
probes <- probes[nzchar(probes)]
map <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                             keys = probes,
                             columns = "SYMBOL",
                             keytype = "PROBEID")
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", ]
map <- map[!duplicated(map$PROBEID), ]
write.csv(map, "F:/坏疽性脓皮病/outputs/mr/gut_skin/probe_to_gene.csv", row.names = FALSE)
cat("Mapped probes:", nrow(map), "of", length(probes), "\n")
