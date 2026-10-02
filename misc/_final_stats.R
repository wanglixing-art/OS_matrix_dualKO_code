suppressMessages(library(data.table))
TM <- "D:/projects/OS_matrix_dualKO/results/tables"
# 1) 超几何：DRG 是否富集于 P2 候选池（299/1501 背景）
for (k in c("RUNX2", "BUB1")) {
  d <- fread(file.path(TM, paste0("P5_DRG_", ifelse(k=="RUNX2","RUNX2_stroma","BUB1_prolif"), ".csv")))
  drg <- d[DRG == TRUE]$gene
  hit <- sum(drg %in% fread(file.path(TM,"P2_candidate_pool_v3_matrixaxis.csv"))$gene)
  q <- phyper(hit - 1, 299, 1501 - 299, length(drg), lower.tail = FALSE)
  cat(sprintf("%s-KO: %d DRG, %d in pool (expect %.1f), hypergeom P = %.4f\n",
              k, length(drg), hit, length(drg) * 299 / 1501, q))
}
# 2) BUB1-KO 全部 17 个 DRG
d2 <- fread(file.path(TM, "P5_DRG_BUB1_prolif.csv"))
cat("\nBUB1-KO DRG (Z>2), all 17:\n")
print(d2[DRG == TRUE, .(gene, Z = round(Z, 2), FC = round(FC, 3))])
# 3) RUNX2 的池内 4 个基因
d1 <- fread(file.path(TM, "P5_DRG_RUNX2_stroma.csv"))
pool <- fread(file.path(TM, "P2_candidate_pool_v3_matrixaxis.csv"))$gene
cat("\nRUNX2-KO DRG in P2 pool:", paste(intersect(d1[DRG==TRUE]$gene, pool), collapse=", "), "\n")
cat("BUB1-KO  DRG in P2 pool:", paste(intersect(d2[DRG==TRUE]$gene, pool), collapse=", "), "\n")
