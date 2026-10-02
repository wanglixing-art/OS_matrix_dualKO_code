# ============================================================
# P5-final: DRG 判定（Z 阈值法，scTenifoldKnk 文献标准口径）
#
# 口径依据（重要）：
#   1) scTenifoldKnk 的 p 值来自 pchisq(FC, df=1) 或 Efron 经验零分布；
#      1,500 基因规模下 BH 校正后仅 KO 靶基因自身能过 FDR<0.05
#      —— 这是该方法的已知特性，文献中普遍如此（非本数据缺陷）。
#   2) 原论文 (Osorio et al., NAR 2022) 及模板文献 (Yang 2026 Molecules)
#      的 DRG 提取均以 **Z 值 / distance 排序** 为准，再做下游富集。
#   → 主口径：Z > 2 为显著差异调控基因 (DRG)；同时报告 FDR 结果作对照。
#
# 复用已保存的 manifoldAlignment（P5_knk_*.rds），无需重跑建网。
# ============================================================
suppressMessages({ library(scTenifoldKnk); library(data.table) })
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); DP <- file.path(ROOT, "data/processed")
LOG <- file.path(ROOT, "logs/P5_final.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
Z_THR <- 2      # 正向阈值（见下）
# 重要：只用正向 Z > 2。
#   实测发现若干基因（ISCU/FXR1/GSTM3 等）distance 下溢到 ~1e-17（浮点下界），
#   经验零分布(mean=0.024,sd=1.028) 会把这类「几乎无位移」判成负向显著
#   (Z=-2.59, FC~1e-22)，属数值伪影而非真实扰动，故负向一侧不纳入 DRG。
DIST_MIN <- 1e-10   # distance 低于此值视为下溢伪影

res_list <- list()
for (tag in c("RUNX2_stroma", "BUB1_prolif")) {
  f <- file.path(DP, paste0("P5_knk_", tag, ".rds"))
  if (!file.exists(f)) { logmsg("MISSING (KO 尚未跑完):", tag); next }
  logmsg("=====", tag, "=====")
  r <- readRDS(f)
  d <- as.data.table(dRegulation(r$manifoldAlignment, empiricalNull = TRUE))
  d[, KO := sub("_.*", "", tag)]
  d[, compartment := ifelse(KO == "RUNX2", "CAF/MSC+Osteoblastic", "Proliferating")]
  d[, underflow := distance < DIST_MIN]
  d[, DRG := (Z > Z_THR) & (!underflow)]
  logmsg("  genes:", nrow(d), "| DRG (Z>", Z_THR, "):", sum(d$DRG),
         "| FDR<0.05:", sum(d$p.adj < 0.05),
         "| underflow-excluded:", sum(d$underflow))
  logmsg("  top 20 by Z:")
  print(head(d[order(-Z), .(gene, distance, Z, FC, p.value, p.adj, DRG)], 20))
  fwrite(d, file.path(TM, paste0("P5_DRG_", tag, ".csv")))
  res_list[[tag]] <- d
}

if (length(res_list) == 2) {
  d1 <- res_list[["RUNX2_stroma"]]; d2 <- res_list[["BUB1_prolif"]]
  s1 <- d1[DRG == TRUE]$gene; s2 <- d2[DRG == TRUE]$gene
  ov <- intersect(s1, s2); only1 <- setdiff(s1, s2); only2 <- setdiff(s2, s1)
  logmsg("===== 双 KO 对比 (Z>", Z_THR, ") =====")
  logmsg("RUNX2-KO DRG:", length(s1), "| BUB1-KO DRG:", length(s2))
  logmsg("shared:", length(ov), "| RUNX2-only:", length(only1), "| BUB1-only:", length(only2))

  cmp <- data.table(category = c("Shared", "RUNX2-specific", "BUB1-specific"),
                    n = c(length(ov), length(only1), length(only2)))
  fwrite(cmp, file.path(TM, "P5_DRG_overlap_summary.csv"))
  fwrite(data.table(gene = ov),    file.path(TM, "P5_DRG_shared.csv"))
  fwrite(data.table(gene = only1), file.path(TM, "P5_DRG_RUNX2_specific.csv"))
  fwrite(data.table(gene = only2), file.path(TM, "P5_DRG_BUB1_specific.csv"))

  w <- merge(d1[, .(gene, Z_RUNX2 = Z, p_RUNX2 = p.value, DRG_RUNX2 = DRG)],
             d2[, .(gene, Z_BUB1 = Z, p_BUB1 = p.value, DRG_BUB1 = DRG)],
             by = "gene", all = TRUE)
  fwrite(w, file.path(TM, "P5_DRG_wide_both_KO.csv"))

  logmsg("shared DRG (Z>2 in both):"); print(ov)
  logmsg("RUNX2-specific (top 20 by Z):")
  print(head(d1[gene %in% only1][order(-Z), .(gene, Z)], 20))
  logmsg("BUB1-specific (top 20 by Z):")
  print(head(d2[gene %in% only2][order(-Z), .(gene, Z)], 20))
}

logmsg("P5-final complete.")
