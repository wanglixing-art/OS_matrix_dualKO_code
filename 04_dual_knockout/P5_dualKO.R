# ============================================================
# P5: scTenifoldKnk 双 KO 对比（核心创新点）
#   KO-1: RUNX2  —— 基质区室（CAF/MSC + Osteoblastic，RUNX2-high）
#   KO-2: BUB1   —— 增殖区室（Proliferating，BUB1-high）
#
# 方法学同源: Yang JH et al., Molecules 2026;31:1901 (CKAP2 模板)
#   模板做法: 单一区室建 scGRN → 虚拟敲除 → DRG → 富集
#   本项目创新: 双区室双 KO 对照架构 —— DRG 重叠/特异性对比
#
# ===== 性能实测（本机 R 4.6.1）=====
#   pcNet 建网耗时 O(n_genes²)（内部 n_triplets = n*(n-1) 三元组构建）：
#     1,001 基因/200 细胞 = 132.5 s ；2,001 = 480.6 s
#   实测 15,109 基因（全基因）在 Step3 静默卡死（2.28 亿三元组）→ 必须降维
#   采用 N_GENES=1500 + N_NET=5 → 预计 5×2×~300s ≈ 50 min
#   Rcpp 后端已确认生效（verbose 输出 "Using compiled Rcpp backend"）
#
# ===== 【重要】本脚本与下游口径分工 =====
#   本脚本（P5_dualKO.R）职责 = **KO 计算引擎**：
#     读标注 Seurat → 分区室取子集 → HVG 降维 → scTenifoldKnk → 存 checkpoint
#     （data/processed/P5_knk_<tag>.rds，含 manifoldAlignment 可复算）
#   本脚本尾部的「双 KO 对比 / 绘图」段为**早期版本**，其 DRG 判定用 p.adj<0.05
#     —— 在 1,500 基因规模下仅 KO 靶基因自身能过，已废弃不用。
#   ✅ 最终 DRG 判定与对比以 **code/P5_final_drg.R** 为准：
#     DRG = (Z > 2) & (distance >= 1e-10)   [单向 Z + 排除浮点下溢伪影]
#   ✅ 富集与可视化以 **code/P5_post.R** 为准。
# ============================================================
suppressMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(scTenifoldKnk); library(ggplot2); library(patchwork)
})
set.seed(20260928)
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
DP <- file.path(ROOT, "data/processed")
LOG <- file.path(ROOT, "logs/P5_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}

HUB_RUNX2 <- "RUNX2"; HUB_BUB1 <- "BUB1"
COMP_STROMA <- c("CAF/MSC", "Osteoblastic")   # RUNX2-high 基质轴
COMP_PROLIF <- "Proliferating"                # BUB1-high 增殖轴
N_GENES   <- 1500     # 高变基因数（实测优化：平衡速度/覆盖）
N_CELLS   <- 1000     # 每区室子采样细胞上限
N_NET     <- 5        # 网络数量（>=5，scTenifoldNet 官方建议）
N_CELLSNET<- 300      # 每次网络子采样细胞数（默认 500，降至 300 提速）

# ---------- 1) 载入 & 区室提取 ----------
logmsg("Loading annotated Seurat object ...")
seu <- readRDS(file.path(DP, "P4_seurat_annotated.rds"))
logmsg("total cells:", ncol(seu))
ct <- as.character(seu$celltype_marker)
logmsg("celltype table:"); print(table(ct))

# ---------- 2) 通用 KO 执行函数 ----------
run_ko <- function(comp_cts, gKO, tag) {
  logmsg("===== KO:", gKO, "| compartment:", paste(comp_cts, collapse = "+"), "=====")
  idx <- which(ct %in% comp_cts)
  logmsg("  compartment cells:", length(idx))
  if (length(idx) > N_CELLS) { set.seed(1); idx <- sample(idx, N_CELLS) }
  sub <- subset(seu, cells = colnames(seu)[idx])
  logmsg("  sampled cells:", ncol(sub))

  cnt <- as.matrix(LayerData(sub, layer = "counts"))
  logmsg("  raw matrix:", nrow(cnt), "genes x", ncol(cnt), "cells")
  # 基因过滤：>=3 细胞检出
  cnt <- cnt[rowSums(cnt > 0) >= 3, , drop = FALSE]
  logmsg("  after expr filter:", nrow(cnt), "genes")
  stopifnot(gKO %in% rownames(cnt))

  # 高变基因降维（必须：全基因规模会导致 pcNet 三元组爆炸）
  libsz <- colSums(cnt)
  cpm <- t(t(cnt) / libsz) * 1e6
  v <- apply(log2(cpm + 1), 1, var)
  hvg <- rownames(cnt)[order(v, decreasing = TRUE)][seq_len(min(N_GENES, nrow(cnt)))]
  genes_keep <- unique(c(hvg, gKO, "RUNX2", "BUB1"))
  cnt <- cnt[genes_keep, , drop = FALSE]
  logmsg("  after HVG reduction:", nrow(cnt), "genes (incl. both hubs)")
  logmsg("  ", gKO, " mean raw count:", round(mean(cnt[gKO, ]), 3))
  gc()   # 保留 sub：QC 后若 hub 被剔除需回取原始计数

  logmsg("  running scTenifoldKnk (nNet=", N_NET, ", nCells/net=", N_CELLSNET, ") ...")
  t0 <- Sys.time()
  # 注意: scQC 内置 removeOutlierCells 在本版本对已过滤矩阵会报
  #   "x must be an array of at least two dimensions" → 自行做等价 QC 后传 qc=FALSE
  # 自行 QC: 库大小 >=1000；基因检出率 >=5%
  lib <- colSums(cnt)
  cnt <- cnt[, lib >= 1000, drop = FALSE]
  cnt <- cnt[rowMeans(cnt > 0) >= 0.05, , drop = FALSE]
  logmsg("  after QC:", nrow(cnt), "genes x", ncol(cnt), "cells")
  # 关键: 被 KO 的 hub 必须存在（区室内表达低会被 QC 剔除）；
  #   同时把另一个 hub 也强制加入（用于观察跨轴效应）
  cnt <- cnt[rowSums(cnt > 0) >= 3, , drop = FALSE]
  force_genes <- c(gKO, HUB_RUNX2, HUB_BUB1)
  add <- setdiff(force_genes, rownames(cnt))
  if (length(add) > 0) {
    raw_back <- as.matrix(LayerData(sub, layer = "counts"))[add, colnames(cnt), drop = FALSE]
    cnt <- rbind(cnt, raw_back)
    logmsg("  force-added hubs (dropped by QC):", paste(add, collapse = ", "))
  }
  stopifnot(gKO %in% rownames(cnt), "RUNX2" %in% rownames(cnt), "BUB1" %in% rownames(cnt))
  logmsg("  final matrix:", nrow(cnt), "genes x", ncol(cnt), "cells")
  rm(sub); gc()

  res <- scTenifoldKnk(
    countMatrix = cnt, gKO = gKO,
    qc = FALSE,
    nc_nNet = N_NET, nc_nCells = N_CELLSNET, nc_nComp = 3,
    nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9,
    td_K = 3, td_maxIter = 1000, td_maxError = 1e-5, td_nDecimal = 3,
    ma_nDim = 2, dr_empiricalNull = FALSE, nCores = 4)
  logmsg("  done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")

  dr <- as.data.table(res$diffRegulation)
  logmsg("  DRG cols:", paste(colnames(dr), collapse = ", "), "| rows:", nrow(dr))
  dr[, KO := gKO][, compartment := paste(comp_cts, collapse = "+")]
  fwrite(dr, file.path(TM, paste0("P5_DRG_", tag, ".csv")))
  pcol <- if ("p.adj" %in% colnames(dr)) "p.adj" else grep("^p", colnames(dr), value = TRUE)[1]
  logmsg("  using p column:", pcol)
  sig <- dr[get(pcol) < 0.05]
  logmsg("  significant DRG:", nrow(sig))
  logmsg("  top 15:")
  print(head(dr[order(get(pcol), -distance)], 15))
  saveRDS(res, file.path(DP, paste0("P5_knk_", tag, ".rds")))
  rm(cnt, res); gc()
  list(dr = dr, sig = sig, pcol = pcol)
}

r1 <- run_ko(COMP_STROMA, HUB_RUNX2, "RUNX2_stroma")
r2 <- run_ko(COMP_PROLIF, HUB_BUB1, "BUB1_prolif")

# ---------- 3) 双 KO 对比 ----------
dr1 <- r1$dr; dr2 <- r2$dr; s1 <- r1$sig; s2 <- r2$sig
logmsg("===== 双 KO 对比 =====")
logmsg("RUNX2-KO DRG:", nrow(dr1), "sig:", nrow(s1))
logmsg("BUB1-KO  DRG:", nrow(dr2), "sig:", nrow(s2))

ov    <- intersect(s1$gene, s2$gene)
only1 <- setdiff(s1$gene, s2$gene)
only2 <- setdiff(s2$gene, s1$gene)
logmsg("shared:", length(ov), "| RUNX2-only:", length(only1), "| BUB1-only:", length(only2))

cmp <- data.table(category = c("Shared DRG", "RUNX2-specific", "BUB1-specific"),
                  n = c(length(ov), length(only1), length(only2)))
fwrite(cmp, file.path(TM, "P5_DRG_overlap_summary.csv"))
fwrite(data.table(gene = ov),    file.path(TM, "P5_DRG_shared.csv"))
fwrite(data.table(gene = only1), file.path(TM, "P5_DRG_RUNX2_specific.csv"))
fwrite(data.table(gene = only2), file.path(TM, "P5_DRG_BUB1_specific.csv"))

allg <- rbind(dr1[, .(gene, KO, distance, Z, p.adj = get(r1$pcol))],
              dr2[, .(gene, KO, distance, Z, p.adj = get(r2$pcol))])
w <- dcast(allg, gene ~ KO, value.var = c("Z", "distance", "p.adj"))
fwrite(w, file.path(TM, "P5_DRG_wide_both_KO.csv"))

logmsg("cross-hit: BUB1 in RUNX2-KO DRG?", HUB_BUB1 %in% s1$gene,
       "| RUNX2 in BUB1-KO DRG?", HUB_RUNX2 %in% s2$gene)

# 富集分析（若 clusterProfiler 可用）
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  logmsg("Running GO enrichment on DRGs ...")
  suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  for (nm in c("shared", "RUNX2_specific", "BUB1_specific")) {
    gs <- switch(nm, shared = ov, RUNX2_specific = only1, BUB1_specific = only2)
    if (length(gs) < 10) { logmsg("  ", nm, ": too few genes (", length(gs), "), skip"); next }
    eg <- tryCatch(bitr(gs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                   error = function(e) NULL)
    if (is.null(eg) || nrow(eg) < 5) next
    ego <- tryCatch(enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
                    error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      fwrite(as.data.frame(ego), file.path(TM, paste0("P5_GO_", nm, ".csv")))
      logmsg("  ", nm, " GO terms:", nrow(as.data.frame(ego)))
    }
  }
} else {
  logmsg("clusterProfiler/org.Hs.eg.db 不可用，跳过富集")
}

# ---------- 4) 图 ----------
pal <- c("Shared DRG" = "#6A3D9A", "RUNX2-specific" = "#E31A1C", "BUB1-specific" = "#1F78B4")

pa <- ggplot(cmp, aes(category, n, fill = category)) +
  geom_col(width = 0.62, color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Dual knockout: DRG overlap", y = "DRGs (p.adj < 0.05)", x = NULL) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

zs <- merge(dr1[, .(gene, Z1 = Z, p1 = get(r1$pcol))],
            dr2[, .(gene, Z2 = Z, p2 = get(r2$pcol))], by = "gene")
zs[, cat := ifelse(gene %in% ov, "Shared",
                   ifelse(gene %in% only1, "RUNX2-specific",
                          ifelse(gene %in% only2, "BUB1-specific", "Other")))]
fwrite(zs, file.path(TM, "P5_DRG_Z_scatter.csv"))
zs[, is_ctrl := gene %in% c(HUB_RUNX2, HUB_BUB1)]
pb <- ggplot(zs, aes(Z1, Z2, color = cat)) +
  geom_point(size = 0.9, alpha = 0.7) +
  geom_point(data = zs[is_ctrl == TRUE], size = 3, shape = 21,
             fill = "yellow", color = "black", stroke = 0.6) +
  geom_text(data = zs[is_ctrl == TRUE], aes(label = gene), vjust = -1.1, size = 3.4, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  scale_color_manual(values = c("Shared" = "#6A3D9A", "RUNX2-specific" = "#E31A1C",
                                "BUB1-specific" = "#1F78B4", "Other" = "grey82"), name = "") +
  labs(title = "Dual KO perturbation landscape",
       subtitle = "Yellow = the two knockout targets themselves",
       x = "RUNX2-KO Z-score (stromal compartment)",
       y = "BUB1-KO Z-score (proliferative compartment)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(FG, "P5_dualKO_overview.png"), pa | pb,
       width = 13, height = 5.5, dpi = 300)
logmsg("saved P5_dualKO_overview.png")

# Top DRG 条形
mk_top <- function(dr, pcol, lab) {
  d <- dr[order(get(pcol), -abs(Z))]
  d <- d[seq_len(min(15, nrow(d)))]
  d[, .(gene, Z, KO = lab)]
}
tt <- rbind(mk_top(dr1, r1$pcol, "RUNX2-KO (stroma)"),
            mk_top(dr2, r2$pcol, "BUB1-KO (prolif)"))
tt[, gene := factor(gene, levels = unique(tt[order(abs(Z)), gene]))]
pc <- ggplot(tt, aes(gene, Z, fill = KO)) +
  geom_col(width = 0.7, color = "grey20", linewidth = 0.25) +
  coord_flip() + facet_wrap(~KO, scales = "free_y") +
  scale_fill_manual(values = c("RUNX2-KO (stroma)" = "#E31A1C",
                               "BUB1-KO (prolif)" = "#1F78B4"), guide = "none") +
  labs(title = "Top 15 perturbed genes per knockout", x = NULL, y = "Z-score") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FG, "P5_topDRG_per_KO.png"), pc, width = 11, height = 6.5, dpi = 300)
logmsg("saved P5_topDRG_per_KO.png")

logmsg("=== P5 SUMMARY ===")
cat("RUNX2-KO (stroma: CAF/MSC+Osteoblastic) DRG:", nrow(dr1), "sig:", nrow(s1), "\n")
cat("BUB1-KO  (prolif: Proliferating)          DRG:", nrow(dr2), "sig:", nrow(s2), "\n")
cat("shared:", length(ov), "| RUNX2-only:", length(only1), "| BUB1-only:", length(only2), "\n")
logmsg("P5 complete.")
