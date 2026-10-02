# ============================================================
# P4g-score: inferCNV 结果的权威打分（独立于 P4g_infercnv.R）
#
# 背景修正：inferCNV 1.28.0 在 HMM 模式下产出的
#   infercnv.*.observations.txt  = HMM 离散状态（1–6 整数），不可用于连续打分
#   expr.infercnv.preliminary.dat = 连续残差（log2FC 校正后，围绕 1.0 波动）← 正确数据源
#   17_HMM_pred*.pred_cnv_genes.dat = 长表（cell_group_name | gene_region | state）← 二级证据
#
# 恶性判定口径（双证据）：
#   A) 连续 CNV 分数 = 每细胞 |残差 - 1| 的均值；阈值 = 参照细胞 95 分位
#   B) HMM 状态偏离度 = 每细胞 HMM 状态偏离中性(3) 的程度，来自 pred_cnv_genes
# ============================================================
suppressMessages({ library(data.table); library(ggplot2) })
set.seed(20260928)
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
INFCNV <- file.path(ROOT, "results/infercnv")
LOG <- file.path(ROOT, "logs/P4g_score.log")
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)

# ---------- 0) 细胞元数据（group/celltype/sample）----------
ann <- fread(file.path(TM, "P4g_infercnv_annotation.csv"))  # cell, group
logmsg("annotation cells:", nrow(ann), "| groups:")
print(table(ann$group))

# 从 Seurat 对象取 celltype / sample / 代理法恶性标签 映射
suppressMessages(library(Seurat))
seu <- readRDS(file.path(ROOT, "data/processed/P4_seurat_annotated.rds"))
md <- data.table(cell = colnames(seu),
                 celltype = as.character(seu$celltype_marker),
                 sample = as.character(seu$sample),
                 malignant_proxy = as.character(seu$malignant_proxy))
rm(seu); gc()
logmsg("metadata cells:", nrow(md), "| proxy malignant:", sum(md$malignant_proxy == "Malignant", na.rm = TRUE))

# ---------- 1) 连续残差矩阵 -> 每细胞 CNV 分数 ----------
f <- file.path(INFCNV, "expr.infercnv.preliminary.dat")
logmsg("reading residual matrix:", f)
mat <- fread(f, sep = "\t", header = TRUE, data.table = FALSE)
genes <- mat[[1]]
mat <- as.matrix(mat[, -1, drop = FALSE])
rownames(mat) <- genes
logmsg("matrix:", nrow(mat), "genes x", ncol(mat), "cells")

# 每细胞 CNV 分数
cnv_score <- colMeans(abs(mat - 1), na.rm = TRUE)
# 每细胞变异幅度（另一种口径：sd 而非 mean abs）
cnv_sd <- apply(mat, 2, sd, na.rm = TRUE)
rm(mat); gc()

ann_map <- setNames(ann$group, ann$cell)
grp <- ann_map[names(cnv_score)]
# 不在 ann 里的细胞视为 ref（不应发生）
grp[is.na(grp)] <- "ref"

ct_map <- setNames(md$celltype, md$cell)
sp_map <- setNames(md$sample, md$cell)

res <- data.table(
  cell      = names(cnv_score),
  cnv_score = as.numeric(cnv_score),
  cnv_sd    = as.numeric(cnv_sd),
  group     = grp,
  celltype  = unname(ct_map[names(cnv_score)]),
  sample    = unname(sp_map[names(cnv_score)])
)

ref_score <- res[group == "ref", cnv_score]
thr95 <- as.numeric(quantile(ref_score, 0.95, na.rm = TRUE))
thr99 <- as.numeric(quantile(ref_score, 0.99, na.rm = TRUE))
logmsg("ref CNV score: median=", round(median(ref_score), 4),
       " mean=", round(mean(ref_score), 4),
       " sd=", round(sd(ref_score), 4))
logmsg("thresholds: 95th=", round(thr95, 4), " 99th=", round(thr99, 4))

res[, malignant_95 := ifelse(cnv_score > thr95, "Malignant", "Non-malignant")]
res[, malignant_99 := ifelse(cnv_score > thr99, "Malignant", "Non-malignant")]
# z-score 口径（相对参照分布）
mu <- mean(ref_score); sg <- sd(ref_score)
res[, cnv_z := (cnv_score - mu) / sg]
res[, malignant_z2 := ifelse(cnv_z > 2, "Malignant", "Non-malignant")]

fwrite(res, file.path(TM, "P4g_percell_cnv_score.csv"))
logmsg("malignant (95th pct):", sum(res$malignant_95 == "Malignant"), "/", nrow(res))
logmsg("malignant (99th pct):", sum(res$malignant_99 == "Malignant"))
logmsg("malignant (z>2)     :", sum(res$malignant_z2 == "Malignant"))

# ---------- 2) 按细胞类型 / 样本汇总 ----------
ct_sum <- res[, .(n = .N,
                  mal95 = sum(malignant_95 == "Malignant"),
                  pct95 = round(100 * mean(malignant_95 == "Malignant"), 1),
                  mean_score = round(mean(cnv_score), 4),
                  mean_z = round(mean(cnv_z), 3)),
              by = celltype][order(-pct95)]
fwrite(ct_sum, file.path(TM, "P4g_malignant_pct_by_celltype.csv"))
logmsg("by celltype:")
print(ct_sum)

grp_sum <- res[, .(n = .N, mal95 = sum(malignant_95 == "Malignant"),
                   pct95 = round(100 * mean(malignant_95 == "Malignant"), 1),
                   mean_score = round(mean(cnv_score), 4)), by = group]
fwrite(grp_sum, file.path(TM, "P4g_malignant_pct_by_group.csv"))
logmsg("by group (ref vs obs):")
print(grp_sum)

sp_sum <- res[, .(n = .N, mal95 = sum(malignant_95 == "Malignant"),
                  pct95 = round(100 * mean(malignant_95 == "Malignant"), 1),
                  mean_score = round(mean(cnv_score), 4)), by = sample][order(sample)]
fwrite(sp_sum, file.path(TM, "P4g_malignant_pct_by_sample.csv"))
logmsg("by sample:")
print(sp_sum)

# obs 内按细胞类型
logmsg("obs-only by celltype:")
print(res[group == "obs", .(n = .N, mal95 = sum(malignant_95 == "Malignant"),
                            pct95 = round(100*mean(malignant_95 == "Malignant"), 1),
                            mean_score = round(mean(cnv_score), 4)), by = celltype])

# ---------- 3) HMM 状态证据（pred_cnv_genes 长表）----------
pf <- file.path(INFCNV, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat")
if (file.exists(pf)) {
  pg <- fread(pf)
  logmsg("pred_cnv_genes rows:", nrow(pg), "| states:", paste(sort(unique(pg$state)), collapse = ","))
  # 每个亚群的（按基因数加权的）状态均值
  cg <- pg[, .(mean_state = mean(state),
               frac_del = mean(state < 3),     # 拷贝丢失
               frac_amp = mean(state > 3),     # 拷贝获得
               n_genes = .N), by = cell_group_name]
  fwrite(cg, file.path(TM, "P4g_hmm_state_by_subcluster.csv"))
  logmsg("HMM state by subcluster (top20 by |state-3|):")
  cg[, dev := abs(mean_state - 3)]
  print(head(cg[order(-dev)], 20))

  # 亚群 -> 细胞 映射
  cgf <- file.path(INFCNV, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings")
  if (file.exists(cgf)) {
    cgmap <- fread(cgf)
    setnames(cgmap, c("cell_group_name", "cell"))
    m <- merge(res, cgmap[, .(cell, cell_group_name)], by = "cell", all.x = TRUE)
    m <- merge(m, cg[, .(cell_group_name, mean_state, frac_del, frac_amp)], by = "cell_group_name", all.x = TRUE)
    m[, hmm_dev := abs(mean_state - 3)]
    fwrite(m, file.path(TM, "P4g_percell_cnv_with_hmm.csv"))

    # HMM 证据按细胞类型
    hmm_ct <- m[!is.na(mean_state), .(n = .N, mean_state = round(mean(mean_state), 3),
                                      frac_amp = round(mean(frac_amp), 3),
                                      frac_del = round(mean(frac_del), 3)),
                by = celltype][order(-mean_state)]
    fwrite(hmm_ct, file.path(TM, "P4g_hmm_state_by_celltype.csv"))
    logmsg("HMM mean state by celltype (obs should be >3 = gain):")
    print(hmm_ct)

    # 两种证据的一致性
    mm <- m[!is.na(mean_state)]
    ct2 <- mm[, .N, by = .(cont = malignant_95, hmm = ifelse(mean_state > 3.2, "Amp", ifelse(mean_state < 2.8, "Del", "Neutral")))]
    fwrite(ct2, file.path(TM, "P4g_evidence_concordance.csv"))
    logmsg("concordance (continuous vs HMM):")
    print(dcast(ct2, cont ~ hmm, value.var = "N", fill = 0))
  }
}

# ---------- 4) 与代理法交叉表 ----------
if ("malignant_proxy" %in% colnames(md)) {
  mg <- merge(res, md[, .(cell, malignant_proxy)], by = "cell", all.x = TRUE)
  ct <- mg[!is.na(malignant_proxy), .N, by = .(inferCNV = malignant_95, proxy = malignant_proxy)]
  fwrite(ct, file.path(TM, "P4g_infercnv_vs_proxy_crosstab.csv"))
  logmsg("inferCNV vs proxy crosstab:")
  print(dcast(ct, inferCNV ~ proxy, value.var = "N", fill = 0))
}

# ---------- 5) 图 ----------
pal_ct <- c("Osteoblastic" = "#D7263D", "Proliferating" = "#F46036",
            "CAF/MSC" = "#2E86AB", "Myeloid" = "#6A994E", "Tcell" = "#BC4B51",
            "Endothelial" = "#8D99AE", "Bcell" = "#7209B7", "NK" = "#F2A900",
            "Osteoclast" = "#4A4E69")

ctp <- ct_sum[celltype %in% names(pal_ct)]
p1 <- ggplot(ctp, aes(x = reorder(celltype, pct95), y = pct95, fill = celltype)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = paste0(pct95, "% (n=", n, ")")), hjust = -0.06, size = 3.2) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = pal_ct, guide = "none") +
  scale_y_continuous(limits = c(0, 125), expand = c(0, 0)) +
  labs(title = "inferCNV malignant fraction by cell type",
       subtitle = paste0("threshold = 95th pct of reference CNV score (", round(thr95, 3), ")"),
       x = NULL, y = "% Malignant") +
  theme_bw(base_size = 12) +
  theme(plot.margin = margin(6, 52, 6, 6), panel.grid.minor = element_blank())

p2 <- ggplot(res, aes(x = reorder(celltype, cnv_score, median), y = cnv_score, fill = group)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  geom_hline(yintercept = thr95, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("obs" = "#D7263D", "ref" = "#2E86AB"), name = "") +
  labs(title = "Per-cell CNV score: obs (osteoblastic lineage) vs ref",
       subtitle = "red dashed = malignant threshold", x = NULL, y = "CNV score") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.minor = element_blank())

p3 <- ggplot(sp_sum, aes(sample, pct95)) +
  geom_col(fill = "#2E86AB", width = 0.68, color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = paste0(pct95, "%\nn=", n)), vjust = -0.35, size = 3.2) +
  scale_y_continuous(limits = c(0, max(sp_sum$pct95) * 1.28), expand = c(0, 0)) +
  labs(title = "Malignant fraction per sample", x = NULL, y = "% Malignant") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

p4 <- ggplot(grp_sum, aes(group, mean_score, fill = group)) +
  geom_col(width = 0.5, color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = round(mean_score, 3)), vjust = -0.4, size = 3.6) +
  scale_fill_manual(values = c("obs" = "#D7263D", "ref" = "#2E86AB"), guide = "none") +
  labs(title = "Mean CNV score: obs vs ref", x = NULL, y = "Mean CNV score") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

library(patchwork)
gg <- (p1 | p2) / (p3 | p4)
ggsave(file.path(FG, "P4g_infercnv_malignancy_summary.png"), gg,
       width = 14, height = 9.5, dpi = 300)
logmsg("saved P4g_infercnv_malignancy_summary.png")

# ---------- 6) 真实 CNV 热图（连续残差，按细胞类型排序）----------
of <- file.path(INFCNV, "expr.infercnv.preliminary.dat")
logmsg("building heatmap from residual matrix ...")
hm <- fread(of, sep = "\t", header = TRUE, data.table = FALSE)
g <- hm[[1]]; hm <- as.matrix(hm[, -1, drop = FALSE]); rownames(hm) <- g
# 只取 obs 细胞 + 按细胞类型抽样，控制尺寸
obs_cells <- res[group == "obs", cell]
obs_cells <- intersect(obs_cells, colnames(hm))
sel <- res[cell %in% obs_cells][order(celltype, -cnv_score)]
set.seed(1)
sel <- sel[, .SD[sample(.N, min(.N, 500))], by = celltype]
cols <- sel$cell
sub <- hm[, cols, drop = FALSE]
ord <- c("Osteoblastic", "Proliferating")
cell_ct_h <- sel$celltype[match(cols, sel$cell)]
o <- order(match(cell_ct_h, ord), -sel$cnv_score[match(cols, sel$cell)])
sub <- sub[, o, drop = FALSE]; cell_ct_h <- cell_ct_h[o]

df <- data.table(gene = rep(rownames(sub), times = ncol(sub)),
                 cell = rep(colnames(sub), each = nrow(sub)),
                 val = as.vector(sub))
df[, gidx := rep(seq_len(nrow(sub)), times = ncol(sub))]
vc <- quantile(abs(df$val - 1), 0.99, na.rm = TRUE)
phm <- ggplot(df, aes(cell, gidx, fill = val)) +
  geom_raster() +
  scale_fill_gradient2(low = "#00008B", mid = "white", high = "#8B0000",
                       midpoint = 1, limits = c(1 - vc, 1 + vc),
                       oob = scales::squish, name = "CNV\nresidual") +
  labs(title = "inferCNV heatmap — osteoblastic lineage (obs cells)",
       subtitle = "blue = copy loss, red = copy gain; genes in chromosome order",
       x = "Cells (grouped by cell type, sorted by CNV score)", y = "Genes") +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
ggsave(file.path(FG, "P4g_infercnv_heatmap_obs.png"), phm,
       width = 14, height = 6, dpi = 300)
logmsg("saved P4g_infercnv_heatmap_obs.png")

# ---------- 7) hub 表达 vs inferCNV 恶性状态 ----------
hb <- FetchData(readRDS(file.path(ROOT, "data/processed/P4_seurat_annotated.rds")),
                vars = c("BUB1", "RUNX2"), layer = "data")
hb$cell <- rownames(hb)
hbt <- as.data.table(hb)
mg2 <- merge(hbt, res[, .(cell, malignant_95, cnv_score, celltype, group)], by = "cell")
mg2 <- mg2[!is.na(malignant_95)]
long <- melt(mg2, id.vars = c("cell", "malignant_95", "group", "celltype"),
             measure.vars = c("BUB1", "RUNX2"), variable.name = "hub", value.name = "expr")
fwrite(mg2, file.path(TM, "P4g_hub_expr_vs_infercnv_malignant.csv"))
ph <- ggplot(long[group == "obs"], aes(malignant_95, expr, fill = malignant_95)) +
  geom_violin(scale = "width", alpha = 0.9, linewidth = 0.3, color = "grey25") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", linewidth = 0.3) +
  facet_wrap(~hub, scales = "free_y") +
  scale_fill_manual(values = c("Malignant" = "#D7263D", "Non-malignant" = "#2E86AB"), name = "") +
  labs(title = "Hub expression vs inferCNV malignancy (obs cells)",
       x = NULL, y = "log-normalized expression") +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(FG, "P4g_hub_vs_infercnv_malignant.png"), ph,
       width = 8, height = 4.8, dpi = 300)
logmsg("saved P4g_hub_vs_infercnv_malignant.png")

logmsg("=== P4g-score SUMMARY ===")
cat("cells:", nrow(res), "| ref:", sum(res$group == "ref"), "| obs:", sum(res$group == "obs"), "\n")
cat("threshold(95pct):", round(thr95, 4), "\n")
cat("malignant:", sum(res$malignant_95 == "Malignant"), "\n")
cat("obs malignant pct:", round(100 * mean(res[group == "obs", malignant_95] == "Malignant"), 1), "%\n")
cat("ref malignant pct:", round(100 * mean(res[group == "ref", malignant_95] == "Malignant"), 1), "%\n")
logmsg("P4g-score complete.")
