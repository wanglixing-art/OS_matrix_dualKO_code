# ============================================================
# P4g-score2: inferCNV 打分（v2，顺序对齐修复版）
#
# v1 失败根因: dat 列名与 ann$cell 的名字匹配链断裂
#   （fread header 首列吞掉第一个细胞名 → colMeans 名字错位 → 全 NA → 兜底全 ref → 阈值被污染）
# v2 修复: inferCNV 保持输入顺序输出，dat 数据列 j ↔ ann 行 j，纯位置对齐 + 断言
# ============================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(20260928)
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
INFCNV <- file.path(ROOT, "results/infercnv")
LOG <- file.path(ROOT, "logs/P4g_score2.log")
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)

# ---------- 0) annotation（inferCNV 输入顺序，txt 版：无表头制表符）----------
ann_txt <- fread(file.path(TM, "P4g_infercnv_annotation.txt"), header = FALSE,
                 col.names = c("cell", "group"))
logmsg("annotation rows:", nrow(ann_txt), "| group counts:")
print(table(ann_txt$group))
stopifnot(nrow(ann_txt) == 10612)

# Seurat 元数据（celltype/sample/proxy）
suppressMessages(library(Seurat))
seu <- readRDS(file.path(ROOT, "data/processed/P4_seurat_annotated.rds"))
md <- data.table(cell = colnames(seu),
                 celltype = as.character(seu$celltype_marker),
                 sample   = as.character(seu$sample),
                 malignant_proxy = as.character(seu$malignant_proxy))
rm(seu); gc()
logmsg("metadata cells:", nrow(md), "| proxy malignant:", sum(md$malignant_proxy == "Malignant", na.rm = TRUE))

# ann$cell 与 md$cell 必须能匹配（同一来源：colnames(sub)）
mfrac <- mean(ann_txt$cell %in% md$cell)
logmsg("ann->md cell match fraction:", round(mfrac, 4))
stopifnot(mfrac > 0.99)
ann <- merge(ann_txt, md, by = "cell", all.x = TRUE)
stopifnot(nrow(ann) == nrow(ann_txt))
# 恢复 inferCNV 输入顺序（merge 会重排）
ann <- ann_txt[, .(cell, group)][ann, on = "cell"]
logmsg("ann merged | head:")
print(head(ann, 3))

# ---------- 1) 读连续残差矩阵，位置对齐 ----------
# v1 教训: fread(header=TRUE) 因 header 行(10,612 token)与数据行(10,613 字段)不对齐，
#   把 SAMD11 数据行误当 header → 全部错位。改用确定性解析：
#   header 行 = 10,612 个细胞名（对应数据行的第 2..10,613 列）；skip=1 读纯数据。
f <- file.path(INFCNV, "expr.infercnv.preliminary.dat")
logmsg("reading residual matrix (deterministic parse) ...")
hdr_line <- readLines(f, n = 1)
cells <- strsplit(hdr_line, "\t", fixed = TRUE)[[1]]
logmsg("header cells:", length(cells))
stopifnot(length(cells) == nrow(ann_txt))
dt <- fread(f, sep = "\t", header = FALSE, skip = 1, showProgress = FALSE)
stopifnot(ncol(dt) == length(cells) + 1)   # 基因名列 + 值列
genes <- as.character(dt[[1]])
mat <- as.matrix(dt[, -1, drop = FALSE])
colnames(mat) <- cells
rownames(mat) <- genes
rm(dt); gc()
logmsg("matrix:", nrow(mat), "genes x", ncol(mat), "cells")
# 顺序断言: dat 列顺序 == annotation 行顺序（inferCNV 保持输入顺序）
stopifnot(all(colnames(mat) == ann_txt$cell))

cnv_score <- colMeans(abs(mat - 1), na.rm = TRUE)
cnv_sd    <- apply(mat, 2, sd, na.rm = TRUE)
rm(mat); gc()

res <- data.table(
  cell      = ann$cell,           # 位置对齐（列 j ↔ 行 j）
  cnv_score = as.numeric(cnv_score),
  cnv_sd    = as.numeric(cnv_sd),
  group     = ann$group,
  celltype  = ann$celltype,
  sample    = ann$sample,
  malignant_proxy = ann$malignant_proxy
)
logmsg("res built | NA check: celltype=", sum(is.na(res$celltype)),
       " sample=", sum(is.na(res$sample)), " group=", sum(is.na(res$group)))
stopifnot(!any(is.na(res$group)), !any(is.na(res$celltype)), !any(is.na(res$sample)))

# ---------- 2) 阈值只基于真参照 ----------
ref_score <- res[group == "ref", cnv_score]
stopifnot(length(ref_score) == sum(res$group == "ref"))
thr95 <- as.numeric(quantile(ref_score, 0.95, na.rm = TRUE))
thr99 <- as.numeric(quantile(ref_score, 0.99, na.rm = TRUE))
mu <- mean(ref_score); sg <- sd(ref_score)
res[, cnv_z := (cnv_score - mu) / sg]
res[, malignant_95 := ifelse(cnv_score > thr95, "Malignant", "Non-malignant")]
res[, malignant_99 := ifelse(cnv_score > thr99, "Malignant", "Non-malignant")]
res[, malignant_z2 := ifelse(cnv_z > 2, "Malignant", "Non-malignant")]
logmsg("ref n=", length(ref_score), " median=", round(median(ref_score), 4),
       " | thr95=", round(thr95, 4), " thr99=", round(thr99, 4))
logmsg("malignant: 95th=", sum(res$malignant_95 == "Malignant"),
       " 99th=", sum(res$malignant_99 == "Malignant"),
       " z>2=", sum(res$malignant_z2 == "Malignant"))

fwrite(res, file.path(TM, "P4g_percell_cnv_score.csv"))

# ---------- 3) 汇总表 ----------
ct_sum <- res[, .(n = .N, mal95 = sum(malignant_95 == "Malignant"),
                  pct95 = round(100 * mean(malignant_95 == "Malignant"), 1),
                  mean_score = round(mean(cnv_score), 4),
                  mean_z = round(mean(cnv_z), 3)), by = celltype][order(-pct95)]
fwrite(ct_sum, file.path(TM, "P4g_malignant_pct_by_celltype.csv"))
logmsg("by celltype:"); print(ct_sum)

grp_sum <- res[, .(n = .N, mal95 = sum(malignant_95 == "Malignant"),
                   pct95 = round(100 * mean(malignant_95 == "Malignant"), 1),
                   mean_score = round(mean(cnv_score), 4)), by = group]
fwrite(grp_sum, file.path(TM, "P4g_malignant_pct_by_group.csv"))
logmsg("by group:"); print(grp_sum)

sp_sum <- res[, .(n = .N, mal95 = sum(malignant_95 == "Malignant"),
                  pct95 = round(100 * mean(malignant_95 == "Malignant"), 1),
                  mean_score = round(mean(cnv_score), 4)), by = sample][order(sample)]
fwrite(sp_sum, file.path(TM, "P4g_malignant_pct_by_sample.csv"))
logmsg("by sample:"); print(sp_sum)

logmsg("obs-only by celltype:")
print(res[group == "obs", .(n = .N, mal95 = sum(malignant_95 == "Malignant"),
                            pct95 = round(100 * mean(malignant_95 == "Malignant"), 1),
                            mean_score = round(mean(cnv_score), 4)), by = celltype])

# ---------- 4) 与代理法交叉表（同 cell 集合内）----------
ct_x <- res[group == "obs", .N, by = .(inferCNV = malignant_95, proxy = malignant_proxy)]
fwrite(ct_x, file.path(TM, "P4g_infercnv_vs_proxy_crosstab.csv"))
logmsg("inferCNV vs proxy crosstab (obs cells):")
print(dcast(ct_x, inferCNV ~ proxy, value.var = "N", fill = 0))

# ---------- 5) HMM 状态证据 ----------
pf <- file.path(INFCNV, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat")
cgf <- file.path(INFCNV, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings")
if (file.exists(pf) && file.exists(cgf)) {
  pg <- fread(pf)
  logmsg("pred_cnv_genes rows:", nrow(pg), "| states:", paste(sort(unique(pg$state)), collapse = ","))
  cg <- pg[, .(mean_state = mean(state), frac_del = mean(state < 3),
               frac_amp = mean(state > 3), n_genes = .N), by = cell_group_name]
  fwrite(cg, file.path(TM, "P4g_hmm_state_by_subcluster.csv"))

  cgmap <- fread(cgf, header = TRUE)
  stopifnot(all(c("cell_group_name", "cell") %in% colnames(cgmap)))
  m <- merge(res, cgmap[, .(cell, cell_group_name)], by = "cell", all.x = TRUE)
  m <- merge(m, cg[, .(cell_group_name, mean_state, frac_del, frac_amp)],
             by = "cell_group_name", all.x = TRUE)
  m[, hmm_dev := abs(mean_state - 3)]
  fwrite(m, file.path(TM, "P4g_percell_cnv_with_hmm.csv"))

  hmm_ct <- m[!is.na(mean_state), .(n = .N, mean_state = round(mean(mean_state), 3),
                                    frac_amp = round(mean(frac_amp), 3),
                                    frac_del = round(mean(frac_del), 3)), by = celltype][order(-mean_state)]
  fwrite(hmm_ct, file.path(TM, "P4g_hmm_state_by_celltype.csv"))
  logmsg("HMM mean state by celltype:"); print(hmm_ct)

  mm <- m[!is.na(mean_state)]
  ct2 <- mm[, .N, by = .(cont = malignant_95,
                         hmm = ifelse(mean_state > 3.2, "Amp", ifelse(mean_state < 2.8, "Del", "Neutral")))]
  fwrite(ct2, file.path(TM, "P4g_evidence_concordance.csv"))
  logmsg("concordance (continuous vs HMM):")
  print(dcast(ct2, cont ~ hmm, value.var = "N", fill = 0))
}

# ---------- 6) 主汇总图 ----------
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

gg <- (p1 | p2) / (p3 | p4)
ggsave(file.path(FG, "P4g_infercnv_malignancy_summary.png"), gg,
       width = 14, height = 9.5, dpi = 300)
logmsg("saved P4g_infercnv_malignancy_summary.png")

# ---------- 7) CNV 热图（复用内存中已有的 res，不再重读大矩阵）----------
logmsg("building heatmap from residual matrix ...")
hdr_line2 <- readLines(f, n = 1)
cells2 <- strsplit(hdr_line2, "\t", fixed = TRUE)[[1]]
dt2 <- fread(f, sep = "\t", header = FALSE, skip = 1, showProgress = FALSE)
hm <- as.matrix(dt2[, -1, drop = FALSE])
gname <- as.character(dt2[[1]])
colnames(hm) <- cells2
rownames(hm) <- gname
rm(dt2); gc()
stopifnot(all(colnames(hm) == ann_txt$cell))

# 只画 obs 细胞（列 1..6000 段内 ann$group=="obs" 的列），按 celltype 分组、组内按分数降序
obs_j <- which(ann$group == "obs")
set.seed(1)
sel_j <- tapply(obs_j, ann$celltype[obs_j], function(jj)
  jj[order(-cnv_score[jj])][seq_len(min(length(jj), 500))])
sel_j <- unlist(sel_j, use.names = FALSE)
sub <- hm[, sel_j, drop = FALSE]
cell_ct_h <- ann$celltype[sel_j]
# 列序：Osteoblastic 在前，Proliferating 在后（组内已按分数排）
o <- order(match(cell_ct_h, c("Osteoblastic", "Proliferating")))
sub <- sub[, o, drop = FALSE]; cell_ct_h <- cell_ct_h[o]
logmsg("heatmap matrix:", nrow(sub), "genes x", ncol(sub), "obs cells")

df <- as.data.table(as.table(sub))
setnames(df, c("gene", "cell", "val"))
df[, gene := as.character(gene)][, cell := as.character(cell)]
df[, gidx := match(gene, rownames(sub))]
# matrix melt 后 cell 是列名因子，顺序即 sub 列序
vc <- quantile(abs(df$val - 1), 0.99, na.rm = TRUE)
phm <- ggplot(df, aes(cell, gidx, fill = val)) +
  geom_raster() +
  scale_fill_gradient2(low = "#00008B", mid = "white", high = "#8B0000",
                       midpoint = 1, limits = c(1 - vc, 1 + vc),
                       oob = scales::squish, name = "CNV\nresidual") +
  labs(title = "inferCNV heatmap: osteoblastic lineage (obs) cells",
       subtitle = "Blue = copy loss, red = copy gain; genes in chromosome order; cells grouped by cell type, sorted by CNV score",
       x = "Obs cells (Osteoblastic | Proliferating)", y = "Genes") +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(FG, "P4g_infercnv_heatmap_obs.png"), phm,
       width = 14, height = 6, dpi = 300)
logmsg("saved P4g_infercnv_heatmap_obs.png")
rm(hm, sub, df); gc()

# ---------- 8) hub 表达 vs inferCNV 恶性 ----------
suppressMessages(library(Seurat))
seu <- readRDS(file.path(ROOT, "data/processed/P4_seurat_annotated.rds"))
hb <- FetchData(seu, vars = c("BUB1", "RUNX2"), layer = "data")
hb$cell <- rownames(hb)
hbt <- as.data.table(hb)
mg2 <- merge(hbt, res[, .(cell, malignant_95, cnv_score, celltype, group)], by = "cell")
logmsg("hub-expr merged rows:", nrow(mg2))
fwrite(mg2, file.path(TM, "P4g_hub_expr_vs_infercnv_malignant.csv"))
long <- melt(mg2, id.vars = c("cell", "malignant_95", "group", "celltype"),
             measure.vars = c("BUB1", "RUNX2"), variable.name = "hub", value.name = "expr")
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

logmsg("=== P4g-score2 SUMMARY ===")
cat("cells:", nrow(res), "| ref:", sum(res$group == "ref"), "| obs:", sum(res$group == "obs"), "\n")
cat("threshold(95pct of ref):", round(thr95, 4), "\n")
cat("malignant:", sum(res$malignant_95 == "Malignant"),
    "| obs pct:", round(100 * mean(res[group == "obs", malignant_95] == "Malignant"), 1), "%\n")
logmsg("P4g-score2 complete.")
