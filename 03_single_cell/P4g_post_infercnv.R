# ============================================================
# P4g-post: inferCNV 正式 CNV 判定 —— 可视化与交叉验证
# 依赖 P4g_infercnv.R 产出的：
#   - results/tables/P4g_percell_cnv_score.csv   (每细胞 CNV 分数 + 恶性判定)
#   - results/infercnv/infercnv.observations.txt (HMM i6 残差矩阵)
#   - results/infercnv/infercnv.references.txt
# 产出：
#   - 主图: CNV 热图（按细胞类型分组）
#   - 恶性比例 by 细胞类型 / by 样本
#   - 与代理法交叉表（验证一致性）
#   - 双 hub 表达 vs 恶性状态
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

ROOT   <- "D:/projects/OS_matrix_dualKO"
TM     <- file.path(ROOT, "results/tables")
FG     <- file.path(ROOT, "results/figures")
INFCNV <- file.path(ROOT, "results/infercnv")
dir.create(FG, showWarnings = FALSE, recursive = TRUE)

LOG <- file.path(ROOT, "logs/P4g_post.log")
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)

HUBS <- c("BUB1", "RUNX2")
OBS_CTS <- c("Osteoblastic", "Proliferating")

# ---------- 1) 读 CNV 分数表 ----------
score_f <- file.path(TM, "P4g_percell_cnv_score.csv")
if (!file.exists(score_f)) stop("缺少 P4g_percell_cnv_score.csv —— P4g 尚未跑完")
res <- fread(score_f)
logmsg("loaded per-cell CNV scores:", nrow(res), "cells")
logmsg("malignant:", sum(res$malignant == "Malignant"))

# ---------- 2) 与代理法交叉表 ----------
proxy_f <- file.path(TM, "P4_hub_expr_by_malignancy.csv")
ann_f   <- file.path(TM, "P4_cluster_annotation_SingleR.csv")

# 从 Seurat 对象取代理法恶性标签（与 P4a-c 口径一致）
seu_f <- file.path(ROOT, "data/processed/P4_seurat_annotated.rds")
if (file.exists(seu_f)) {
  suppressPackageStartupMessages(library(Seurat))
  seu <- readRDS(seu_f)
  md <- seu@meta.data
  md$cell <- rownames(md)
  keep <- c("cell", "malignant_proxy", "celltype", "orig.ident")
  keep <- intersect(keep, colnames(md))
  mdt <- as.data.table(md[, keep, drop = FALSE])
  if ("malignant_proxy" %in% colnames(mdt)) {
    mg <- merge(res, mdt[, .(cell, malignant_proxy)], by = "cell", all.x = TRUE)
    ct <- mg[!is.na(malignant_proxy),
             .N, by = .(inferCNV = malignant, proxy = malignant_proxy)]
    fwrite(ct, file.path(TM, "P4g_infercnv_vs_proxy_crosstab.csv"))
    logmsg("crosstab (inferCNV vs proxy):")
    print(ct)
    # 一致性指标
    tab <- dcast(ct, inferCNV ~ proxy, value.var = "N", fill = 0)
    logmsg("crosstab wide:")
    print(tab)
  } else {
    logmsg("WARN: Seurat meta 无 malignant_proxy 列，跳过交叉表")
  }
  rm(seu, md); gc()
} else {
  logmsg("WARN: 未找到 P4_seurat_annotated.rds，跳过交叉表")
}

# ---------- 3) 恶性比例 by 细胞类型 / by 样本 ----------
ct_sum <- res[, .(n = .N,
                  malignant = sum(malignant == "Malignant"),
                  pct = round(100 * mean(malignant == "Malignant"), 1)),
              by = celltype][order(-pct)]
fwrite(ct_sum, file.path(TM, "P4g_malignant_pct_by_celltype.csv"))
logmsg("malignant pct by celltype:")
print(ct_sum)

sp_sum <- res[, .(n = .N,
                  malignant = sum(malignant == "Malignant"),
                  pct = round(100 * mean(malignant == "Malignant"), 1)),
              by = sample][order(sample)]
fwrite(sp_sum, file.path(TM, "P4g_malignant_pct_by_sample.csv"))
logmsg("malignant pct by sample:")
print(sp_sum)

# 恶性细胞在 obs/ref 分组的分布（应几乎全部落在 obs）
if ("group" %in% colnames(res)) {
  grp <- res[, .(n = .N, malignant = sum(malignant == "Malignant"),
                 pct = round(100 * mean(malignant == "Malignant"), 1)), by = group]
  logmsg("by group (ref vs obs):")
  print(grp)
}

# ---------- 4) 图：恶性比例条形图 ----------
pal_ct <- c("Osteoblastic" = "#D7263D", "Proliferating" = "#F46036",
            "CAF/MSC" = "#2E86AB", "Myeloid" = "#6A994E", "Tcell" = "#BC4B51",
            "Endothelial" = "#8D99AE", "Bcell" = "#7209B7", "NK" = "#F2A900",
            "Osteoclast" = "#4A4E69", "ref" = "#B0B0B0")

p1 <- ggplot(ct_sum, aes(x = reorder(celltype, pct), y = pct, fill = celltype)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = paste0(pct, "%\n(n=", n, ")")), hjust = -0.08, size = 3) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = pal_ct, guide = "none") +
  scale_y_continuous(limits = c(0, 118), expand = c(0, 0)) +
  labs(title = "inferCNV malignant fraction by cell type",
       subtitle = paste0("threshold = 95th pct of reference CNV score"),
       x = NULL, y = "% Malignant") +
  theme_bw(base_size = 12) +
  theme(plot.margin = margin(6, 48, 6, 6),
        panel.grid.minor = element_blank())

p2 <- ggplot(sp_sum, aes(x = sample, y = pct, fill = sample)) +
  geom_col(width = 0.7, color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(limits = c(0, max(sp_sum$pct) * 1.25), expand = c(0, 0)) +
  labs(title = "Malignant fraction per sample", x = NULL, y = "% Malignant") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

p3 <- ggplot(res, aes(x = reorder(celltype, cnv_score, median), y = cnv_score, fill = malignant)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.3, alpha = 0.9) +
  scale_fill_manual(values = c("Malignant" = "#D7263D", "Non-malignant" = "#2E86AB")) +
  labs(title = "Per-cell CNV score distribution", x = NULL, y = "CNV score") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())

gg <- (p1 | p2) / (p3) + plot_layout(heights = c(1, 1))
ggsave(file.path(FG, "P4g_infercnv_malignancy_summary.png"), gg,
       width = 13, height = 10, dpi = 300)
logmsg("saved P4g_infercnv_malignancy_summary.png")

# ---------- 5) CNV 热图（obs 细胞，按细胞类型排序） ----------
hm_cells <- res[malignant == "Malignant" | celltype %in% OBS_CTS]
# 每个细胞类型最多取 400 细胞，控制位图尺寸
set.seed(42)
hm_sel <- hm_cells[, .SD[sample(.N, min(.N, 400))], by = celltype]
logmsg("heatmap cells:", nrow(hm_sel))

obs_f <- file.path(INFCNV, "infercnv.observations.txt")
if (file.exists(obs_f)) {
  logmsg("reading observation matrix (may be large) ...")
  m <- as.matrix(fread(obs_f, sep = "\t", header = TRUE), rownames = 1)
  cols <- intersect(hm_sel$cell, colnames(m))
  logmsg("heatmap matrix:", length(cols), "cells x", nrow(m), "genes")
  sub <- m[, cols, drop = FALSE]
  # 按细胞类型的（良性 -> 恶性）顺序排列
  ord_ct <- c("CAF/MSC", "Myeloid", "Tcell", "Endothelial", "Bcell", "NK",
              "Osteoclast", "Osteoblastic", "Proliferating")
  cell_ct <- hm_sel$celltype[match(cols, hm_sel$cell)]
  o <- order(match(cell_ct, ord_ct), cell_ct)
  sub <- sub[, o, drop = FALSE]
  cell_ct <- cell_ct[o]

  # 行=基因 需按染色体物理顺序（gene_order 已排序，直接保留）
  df <- as.data.table(as.table(sub))
  setnames(df, c("gene", "cell", "val"))
  df[, cell := factor(cell, levels = colnames(sub))]
  df[, gene_idx := match(gene, rownames(sub))]

  # 用 tile 绘制（比 pheatmap 更快且可控）
  ncell <- length(cols); ngene <- nrow(sub)
  # 降采样基因轴以减少内存（每 1 基因取 1，最多 20000）
  if (ngene > 20000) {
    keep_idx <- round(seq(1, ngene, length.out = 20000))
    df <- df[gene_idx %in% keep_idx]
    ngene <- 20000
  }
  vmax <- quantile(abs(df$val - 1), 0.995, na.rm = TRUE)
  p_hm <- ggplot(df, aes(x = cell, y = gene_idx, fill = val)) +
    geom_raster() +
    scale_fill_gradient2(low = "#00008B", mid = "white", high = "#8B0000",
                         midpoint = 1, limits = c(1 - vmax, 1 + vmax),
                         oob = scales::squish, name = "CNV") +
    labs(title = "inferCNV heatmap (HMM i6) — malignant & osteoblastic/proliferating cells",
         subtitle = "blue = copy loss, red = copy gain",
         x = "Cells (grouped by cell type)", y = "Genes (chromosome order)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank())

  ggsave(file.path(FG, "P4g_infercnv_heatmap_obs.png"), p_hm,
         width = 14, height = 6, dpi = 300)
  logmsg("saved P4g_infercnv_heatmap_obs.png")

  # 细胞类型分隔线注释表
  blk <- data.table(celltype = cell_ct)[, .(n = .N), by = celltype]
  blk[, cum := cumsum(n)]
  fwrite(blk, file.path(TM, "P4g_heatmap_celltype_blocks.csv"))

  rm(m, sub, df); gc()
} else {
  logmsg("WARN: infercnv.observations.txt 不存在，跳过热图")
}

# ---------- 6) 双 hub 表达 vs inferCNV 恶性状态 ----------
hub_f <- file.path(TM, "P4_hub_expr_by_malignancy.csv")
if (file.exists(seu_f) && file.exists(hub_f)) {
  suppressPackageStartupMessages(library(Seurat))
  seu <- readRDS(seu_f)
  # 取双 hub 表达
  hb <- FetchData(seu, vars = HUBS, layer = "data")
  hb$cell <- rownames(hb)
  hbt <- as.data.table(hb)
  mg <- merge(hbt, res[, .(cell, malignant, celltype)], by = "cell", all.x = TRUE)
  mg <- mg[!is.na(malignant)]
  long <- melt(mg, id.vars = c("cell", "malignant", "celltype"),
               measure.vars = HUBS, variable.name = "hub", value.name = "expr")
  fwrite(mg, file.path(TM, "P4g_hub_expr_vs_infercnv_malignant.csv"))
  ph <- ggplot(long, aes(x = malignant, y = expr, fill = malignant)) +
    geom_violin(scale = "width", alpha = 0.9, linewidth = 0.3, color = "grey25") +
    geom_boxplot(width = 0.16, outlier.size = 0.2, fill = "white", linewidth = 0.3) +
    facet_wrap(~hub, scales = "free_y") +
    scale_fill_manual(values = c("Malignant" = "#D7263D", "Non-malignant" = "#2E86AB")) +
    labs(title = "Hub gene expression vs inferCNV malignancy",
         x = NULL, y = "log-normalized expression") +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  ggsave(file.path(FG, "P4g_hub_vs_infercnv_malignant.png"), ph,
         width = 8, height = 4.5, dpi = 300)
  logmsg("saved P4g_hub_vs_infercnv_malignant.png")
  rm(seu, hb, hbt, mg, long); gc()
}

logmsg("=== P4g-post complete ===")
