# =============================================================================
# P7b_fig1_design.R  -  主图 1：研究设计与数据流总览（流程图）
# -----------------------------------------------------------------------------
# 输出: results/figures/F1_design.pdf / .png （矢量优先，17cm x 22cm 纵向）
# 布局引擎: 坐标系 1 unit = 1 mm（xlim 0-170, ylim 0-220），盒宽/字号经过匹配，
#          文本行宽 <= 盒宽 - 6mm，杜绝溢出。
# 注: 不使用 em dash，统一 ASCII "-"（本机设备缺 Unicode 字形，见 P4e 经验）
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

out_dir <- "D:/projects/OS_matrix_dualKO/results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- 配色
COL_DATA   <- "#4C72B0"   # 数据源（蓝）
COL_POOL   <- "#55A868"   # 候选池（绿）
COL_HUB    <- "#C44E52"   # 枢纽（红）
COL_SC     <- "#8172B2"   # 单细胞（紫）
COL_KO     <- "#E8D28A"   # 虚拟敲除（黄）
COL_SIG    <- "#64B5CD"   # 签名（青）
COL_VAL    <- "#DD8452"   # 验证（橙）
COL_MICRO  <- "#937860"   # 微环境（棕）
COL_ANCHOR <- "#DA8BC3"   # 转化锚点（粉）

FS <- 2.55   # 字号（ggplot size, mm）
CH <- 1.30   # 估平均字符宽 mm（该字号下）

box <- function(cx, cy, w, h, lines, fill, txt = "white", fs = FS) {
  lab <- paste(lines, collapse = "\n")
  maxch <- max(nchar(lines))
  stopifnot(maxch * CH <= w - 4)   # 文本必须放得下
  list(
    annotate("rect", xmin = cx - w/2, xmax = cx + w/2,
             ymin = cy - h/2, ymax = cy + h/2,
             fill = fill, colour = "white", linewidth = 0.6),
    annotate("text", x = cx, y = cy, label = lab, colour = txt,
             size = fs, lineheight = 1.25)
  )
}

arrowv <- function(x, y1, y2, x2 = x) {
  annotate("segment", x = x, xend = x2, y = y1, yend = y2,
           arrow = arrow(length = unit(0.13, "cm"), type = "closed"),
           colour = "grey30", linewidth = 0.5)
}

p <- ggplot() + xlim(0, 170) + ylim(0, 220) + theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

# 列中心: c1=30, c2=85, c3=140（侧盒宽 50）; 宽盒跨 c1-c2: 中心 57.5 宽 105
C1 <- 30; C2 <- 85; C3 <- 140
W  <- 50; WW <- 105; WF <- 160
CXW <- 57.5; CXF <- 85

# ---------------- Row 1 (y 195-215): 数据源
p <- p + box(C1, 205, W, 20, c("Bulk transcriptome",
                               "TARGET-OS n = 88",
                               "(EFS 85 / OS 86)",
                               "Masked MAF n = 135"), COL_DATA)
p <- p + box(C2, 205, W, 20, c("Validation cohorts",
                               "GSE42352 84T+15N",
                               "GSE21257 53 / GSE39055 36",
                               "GSE33382 53 / GSE87624 33"), COL_DATA)
p <- p + box(C3, 205, W, 20, c("Single-cell RNA-seq",
                               "GSE162454 (10x)",
                               "6 samples",
                               "42,784 cells after QC"), COL_SC)

# ---------------- Row 2 (y 168-188)
p <- p + box(CXW, 178, WW, 20, c("limma differential expression + WGCNA",
                                 "556 up / 484 down (|log2FC| > 0.585)",
                                 "21 modules; 7 outcome-associated"), COL_POOL)
p <- p + box(C3, 178, W, 20, c("scRNA-seq processing",
                               "Harmony integration",
                               "SingleR annotation",
                               "9 cell types"), COL_SC)
p <- p + arrowv(C1, 195, 188.5) + arrowv(C2, 195, 188.5) + arrowv(C3, 195, 188.5)

# ---------------- Row 3 (y 141-161)
p <- p + box(CXW, 151, WW, 20, c("Candidate pool",
                                 "299 tumour-upregulated genes",
                                 "intersecting outcome modules"), COL_POOL)
p <- p + box(C3, 151, W, 20, c("inferCNV malignancy",
                               "95th-pct reference threshold",
                               "(thr = 0.050)"), COL_SC)
p <- p + arrowv(CXW, 168, 161.5) + arrowv(C3, 168, 161.5)

# ---------------- Row 4 (y 114-134)
p <- p + box(CXW, 124, WW, 20, c("STRING PPI (166 nodes / 460 edges)",
                                 "+ cytoHubba centrality + Cox",
                                 "two orthogonal hubs nominated"), COL_HUB)
p <- p + box(C3, 124, W, 20, c("Compartment mapping",
                               "RUNX2: CAF/MSC + osteoblast",
                               "BUB1: proliferating cells"),
             COL_SC)
p <- p + arrowv(CXW, 141, 134.5) + arrowv(C3, 141, 134.5)

# ---------------- Row 5 (y 87-107): 双枢纽 + CellChat
p <- p + box(C1, 97, W, 20, c("BUB1 (hub 1)",
                              "spindle checkpoint",
                              "HR 1.39/SD, p = 0.0385"), COL_HUB)
p <- p + box(C2, 97, W, 20, c("RUNX2 (hub 2)",
                              "osteogenic master TF",
                              "HR 1.55/SD, p = 0.0176"), COL_HUB)
p <- p + box(C3, 97, W, 20, c("CellChat",
                              "85 pathways / 2,104 LR pairs",
                              "CAF/MSC top sender 6/6",
                              "COLLAGEN 2.3x runner-up"), COL_SC)
p <- p + arrowv(C1, 114, 107.5) + arrowv(C2, 114, 107.5) + arrowv(C3, 114, 107.5)

# ---------------- Row 6 (y 60-80): 双 KO（核心）
p <- p + box(CXW, 70, WW, 20, c("scTenifoldKnk parallel virtual KO",
                                "RUNX2-KO (stroma): 20 DRGs - TGF-beta / AP-1 / SLRP",
                                "BUB1-KO (prolif.): 17 DRGs - myeloid / phagocytosis / compl.",
                                "shared DRG = 0 (orthogonal)"), COL_KO, txt = "grey15")
p <- p + box(C3, 70, W, 20, c("Microenvironment",
                              "ESTIMATE stromal rho = -0.55",
                              "TMB rho = 0.308 (n = 72)"), COL_MICRO)
p <- p + arrowv(C1, 87, 80.5) + arrowv(C2, 87, 80.5) + arrowv(C3, 87, 80.5)

# ---------------- Row 7 (y 33-53): 签名 + 锚点
p <- p + box(CXW, 43, WW, 20, c("LASSO-Cox signature (14 genes)",
                                "train: C-index 0.796, HR 5.43/SD",
                                "p = 1.02e-11; AUC 0.847-0.860"), COL_SIG, txt = "grey15")
p <- p + box(C3, 43, W, 20, c("Translational anchors",
                              "BUB1: BAY 1816032",
                              "RUNX2: ChIP 2,339 sites"), COL_ANCHOR, txt = "grey15")
p <- p + arrowv(CXW, 60, 53.5) + arrowv(C3, 60, 53.5)

# ---------------- Row 8 (y 6-26): 外部验证
p <- p + box(CXF, 16, WF, 20, c("External validation (frozen coefficients)",
                                "GSE39055: C-index 0.696, HR 2.69, p = 0.034",
                                "metastasis: GSE21257 p = 0.0021 / GSE33382 p = 0.0017",
                                "GSE87624: direction consistent (n = 9 met, exploratory)"),
             COL_VAL)
p <- p + arrowv(CXW, 33, 26.5) + arrowv(C3, 33, 26.5)

ggsave(file.path(out_dir, "F1_design.pdf"), p,
       width = 17, height = 22, units = "cm", device = cairo_pdf)
ggsave(file.path(out_dir, "F1_design.png"), p,
       width = 17, height = 22, units = "cm", dpi = 300)

cat("F1_design.pdf / .png written to", out_dir, "\n")
