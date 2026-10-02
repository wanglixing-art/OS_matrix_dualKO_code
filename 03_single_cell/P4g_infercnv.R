# ============================================================
# P4g: inferCNV 正式 CNV 恶性判定（替代成骨双评分代理法）
# 数据: GSE162454 (6 样本 42,784 细胞)
# 参照: 非成骨谱系（Tcell/Bcell/NK/Myeloid/Endothelial/Osteoclast）
# 目标: Osteoblastic + Proliferating（成骨谱系，疑似恶性）
# 方法: inferCNV 1.28.0, HMM i6, 降采样控制规模
# ============================================================
suppressMessages({
  library(infercnv); library(Seurat); library(data.table)
})
set.seed(20260928)
PROJ <- "D:/projects/OS_matrix_dualKO"
RAW <- file.path(PROJ, "data/raw"); TM <- file.path(PROJ, "results/tables")
FG <- file.path(PROJ, "results/figures")
OUTDIR <- file.path(PROJ, "results/infercnv")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

GENE_ORDER <- file.path(RAW, "infercnv_ref/gene_order_hg38.txt")
MAX_REF_PER_CT <- 800      # 每类参照细胞上限（降采样，控内存）
MAX_OBS_PER_CT <- 3000     # 每个目标细胞类型上限

# ---------- 1) 载入并抽样 ----------
logmsg("Loading annotated Seurat object ...")
merged <- readRDS(file.path(PROJ, "data/processed/P4_seurat_annotated.rds"))
logmsg("total cells:", ncol(merged))

ct <- as.character(merged$celltype_marker)
OBS_CTS <- c("Osteoblastic","Proliferating")
REF_CTS <- c("Tcell","Bcell","NK","Myeloid","Endothelial","Osteoclast")

ref_idx <- unlist(lapply(REF_CTS, function(x) {
  i <- which(ct == x)
  if (length(i) > MAX_REF_PER_CT) sample(i, MAX_REF_PER_CT) else i
}))
obs_idx <- unlist(lapply(OBS_CTS, function(x) {
  i <- which(ct == x)
  if (length(i) > MAX_OBS_PER_CT) sample(i, MAX_OBS_PER_CT) else i
}))
keep_idx <- c(ref_idx, obs_idx)
logmsg("sampled cells:", length(keep_idx),
       "| ref:", length(ref_idx), "| obs:", length(obs_idx))

sub <- subset(merged, cells = colnames(merged)[keep_idx])
sub_ct <- as.character(sub$celltype_marker)
names(sub_ct) <- colnames(sub)
sub_sample <- as.character(sub$sample)
names(sub_sample) <- colnames(sub)
rm(merged); invisible(gc())

# ---------- 2) 构建计数矩阵 + 注释 ----------
logmsg("Building counts matrix ...")
counts <- as.matrix(LayerData(sub, layer = "counts"))
logmsg("matrix:", nrow(counts), "genes x", ncol(counts), "cells")

# 注释：观测组（成骨谱系）vs 参照组
ann <- data.frame(
  cell = colnames(sub),
  group = ifelse(as.character(sub$celltype_marker) %in% OBS_CTS, "obs", "ref"),
  stringsAsFactors = FALSE
)
logmsg("annotation:"); print(table(ann$group))
# 必须制表符分隔且无表头（inferCNV 要求）；此前用 fwrite 默认逗号分隔导致解析失败
ann_file <- file.path(TM, "P4g_infercnv_annotation.txt")
fwrite(ann, ann_file, sep = "\t", col.names = FALSE, quote = FALSE)

rm(sub); invisible(gc())

# ---------- 3) 运行 inferCNV ----------
logmsg("Creating infercnv object ...")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file   = ann_file,
  delim              = "\t",
  gene_order_file    = GENE_ORDER,
  ref_group_names    = "ref"
)
rm(counts); invisible(gc())

logmsg("Running inferCNV (HMM i6, threads=4) ...")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff        = 0.1,          # 单细胞数据常规阈值
  out_dir       = OUTDIR,
  cluster_by_groups = FALSE,
  denoise       = TRUE,
  HMM           = TRUE,
  HMM_type      = "i6",
  analysis_mode = "subclusters",
  num_threads   = 4,
  plot_steps    = FALSE,
  no_plot       = FALSE,        # 保留 CNV 热图（论文主图）
  sd_amplifier  = 1.5,
  noise_logistic= TRUE,
  write_expr_matrix = TRUE,     # 输出残差矩阵供打分
  save_rds      = TRUE
)
logmsg("inferCNV run completed")

# ---------- 4) 每细胞 CNV 分数 -> 恶性判定 ----------
logmsg("Extracting per-cell CNV scores ...")
# 用最细的 HMM 残差/表达矩阵：从 out_dir 读取
expr_file <- file.path(OUTDIR, "infercnv.references.txt")
obs_file  <- file.path(OUTDIR, "infercnv.observations.txt")
if (!file.exists(obs_file)) {
  obs_file <- file.path(OUTDIR, "expr.infercnv.dat")
}
logmsg("obs file:", obs_file, "exists:", file.exists(obs_file))

mat <- as.matrix(fread(obs_file, sep = "\t", header = TRUE), rownames = 1)
ref_mat <- as.matrix(fread(expr_file, sep = "\t", header = TRUE), rownames = 1)

# CNV 分数定义：每细胞 |修饰表达 - 1| 的均值（偏离中性的程度）
cnv_score <- colMeans(abs(mat - 1), na.rm = TRUE)
ref_score <- colMeans(abs(ref_mat - 1), na.rm = TRUE)

# 阈值：参照细胞分数的 95 分位作为恶性判定阈值（保守）
thr <- as.numeric(quantile(ref_score, 0.95, na.rm = TRUE))
logmsg("ref CNV score median:", round(median(ref_score),4),
       "| 95th pct threshold:", round(thr,4))

res <- data.table(cell = names(cnv_score), cnv_score = as.numeric(cnv_score))
res[, malignant := ifelse(cnv_score > thr, "Malignant", "Non-malignant")]
res[, celltype := sub_ct[cell]]
res[, sample := sub_sample[cell]]
fwrite(res, file.path(TM, "P4g_percell_cnv_score.csv"))

logmsg("malignant cells:", sum(res$malignant == "Malignant"), "/", nrow(res))
logmsg("by celltype:")
print(res[, .(n = .N, malignant = sum(malignant == "Malignant"),
              pct = round(100*mean(malignant == "Malignant"),1)), by = celltype])
logmsg("by sample (obs only):")
print(res[celltype %in% OBS_CTS,
          .(n = .N, malignant = sum(malignant == "Malignant"),
            pct = round(100*mean(malignant == "Malignant"),1)), by = sample][order(sample)])

logmsg("=== P4g SUMMARY ===")
cat("ref cells:", length(ref_score), "| obs cells:", length(cnv_score), "\n")
cat("threshold:", round(thr,4), "\n")
cat("malignant:", sum(res$malignant == "Malignant"), "\n")
logmsg("P4g complete.")
