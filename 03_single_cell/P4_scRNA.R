# ============================================================
# P4a-c: 单细胞 QC + 整合 + 注释 + 恶性判定 + 双 hub 细胞定位
# 数据: GSE162454 (6 样本) 主数据
# 模板: Molecules 2026;31:1901 方法学（Seurat + Harmony + SingleR + inferCNV 思路）
# ============================================================
suppressMessages({
  library(Seurat); library(harmony); library(SingleR); library(ggplot2)
  library(patchwork); library(data.table)
})
set.seed(20260928)
PROJ <- "D:/projects/OS_matrix_dualKO"
RAW <- file.path(PROJ, "data/raw"); TM <- file.path(PROJ, "results/tables")
FG <- file.path(PROJ, "results/figures")
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
HUBS <- c("BUB1", "RUNX2")

# ---------- 读入 6 个样本 ----------
logmsg("Reading GSE162454 (6 samples) ...")
edir <- file.path(RAW, "GSE162454/extracted")
samples <- list.files(edir, pattern = "_matrix.mtx.gz$")
samples <- sub("_matrix.mtx.gz$", "", samples)
# 组织样本命名: OS_1..OS_6，对应的临床归属（GEO 元数据）
meta_map <- data.frame(
  sample = c("GSM4952363_OS_1","GSM4952364_OS_2","GSM4952365_OS_3",
             "GSM5155198_OS_4","GSM5155199_OS_5","GSM5155200_OS_6"),
  patient = c("OS1","OS2","OS3","OS4","OS5","OS6"))
logmsg("samples found:", paste(samples, collapse = ", "))

objs <- list()
for (s in samples) {
  # 文件带 GSM 前缀，Read10X 只认标准名，故用 ReadMtx 显式指定
  mtx <- ReadMtx(mtx = file.path(edir, paste0(s, "_matrix.mtx.gz")),
                 cells = file.path(edir, paste0(s, "_barcodes.tsv.gz")),
                 features = file.path(edir, paste0(s, "_features.tsv.gz")),
                 feature.column = 2)
  o <- CreateSeuratObject(counts = mtx, project = s, min.cells = 3, min.features = 200)
  o$orig_sample <- s
  o$sample_short <- sub("^GSM[0-9]+_", "", s)
  objs[[s]] <- o
  logmsg("  ", s, ":", ncol(o), "cells")
}

# 合并
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = samples, project = "OS_scRNA")
# Seurat v5: merge 后 counts 分层保存，必须 JoinLayers 合并为单层，
# 否则后续 LayerData(x, layer="data") 会命中多个 data.XXX 层而报错
merged <- JoinLayers(merged)
merged$sample <- merged$sample_short
logmsg("merged:", ncol(merged), "cells x", nrow(merged), "genes")
logmsg("layers after JoinLayers:", paste(Layers(merged[["RNA"]]), collapse = ", "))

# ---------- QC ----------
logmsg("QC ...")
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")
qc_before <- ncol(merged)
p1 <- VlnPlot(merged, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              group.by = "sample", ncol = 3, pt.size = 0)
png(file.path(FG, "P4_QC_violin_before.png"), 2200, 900, res = 140); print(p1); dev.off()
# 模板口径: nFeature 200-6000, mt<20%
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
logmsg("QC:", qc_before, "->", ncol(merged), "cells")
p2 <- VlnPlot(merged, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              group.by = "sample", ncol = 3, pt.size = 0)
png(file.path(FG, "P4_QC_violin_after.png"), 2200, 900, res = 140); print(p2); dev.off()

# ---------- 归一化 + 整合 ----------
logmsg("Normalize + Harmony ...")
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 3000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 40, verbose = FALSE)
merged <- RunHarmony(merged, group.by.vars = "sample", verbose = FALSE, plot_convergence = FALSE)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 1.0, verbose = FALSE)
logmsg("clusters:", length(levels(merged$seurat_clusters)))

png(file.path(FG, "P4_UMAP_clusters.png"), 1400, 1200, res = 150)
print(DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("OS scRNA clusters (Harmony)"))
dev.off()
png(file.path(FG, "P4_UMAP_samples.png"), 1400, 1200, res = 150)
print(DimPlot(merged, reduction = "umap", group.by = "sample") + ggtitle("By sample"))
dev.off()

# ---------- 注释: marker-based + SingleR 双轨（模板用 SingleR）----------
logmsg("Annotation (marker + SingleR) ...")
mk <- list(
  "Osteoblastic"    = c("ALPL","RUNX2","SP7","IBSP","BGLAP","SPP1","COL1A1"),
  "Proliferating"   = c("MKI67","TOP2A","BUB1","CCNB1","UBE2C","CDC20"),
  "CAF/MSC"         = c("PDGFRA","LUM","DCN","COL3A1","COL5A1","ASPN","POSTN"),
  "Myeloid"         = c("LYZ","CD68","CD14","C1QA","C1QB","AIF1","ITGAM","TYROBP"),
  "Tcell"           = c("CD3D","CD3E","CD2","CD8A","IL7R","TRAC"),
  "Bcell"           = c("CD79A","MS4A1","MZB1","IGHG1"),
  "Endothelial"     = c("PECAM1","VWF","CLDN5","RAMP2"),
  "Osteoclast"      = c("CTSK","ACP5","MMP9","TRAP"),
  "NK"              = c("NKG7","GNLY","KLRD1")
)
# 各 cluster 平均表达 marker -> 每列 z-score -> 取最高分类型
avg <- AverageExpression(merged, features = unique(unlist(mk)), group.by = "seurat_clusters",
                         assays = "RNA", slot = "data")$RNA
avg <- as.matrix(avg)
avg_z <- t(scale(t(avg)))            # 基因方向 z-score
avg_z[is.na(avg_z)] <- 0
score <- sapply(mk, function(g) {
  gg <- intersect(g, rownames(avg_z))
  if (length(gg) == 0) return(rep(NA_real_, ncol(avg_z)))
  colMeans(avg_z[gg, , drop = FALSE])
})
score <- as.matrix(score)
ann <- colnames(score)[apply(score, 1, which.max)]
names(ann) <- colnames(avg)
# 关键修正：AverageExpression 的 group.by 列名会被 Seurat 加前缀 "g"
# （因 seurat_clusters 以数字开头，日志: "appending g"），
# 而 merged$seurat_clusters 的值是纯数字 "0","1",...。
# 必须去掉前缀，否则 ann[cl] 全部命中 NA -> 所有细胞变 Unknown。
names(ann) <- sub("^g(?=[0-9])", "", names(ann), perl = TRUE)
logmsg("marker annotation:"); print(table(ann))
fwrite(data.frame(cluster = names(ann), celltype_marker = ann),
       file.path(TM, "P4_cluster_annotation_marker.csv"))
# 显式按 cluster 名映射，避免 NA 索引导致赋值失败
cl <- as.character(merged$seurat_clusters)
merged$celltype_marker <- unname(ann[cl])
n_unknown <- sum(is.na(merged$celltype_marker))
if (n_unknown > 0) logmsg("WARNING: unmapped clusters -> Unknown:", n_unknown, "cells")
merged$celltype_marker[is.na(merged$celltype_marker)] <- "Unknown"
logmsg("celltype_marker distribution:")
print(table(merged$celltype_marker))

png(file.path(FG, "P4_marker_dotplot.png"), 2200, 1200, res = 140)
print(DotPlot(merged, features = unique(unlist(mk)), group.by = "seurat_clusters") +
        RotatedAxis() + ggtitle("Marker expression across clusters"))
dev.off()

# SingleR 参考集（celldex 若不可用则跳过，marker 注释为准）
singler_ok <- requireNamespace("celldex", quietly = TRUE)
if (singler_ok) {
  logmsg("SingleR with celldex::HumanPrimaryCellAtlasData ...")
  ref <- try(celldex::HumanPrimaryCellAtlasData(), silent = TRUE)
  if (!inherits(ref, "try-error")) {
    pred <- SingleR(test = GetAssayData(merged, layer = "data"), ref = ref,
                    labels = ref$label.main, clusters = merged$seurat_clusters)
    sr <- data.frame(cluster = rownames(pred), celltype_SingleR = pred$labels,
                     pruned = pred$pruned.labels)
    fwrite(sr, file.path(TM, "P4_cluster_annotation_SingleR.csv"))
    logmsg("SingleR labels:"); print(table(pred$labels))
    merged$celltype_SingleR <- pred$labels[as.character(merged$seurat_clusters)]
  }
} else logmsg("celldex unavailable -> marker annotation only")

png(file.path(FG, "P4_UMAP_celltypes.png"), 1500, 1200, res = 150)
print(DimPlot(merged, reduction = "umap", group.by = "celltype_marker", label = TRUE, repel = TRUE) +
        ggtitle("Cell types (marker-based)"))
dev.off()

# ---------- 恶性判定（inferCNV 不可用 -> 用成骨谱系 + 增殖双指标近似）----------
logmsg("Malignancy proxy (osteoblastic lineage + proliferation) ...")
merged <- AddModuleScore(merged, features = list(c("ALPL","RUNX2","SP7","IBSP","BGLAP","SPP1","COL1A1")),
                         name = "OsteoScore", seed = 1)
merged <- AddModuleScore(merged, features = list(c("BUB1","TOP2A","MKI67","CCNB1","UBE2C","CDC20","CENPF")),
                         name = "ProlifScore", seed = 1)
# 恶性候选 = 成骨谱系细胞且 OsteoScore/ProlifScore 均高于该类内部中位
osteo_idx <- which(merged$celltype_marker %in% c("Osteoblastic","Proliferating"))
osteo_cells <- colnames(merged)[osteo_idx]
logmsg("osteoblastic+Proliferating lineage cells:", length(osteo_cells))
osteo_osteo <- median(merged$OsteoScore1[osteo_idx])
osteo_prol  <- median(merged$ProlifScore1[osteo_idx])
merged$malignant_proxy <- ifelse(seq_len(ncol(merged)) %in% osteo_idx &
                                 merged$OsteoScore1 > osteo_osteo &
                                 merged$ProlifScore1 > osteo_prol, "Malignant", "Non-malignant")
logmsg("malignant proxy cells:", sum(merged$malignant_proxy == "Malignant"),
       "/", ncol(merged))

# ---------- 双 hub 细胞定位 ----------
logmsg("Dual-hub cell localization ...")
# Seurat v5: 用 LayerData 取 data 层，避免多 layer 报错
hub_mat <- as.matrix(LayerData(merged, layer = "data")[HUBS, , drop = FALSE])
hub_expr <- t(hub_mat)
hub_df <- data.frame(hub_expr, celltype = merged$celltype_marker,
                     malignant = merged$malignant_proxy, sample = merged$sample)
agg_ct <- aggregate(cbind(BUB1, RUNX2) ~ celltype, data = hub_df, FUN = mean)
agg_ct$n_cells <- as.vector(table(hub_df$celltype)[agg_ct$celltype])
agg_ct <- agg_ct[order(-agg_ct$RUNX2), ]
print(agg_ct, row.names = FALSE)
fwrite(agg_ct, file.path(TM, "P4_hub_expr_by_celltype.csv"))

agg_m <- aggregate(cbind(BUB1, RUNX2) ~ malignant, data = hub_df, FUN = mean)
print(agg_m, row.names = FALSE)
fwrite(agg_m, file.path(TM, "P4_hub_expr_by_malignancy.csv"))

png(file.path(FG, "P4_hub_featureplot.png"), 1800, 800, res = 140)
print(FeaturePlot(merged, features = HUBS, reduction = "umap", ncol = 2, order = TRUE))
dev.off()
png(file.path(FG, "P4_hub_violin_celltype.png"), 1800, 800, res = 140)
print(VlnPlot(merged, features = HUBS, group.by = "celltype_marker", pt.size = 0, ncol = 2))
dev.off()

# hub-high vs hub-low 分层比例（模板"按 hub 表达分群"）
merged$RUNX2_expr <- as.numeric(LayerData(merged, layer="data")["RUNX2", ])
merged$BUB1_expr <- as.numeric(LayerData(merged, layer="data")["BUB1", ])
merged$hub_group <- ifelse(merged$RUNX2_expr > median(merged$RUNX2_expr), "RUNX2_high", "RUNX2_low")
prop_tab <- as.data.frame.matrix(table(merged$celltype_marker, merged$hub_group))
prop_tab$celltype <- rownames(prop_tab)
fwrite(prop_tab, file.path(TM, "P4_hub_group_celltype_proportions.csv"))
print(prop_tab, row.names = FALSE)

saveRDS(merged, file.path(PROJ, "data/processed/P4_seurat_annotated.rds"))
logmsg("=== P4a-c SUMMARY ===")
cat("cells:", ncol(merged), "| clusters:", length(levels(merged$seurat_clusters)), "\n")
cat("cell types:", paste(names(table(merged$celltype_marker)), table(merged$celltype_marker), sep = ":", collapse = " | "), "\n")
cat("malignant proxy:", sum(merged$malignant_proxy == "Malignant"), "\n")
logmsg("P4a-c complete.")
