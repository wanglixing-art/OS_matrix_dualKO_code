# ============================================================
# P4d: CellChat 细胞通讯分析（模板: Molecules 2026 CKAP2 方法学）
# 输入: P4_seurat_annotated.rds (P4a-c 产物)
# 输出: 通讯网络图 + 配体-受体对表 + hub 相关区室信号
# ============================================================
suppressMessages({
  library(CellChat); library(Seurat); library(ggplot2)
  library(patchwork); library(data.table); library(NMF)
})
set.seed(20260928)
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

logmsg("Loading annotated Seurat object ...")
merged <- readRDS(file.path(PROJ, "data/processed/P4_seurat_annotated.rds"))
logmsg("cells:", ncol(merged), "| celltypes:", length(unique(merged$celltype_marker)))

# ---------- 构建 CellChat 对象（用手动矩阵规避 Seurat v5 层兼容问题）----------
logmsg("Extracting data matrix ...")
data.input <- LayerData(merged, layer = "data")   # 单层（P4a 已 JoinLayers）
meta <- data.frame(labels = merged$celltype_marker,
                   malignant = merged$malignant_proxy,
                   row.names = colnames(merged))
logmsg("matrix:", nrow(data.input), "genes x", ncol(data.input), "cells")

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat@DB <- CellChatDB.human   # v2 数据库（Secreted + ECM-Receptor 全量）
logmsg("DB interactions:", nrow(cellchat@DB$interactions))

# ---------- 标准流程 ----------
cellchat <- subsetData(cellchat)
# 关键：subsetData 后 data.signaling 已含 DB 基因子集，
# 全量矩阵 data.input 与 Seurat 对象不再需要——必须释放，
# 否则 42k 细胞下（Seurat ~10GB + 全矩阵 8.4GB + worker 副本）
# 会 OOM 静默崩溃（上一轮在 38% 处无报错死掉）
rm(data.input, merged, meta); invisible(gc())
future::plan("multisession", workers = 2)
# presto 在本机 R 4.6.1 下 lazy-load 静默失败，按官方提示退回标准 Wilcoxon
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
logmsg("OverExpressedInteractions done")
invisible(gc())

cellchat <- computeCommunProb(cellchat)
logmsg("computeCommunProb done")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
logmsg("aggregateNet done")

# 关键：推断完成立即落盘，防止后续绘图环节失败丢失全部计算
saveRDS(cellchat, file.path(PROJ, "data/processed/P4_cellchat.rds"))
logmsg("cellchat RDS saved (post-inference checkpoint)")

groupSize <- as.numeric(table(cellchat@idents))
logmsg("group sizes:"); print(groupSize)

# ---------- 网络级输出 ----------
mat_count <- cellchat@net$count
mat_weight <- cellchat@net$weight
fwrite(as.data.frame(mat_count) %>% tibble::rownames_to_column("source"),
       file.path(TM, "P4_cellchat_net_count.csv"), row.names = FALSE)
fwrite(as.data.frame(mat_weight) %>% tibble::rownames_to_column("source"),
       file.path(TM, "P4_cellchat_net_weight.csv"), row.names = FALSE)

try({
  png(file.path(FG, "P4_cellchat_circle_count.png"), 1500, 1400, res = 150)
  par(mfrow = c(1,1), xpd = TRUE)
  netVisual_circle(mat_count, weight.scale = TRUE, label.edge = FALSE,
                   edge.weight.max = max(mat_count), title.name = "Interaction counts")
  dev.off()
}, silent = TRUE)
try({
  png(file.path(FG, "P4_cellchat_circle_weight.png"), 1500, 1400, res = 150)
  par(mfrow = c(1,1), xpd = TRUE)
  netVisual_circle(mat_weight, weight.scale = TRUE, label.edge = FALSE,
                   edge.weight.max = max(mat_weight), title.name = "Interaction weights (strength)")
  dev.off()
}, silent = TRUE)

try({
  png(file.path(FG, "P4_cellchat_heatmap_count.png"), 1500, 1300, res = 150)
  print(netVisual_heatmap(cellchat, measure = "count",
        color.heatmap = "Blues") + ggtitle("Number of interactions"))
  dev.off()
}, silent = TRUE)
try({
  png(file.path(FG, "P4_cellchat_heatmap_weight.png"), 1500, 1300, res = 150)
  print(netVisual_heatmap(cellchat, measure = "weight",
        color.heatmap = "Reds") + ggtitle("Interaction strength"))
  dev.off()
}, silent = TRUE)

# ---------- 信号角色（outgoing/incoming，直接从 net 矩阵计算，稳健不依赖内部 slot 结构）----------
out_w <- rowSums(cellchat@net$weight, na.rm = TRUE)
in_w  <- colSums(cellchat@net$weight, na.rm = TRUE)
out_n <- rowSums(cellchat@net$count,  na.rm = TRUE)
in_n  <- colSums(cellchat@net$count,  na.rm = TRUE)
sig_roles <- data.frame(
  celltype   = names(out_w),
  out_weight = as.numeric(out_w),
  in_weight  = as.numeric(in_w[match(names(out_w), names(in_w))]),
  out_count  = as.numeric(out_n),
  in_count   = as.numeric(in_n[match(names(out_n), names(in_n))])
)
fwrite(sig_roles, file.path(TM, "P4_cellchat_signaling_roles.csv"))
print(sig_roles)

try({
  png(file.path(FG, "P4_cellchat_roles_scatter.png"), 1400, 1200, res = 150)
  print(netAnalysis_signalingRole_scatter(cellchat,
        label = c("CAF/MSC","Osteoblastic","Proliferating")))
  dev.off()
}, silent = TRUE)

# ---------- 显著通路 ----------
pathways <- cellchat@netP$pathways
logmsg("major signaling pathways (", length(pathways), "):")
print(pathways)
fw <- data.frame(pathway = pathways)
fwrite(fw, file.path(TM, "P4_cellchat_pathways.csv"))

# ---------- hub 区室相关 L-R 对 ----------
# RUNX2-high 区室 = CAF/MSC + Osteoblastic; BUB1-high 区室 = Proliferating
allLR <- subsetCommunication(cellchat)   # 全部 L-R 对
# 兼容不同版本列名（pathway_name vs pathway）
pw_col <- intersect(c("pathway_name","pathway"), colnames(allLR))[1]
logmsg("LR table columns:", paste(colnames(allLR), collapse = ", "))
fwrite(allLR, file.path(TM, "P4_cellchat_all_LR_pairs.csv"))
logmsg("total LR pairs:", nrow(allLR))

hub_sources <- allLR[allLR$source %in% c("CAF/MSC","Osteoblastic","Proliferating"), ]
hub_sources <- hub_sources[order(-hub_sources$prob), ]
fwrite(hub_sources, file.path(TM, "P4_cellchat_LR_from_hub_compartments.csv"))
logmsg("LR pairs from hub compartments:", nrow(hub_sources))
logmsg("Top20 strongest signals from hub compartments:")
print(head(hub_sources[, c("source","target", pw_col, "ligand","receptor","prob")], 20))

# ---------- 关键通路 bubble 图（hub 区室为 source）----------
key_pw <- unique(hub_sources[[pw_col]][hub_sources$prob > quantile(hub_sources$prob, 0.9)])
key_pw <- head(setdiff(key_pw, c(NA, "")), 6)
logmsg("key pathways for bubble:", paste(key_pw, collapse = ", "))
try({
  if (length(key_pw) >= 1) {
    png(file.path(FG, "P4_cellchat_bubble_hub_sources.png"), 2400, 1400, res = 140)
    p <- netVisual_bubble(cellchat, sources.use = c("CAF/MSC","Osteoblastic","Proliferating"),
                          signaling = key_pw, remove.isolate = FALSE)
    print(p)
    dev.off()
  }
}, silent = TRUE)

# ---------- CAF/MSC + Osteoblastic 为 receiver 的信号 ----------
hub_recv <- allLR[allLR$target %in% c("CAF/MSC","Osteoblastic","Proliferating"), ]
hub_recv <- hub_recv[order(-hub_recv$prob), ]
fwrite(hub_recv, file.path(TM, "P4_cellchat_LR_to_hub_compartments.csv"))
logmsg("LR pairs to hub compartments:", nrow(hub_recv))

logmsg("=== P4d SUMMARY ===")
cat("cells: 42784 | pathways:", length(pathways), "\n")
cat("LR pairs total:", nrow(allLR),
    "| from hub comps:", nrow(hub_sources),
    "| to hub comps:", nrow(hub_recv), "\n")
logmsg("P4d complete.")
