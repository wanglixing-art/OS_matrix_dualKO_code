# ============================================================
# P4e: CellChat 6 样本 per-sample 对比分析（CellChat comparison 模式）
# 输入: P4_seurat_annotated.rds（含 6 样本）
# 输出: 各样本独立 CellChat 对象 + merged 对比 + 差异通路排名
# ============================================================
suppressMessages({
  library(CellChat); library(Seurat); library(ggplot2)
  library(patchwork); library(data.table)
})
set.seed(20260928)
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
HUBS <- c("BUB1", "RUNX2")

logmsg("Loading annotated Seurat object ...")
merged <- readRDS(file.path(PROJ, "data/processed/P4_seurat_annotated.rds"))
samples <- sort(unique(as.character(merged$sample)))
logmsg("samples:", paste(samples, collapse = ", "))

data.all <- LayerData(merged, layer = "data")
meta.all <- data.frame(labels   = merged$celltype_marker,
                       samples  = merged$sample,
                       row.names = colnames(merged))
rm(merged); invisible(gc())

# ---------- 逐样本推断 ----------
object.list <- list()
for (s in samples) {
  logmsg("===== sample:", s, "=====")
  idx <- which(meta.all$samples == s)
  ct_tab <- table(meta.all$labels[idx])
  # 过滤掉细胞数 <10 的类型（CellChat 要求每组 >=10 细胞）
  keep_ct <- names(ct_tab)[ct_tab >= 10]
  keep_idx <- idx[meta.all$labels[idx] %in% keep_ct]
  logmsg("  cells:", length(keep_idx), "| celltypes kept:", length(keep_ct),
         "|", paste(paste0(keep_ct, "(", ct_tab[keep_ct], ")"), collapse = " "))

  m <- data.all[, keep_idx]
  mt <- meta.all[keep_idx, , drop = FALSE]
  cc <- createCellChat(object = m, meta = mt, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  object.list[[s]] <- cc
  logmsg("  ", s, " inference done | pathways:", length(cc@netP$pathways))
  saveRDS(cc, file.path(PROJ, paste0("data/processed/P4e_cellchat_", s, ".rds")))
  invisible(gc())
}
rm(data.all, meta.all); invisible(gc())

# ---------- 合并对比 ----------
logmsg("Merging per-sample CellChat objects ...")
# lift 到共同基因空间（不同样本识别的信号基因取交集，CellChat 官方推荐）
object.list <- liftCellChat(object.list, group = "0")
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))
logmsg("merged group size:", length(cellchat.merged@net), "samples")

# ---------- 各样本总体交互强度 ----------
gg1 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = c(1,2,3,4,5,6))
gg2 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = c(1,2,3,4,5,6),
                           measure = "weight")
png(file.path(FG, "P4e_compareInteractions.png"), 1600, 800, res = 140)
print(gg1 + gg2)
dev.off()

# ---------- 差异通路（各样本 sender 强度矩阵）----------
# 用 subsetCommunication 的 pathway_name 逐样本汇总通路权重
all_pw <- unique(unlist(lapply(object.list, function(x) x@netP$pathways)))
pw_mat <- matrix(NA_real_, nrow = length(all_pw), ncol = length(object.list),
                 dimnames = list(all_pw, names(object.list)))
for (s in names(object.list)) {
  cc <- object.list[[s]]
  lr_map <- subsetCommunication(cc)          # 含 pathway_name
  agg <- as.data.table(lr_map)[, .(weight = sum(prob, na.rm = TRUE)), by = pathway_name]
  pw_mat[agg$pathway_name, s] <- agg$weight
}
pw_dt <- as.data.frame(pw_mat)
pw_dt$pathway <- rownames(pw_dt)
pw_dt$mean_weight <- rowMeans(pw_mat, na.rm = TRUE)
pw_dt$sd_weight <- apply(pw_mat, 1, sd, na.rm = TRUE)
pw_dt <- pw_dt[order(-pw_dt$mean_weight), ]
fwrite(pw_dt, file.path(TM, "P4e_pathway_weight_by_sample.csv"))
logmsg("Top20 pathways by mean weight across samples:")
print(head(pw_dt[, c("pathway", samples, "mean_weight")], 20))

# ---------- hub 区室 per-sample 信号一致性 ----------
hub_ct <- c("CAF/MSC","Osteoblastic","Proliferating")
consist <- list()
for (s in names(object.list)) {
  cc <- object.list[[s]]
  w <- cc@net$weight
  out_w <- rowSums(w, na.rm = TRUE)
  for (ct in intersect(hub_ct, names(out_w))) {
    consist[[length(consist)+1]] <- data.frame(
      sample = s, celltype = ct, out_weight = as.numeric(out_w[ct]))
  }
}
consist_dt <- rbindlist(consist)
fwrite(consist_dt, file.path(TM, "P4e_hub_outgoing_by_sample.csv"))
logmsg("hub compartment outgoing signal per sample:")
print(dcast(consist_dt, celltype ~ sample, value.var = "out_weight"))

# hub 区室是否在各样本中稳定为最强发送者
rank_dt <- consist_dt[order(sample, -out_weight)]
rank_dt[, rank := seq_len(.N), by = sample]
fwrite(rank_dt, file.path(TM, "P4e_hub_outgoing_rank_by_sample.csv"))

png(file.path(FG, "P4e_hub_outgoing_by_sample.png"), 1700, 1000, res = 150)
p <- ggplot(consist_dt, aes(sample, out_weight, fill = celltype)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c("CAF/MSC" = "#E31A1C",
                               "Osteoblastic" = "#FF7F00",
                               "Proliferating" = "#1F78B4")) +
  labs(x = NULL, y = "Outgoing signaling weight",
       title = "Hub-compartment signaling across the 6 samples", fill = "") +
  theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))
print(p); dev.off()

# ---------- RUNX2-high 基质 vs BUB1-high 增殖：各样本对比 ----------
# 关键叙事：RUNX2 区室（基质）发送强度是否稳定高于 BUB1 区室（增殖）
pivot <- dcast(consist_dt, sample ~ celltype, value.var = "out_weight")
if (all(c("CAF/MSC","Proliferating") %in% colnames(pivot))) {
  pivot[, stroma_vs_prolif := CAF/MSC / Proliferating]
  fwrite(pivot, file.path(TM, "P4e_stroma_vs_proliferation_ratio.csv"))
  logmsg("stroma(CAF/MSC) / proliferation outgoing ratio per sample:")
  print(pivot[, .(sample, `CAF/MSC`, Proliferating, stroma_vs_prolif)])
}

saveRDS(cellchat.merged, file.path(PROJ, "data/processed/P4e_cellchat_merged.rds"))
logmsg("=== P4e SUMMARY ===")
cat("samples:", length(object.list), "\n")
cat("pathways union:", length(all_pw), "\n")
logmsg("P4e complete.")
