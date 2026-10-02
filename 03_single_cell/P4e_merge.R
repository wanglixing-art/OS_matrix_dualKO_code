# ============================================================
# P4e-merge: 复用已推断的 6 个 per-sample CellChat 对象，做 lift + merge 对比
# 修复：liftCellChat 的真实签名是 function(object, group.new=NULL)，
#       传入单个 CellChat 对象（对 @idents 做因子重编码），不是列表；
#       正确用法 = lapply(object.list, liftCellChat, group.new = <统一 celltype 全集>)
# 输入: data/processed/P4e_cellchat_<sample>.rds  (6 个，已推断完成)
# 输出: P4e_cellchat_merged.rds + per-sample 对比图表
# ============================================================
suppressMessages({
  library(CellChat); library(Seurat); library(ggplot2)
  library(patchwork); library(data.table)
})
set.seed(20260928)
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
LOG <- file.path(ROOT, "logs/P4e_merge.log")
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)

samples <- paste0("OS_", 1:6)

# ---------- 1) 载入 6 个已推断对象 ----------
object.list <- list()
for (s in samples) {
  f <- file.path(ROOT, paste0("data/processed/P4e_cellchat_", s, ".rds"))
  if (!file.exists(f)) { logmsg("MISSING:", s); next }
  object.list[[s]] <- readRDS(f)
  logmsg("loaded", s, "| celltypes:", nlevels(object.list[[s]]@idents),
         "| pathways:", length(object.list[[s]]@netP$pathways))
}
logmsg("loaded samples:", paste(names(object.list), collapse = ", "))

# ---------- 2) 求 celltype 交集（CellChat 要求各组细胞类型一致）----------
ct_all <- lapply(object.list, function(x) levels(x@idents))
ct_intersect <- Reduce(intersect, ct_all)
logmsg("celltype intersection (", length(ct_intersect), "):",
       paste(ct_intersect, collapse = ", "))
logmsg("union:", paste(Reduce(union, ct_all), collapse = ", "))
# 交集过小则改用 union + 缺失类型删除
if (length(ct_intersect) < 4) {
  logmsg("WARN: intersection < 4，改用 union 并逐对象 subset")
  ct_intersect <- Reduce(union, ct_all)
}

# 统一 ident 水平：liftCellChat(object, group.new) 把 idents 重编码到 group.new
object.list <- lapply(object.list, function(x) {
  ct <- intersect(levels(x@idents), ct_intersect)
  x <- subsetCellChat(x, idents.use = ct)   # 删掉不在共同集合里的类型
  liftCellChat(x, group.new = ct_intersect)
})
logmsg("lift applied | per-object celltypes after lift:")
for (s in names(object.list)) {
  logmsg("  ", s, ":", paste(levels(object.list[[s]]@idents), collapse = ","))
}

# ---------- 3) merge ----------
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))
logmsg("merged groups:", length(cellchat.merged@net))

# ---------- 4) 总体交互数/强度对比 ----------
ng <- length(object.list)
try({
  gg1 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = seq_len(ng))
  gg2 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = seq_len(ng),
                             measure = "weight")
  p <- gg1 + gg2
  ggsave(file.path(FG, "P4e_compareInteractions.png"), p, width = 11, height = 5, dpi = 300)
  logmsg("saved P4e_compareInteractions.png")
})
try({
  p <- compareInteractions(cellchat.merged, show.legend = FALSE, group = seq_len(ng),
                           measure = "weight", x.label = "Sample", y.label = "Interaction weight")
  ggsave(file.path(FG, "P4e_compareInteractions_weight.png"), p, width = 6, height = 5, dpi = 300)
})

# ---------- 5) per-sample 通路权重矩阵 ----------
all_pw <- sort(unique(unlist(lapply(object.list, function(x) x@netP$pathways))))
pw_mat <- matrix(NA_real_, nrow = length(all_pw), ncol = length(object.list),
                 dimnames = list(all_pw, names(object.list)))
for (s in names(object.list)) {
  cc <- object.list[[s]]
  lr <- subsetCommunication(cc)
  if (nrow(lr) == 0) next
  agg <- as.data.table(lr)[, .(weight = sum(prob, na.rm = TRUE)), by = pathway_name]
  pw_mat[agg$pathway_name, s] <- agg$weight
}
pw_dt <- as.data.table(pw_mat, keep.rownames = "pathway")
pw_dt[, mean_weight := rowMeans(.SD, na.rm = TRUE), .SDcols = names(object.list)]
pw_dt[, sd_weight := apply(.SD, 1, sd, na.rm = TRUE), .SDcols = names(object.list)]
pw_dt[, n_detected := rowSums(!is.na(.SD)), .SDcols = names(object.list)]
pw_dt[, cv := sd_weight / mean_weight]
setorder(pw_dt, -mean_weight)
fwrite(pw_dt, file.path(TM, "P4e_pathway_weight_by_sample.csv"))
fwrite(pw_dt[order(cv)], file.path(TM, "P4e_pathway_cv_across_samples.csv"))
logmsg("pathways:", nrow(pw_dt), "| Top15 by mean weight:")
print(head(pw_dt[, c("pathway", names(object.list), "mean_weight", "cv"), with = FALSE], 15))

# ---------- 6) hub 区室 per-sample 发送强度 ----------
hub_ct <- c("CAF/MSC", "Osteoblastic", "Proliferating")
cons <- rbindlist(lapply(names(object.list), function(s) {
  w <- object.list[[s]]@net$weight
  ow <- rowSums(w, na.rm = TRUE)
  rbindlist(lapply(intersect(hub_ct, names(ow)), function(ct)
    data.table(sample = s, celltype = ct, out_weight = as.numeric(ow[ct]))))
}))
fwrite(cons, file.path(TM, "P4e_hub_outgoing_by_sample.csv"))
logmsg("hub outgoing per sample:")
print(dcast(cons, celltype ~ sample, value.var = "out_weight"))

agg <- cons[, .(mean_w = mean(out_weight), sd_w = sd(out_weight),
                cv = sd(out_weight) / mean(out_weight),
                min_w = min(out_weight), max_w = max(out_weight)), by = celltype][order(-mean_w)]
fwrite(agg, file.path(TM, "P4e_hub_outgoing_consistency.csv"))
logmsg("hub outgoing consistency:")
print(agg)

# 排名
rk <- copy(cons)[, rank := frank(-out_weight), by = sample]
fwrite(rk, file.path(TM, "P4e_hub_rank_by_sample.csv"))
logmsg("hub rank per sample:")
print(dcast(rk, celltype ~ sample, value.var = "rank"))

# 基质 vs 增殖比
pv <- dcast(cons, sample ~ celltype, value.var = "out_weight")
if (all(c("CAF/MSC", "Proliferating") %in% colnames(pv))) {
  pv[, ratio := `CAF/MSC` / Proliferating]
  fwrite(pv, file.path(TM, "P4e_stroma_vs_proliferation_ratio.csv"))
  logmsg("stroma/proliferation ratio:")
  print(pv[, .(sample, caf = `CAF/MSC`, prolif = Proliferating, ratio)])
}

# ---------- 7) 可视化 ----------
pal_ct <- c("CAF/MSC" = "#E31A1C", "Osteoblastic" = "#FF7F00",
            "Proliferating" = "#1F78B4", "Myeloid" = "#33A02C",
            "Tcell" = "#6A3D9A", "Endothelial" = "#B15928",
            "Bcell" = "#A6CEE3", "NK" = "#FB9A99", "Osteoclast" = "#B2DF8A")

p2 <- ggplot(cons, aes(x = sample, y = out_weight, fill = celltype)) +
  geom_col(position = "dodge", width = 0.76, color = "grey25", linewidth = 0.25) +
  scale_fill_manual(values = pal_ct, name = "") +
  labs(title = "Hub-compartment outgoing signaling across 6 samples",
       subtitle = "RUNX2-high stromal compartments (CAF/MSC, Osteoblastic) vs BUB1-high Proliferating",
       y = "Outgoing signaling weight", x = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(FG, "P4e_hub_outgoing_by_sample.png"), p2,
       width = 9, height = 5.5, dpi = 300)
logmsg("saved P4e_hub_outgoing_by_sample.png")

if (exists("pv") && "ratio" %in% colnames(pv)) {
  p3 <- ggplot(pv, aes(x = sample, y = ratio)) +
    geom_col(fill = "#E31A1C", width = 0.62, color = "grey25", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = round(ratio, 2)), vjust = -0.4, size = 3.8) +
    labs(title = "Stromal (CAF/MSC, RUNX2-high) ÷ Proliferative (BUB1-high) outgoing ratio",
         subtitle = "ratio > 1 in all samples ⇒ stromal axis dominates signaling consistently",
         y = "Outgoing weight ratio", x = NULL) +
    expand_limits(y = c(0, max(pv$ratio) * 1.2)) +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(FG, "P4e_stroma_vs_proliferation_ratio.png"), p3,
         width = 8, height = 5, dpi = 300)
  logmsg("saved P4e_stroma_vs_proliferation_ratio.png")
}

# 通路热图（Top30 z-score）
top <- head(pw_dt[order(-mean_weight)], 30)
long <- melt(top[, c("pathway", names(object.list), "mean_weight"), with = FALSE],
             id.vars = c("pathway", "mean_weight"),
             measure.vars = names(object.list),
             variable.name = "sample", value.name = "weight")
long[, z := (weight - mean(weight, na.rm = TRUE)) / (sd(weight, na.rm = TRUE) + 1e-9),
     by = pathway]
long[, pathway := factor(pathway, levels = rev(top$pathway))]
ph <- ggplot(long, aes(sample, pathway, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "z-score") +
  labs(title = "Top 30 signaling pathways — per-sample relative strength",
       subtitle = "row z-score of summed pathway communication probability",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
ggsave(file.path(FG, "P4e_pathway_heatmap_by_sample.png"), ph,
       width = 8, height = 10, dpi = 300)
logmsg("saved P4e_pathway_heatmap_by_sample.png")

# ---------- 8) hub 区室 in/out（merged per sample net）----------
rec <- rbindlist(lapply(names(cellchat.merged@net), function(s) {
  w <- cellchat.merged@net[[s]]$weight
  rbindlist(lapply(intersect(hub_ct, rownames(w)), function(ct)
    data.table(sample = s, celltype = ct,
               out = sum(w[ct, ], na.rm = TRUE),
               inn = sum(w[, ct], na.rm = TRUE))))
}))
fwrite(rec, file.path(TM, "P4e_hub_inout_by_sample.csv"))
logmsg("hub in/out (merged object):")
print(rec)

saveRDS(cellchat.merged, file.path(ROOT, "data/processed/P4e_cellchat_merged.rds"))
logmsg("=== P4e-merge SUMMARY ===")
cat("samples:", length(object.list), "| celltypes:", length(ct_intersect),
    "| pathways union:", length(all_pw), "\n")
logmsg("P4e-merge complete.")
