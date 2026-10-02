# ============================================================
# P4d-replot: 从已保存的 P4_cellchat.rds 重绘失败/空白的图
# 无需重跑推断（57 分钟置换检验结果已在 checkpoint 中）
# ============================================================
suppressMessages({ library(CellChat); library(ggplot2); library(data.table) })
PROJ <- "D:/projects/OS_matrix_dualKO"
FG <- file.path(PROJ, "results/figures"); TM <- file.path(PROJ, "results/tables")
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

logmsg("Loading P4_cellchat.rds ...")
cellchat <- readRDS(file.path(PROJ, "data/processed/P4_cellchat.rds"))

# centrality（roles scatter 依赖；纯矩阵运算，秒级）
cellchat <- netAnalysis_computeCentrality(cellchat)
logmsg("centrality computed")

# ---------- heatmap：直接调用（pheatmap/ggplot 都会画到当前设备） ----------
logmsg("heatmap count ...")
png(file.path(FG, "P4_cellchat_heatmap_count.png"), 1500, 1300, res = 150)
r1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
if (!is.null(r1)) print(r1)   # pheatmap S4 对象，print/show 即绘制
dev.off()

logmsg("heatmap weight ...")
png(file.path(FG, "P4_cellchat_heatmap_weight.png"), 1500, 1300, res = 150)
r2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
if (!is.null(r2)) print(r2)
dev.off()

# ---------- 信号角色 scatter（自绘：out vs in weight，避免 ggrepel 兼容问题） ----------
logmsg("roles scatter (custom ggplot) ...")
roles <- fread(file.path(TM, "P4_cellchat_signaling_roles.csv"))
hub_ct <- c("CAF/MSC","Osteoblastic","Proliferating")
roles[, hub := ifelse(celltype %in% hub_ct, "Hub compartment", "Other")]
p <- ggplot(roles, aes(out_weight, in_weight, color = hub, label = celltype)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(vjust = -0.9, size = 3.6, show.legend = FALSE) +
  scale_color_manual(values = c("Hub compartment" = "#E31A1C", "Other" = "#6B7A8F")) +
  geom_hline(yintercept = median(roles$in_weight), linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = median(roles$out_weight), linetype = "dashed", color = "grey60") +
  labs(x = "Outgoing signaling weight (sender)",
       y = "Incoming signaling weight (receiver)",
       title = "OS scRNA-seq: signaling sender vs receiver roles",
       color = "") +
  theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))
png(file.path(FG, "P4_cellchat_roles_scatter.png"), 1500, 1200, res = 150)
print(p); dev.off()

# ---------- 2D 视觉汇总（senders/receivers 排名条形图，审稿常用） ----------
logmsg("sender/receiver rank bars (custom ggplot) ...")
r1d <- melt(roles, id.vars = c("celltype","hub"),
            measure.vars = c("out_weight","in_weight"),
            variable.name = "role", value.name = "weight")
r1d[, role := ifelse(role == "out_weight", "Sender (outgoing)", "Receiver (incoming)")]
r1d[, celltype := factor(celltype, levels = roles[order(-out_weight), celltype])]
p2 <- ggplot(r1d, aes(celltype, weight, fill = role)) +
  geom_col(position = "dodge", width = 0.75, alpha = 0.92) +
  scale_fill_manual(values = c("Sender (outgoing)" = "#E31A1C",
                               "Receiver (incoming)" = "#1F78B4")) +
  labs(x = NULL, y = "Aggregate signaling weight",
       title = "Senders vs receivers across cell types", fill = "") +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                                   plot.title = element_text(hjust = 0.5))
png(file.path(FG, "P4_cellchat_rank_weight.png"), 1600, 1100, res = 150)
print(p2); dev.off()

logmsg("replot done.")
ls_out <- list.files(FG, pattern = "P4_cellchat")
logmsg("figures now:", paste(ls_out, collapse = ", "))
