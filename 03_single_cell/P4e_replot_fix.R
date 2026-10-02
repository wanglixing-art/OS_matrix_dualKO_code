# ============================================================
# P4e-replot: 修复两张图的 Unicode 字符 bug + NA 显式标注
# 修复点:
#   1. "÷" -> "/"  (Windows PNG 设备缺字形显示为 C7)
#   2. "⇒" -> ":"  (显示为 b..)
#   3. "—" -> "-"  (被截断)
#   4. 通路热图 NA -> 浅灰 + caption "grey = not detected"
# 数据源: 已保存的 results/tables/*.csv (无重算)
# ============================================================
suppressMessages({ library(data.table); library(ggplot2) })
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")

pal_ct <- c("CAF/MSC" = "#E31A1C", "Osteoblastic" = "#FF7F00", "Proliferating" = "#1F78B4")

# ---------- 图1: hub outgoing ----------
cons <- fread(file.path(TM, "P4e_hub_outgoing_by_sample.csv"))
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
cat("saved hub_outgoing\n")

# ---------- 图2: stroma/proliferation ratio (ASCII only) ----------
pv <- fread(file.path(TM, "P4e_stroma_vs_proliferation_ratio.csv"))
p3 <- ggplot(pv, aes(x = sample, y = ratio)) +
  geom_col(fill = "#E31A1C", width = 0.62, color = "grey25", linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = round(ratio, 2)), vjust = -0.4, size = 3.8) +
  labs(title = "Stromal (CAF/MSC, RUNX2-high) vs Proliferative (BUB1-high) outgoing ratio",
       subtitle = "Ratio > 1 in all 6 samples: stromal axis dominates signaling consistently",
       y = "Outgoing weight ratio (CAF/MSC / Proliferating)", x = NULL) +
  expand_limits(y = c(0, max(pv$ratio) * 1.2)) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FG, "P4e_stroma_vs_proliferation_ratio.png"), p3,
       width = 8, height = 5, dpi = 300)
cat("saved stroma_vs_prolif\n")

# ---------- 图3: pathway heatmap (NA = grey85 + caption) ----------
pw <- fread(file.path(TM, "P4e_pathway_weight_by_sample.csv"))
samples <- grep("^OS_", colnames(pw), value = TRUE)
top <- head(pw[order(-mean_weight)], 30)
long <- melt(top, id.vars = c("pathway", "mean_weight", "sd_weight"),
             measure.vars = samples, variable.name = "sample", value.name = "weight")
long[, z := (weight - mean(weight, na.rm = TRUE)) / (sd(weight, na.rm = TRUE) + 1e-9),
     by = pathway]
long[, pathway := factor(pathway, levels = rev(top$pathway))]
ph <- ggplot(long, aes(sample, pathway, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "z-score", na.value = "grey85") +
  labs(title = "Top 30 signaling pathways: per-sample relative strength",
       subtitle = "Row z-score of summed pathway communication probability",
       caption = "Grey = pathway not detected in that sample (NA)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.caption = element_text(color = "grey40", size = 9, hjust = 0))
ggsave(file.path(FG, "P4e_pathway_heatmap_by_sample.png"), ph,
       width = 8, height = 10, dpi = 300)
cat("saved pathway_heatmap\n")

cat("=== P4e-replot done ===\n")
