# ============================================================
# P4e-post: per-sample CellChat 对比 —— 可视化增强 + 一致性检验
# 依赖 P4e_cellchat_percellchat.R 产出的：
#   data/processed/P4e_cellchat_<sample>.rds
#   data/processed/P4e_cellchat_merged.rds
#   results/tables/P4e_pathway_weight_by_sample.csv
#   results/tables/P4e_hub_outgoing_by_sample.csv
# 产出：
#   - per-sample 通路权重热图（Top 30）
#   - hub 区室发送强度 per-sample 一致性（含 CV）
#   - 基质/增殖比 per-sample 图
#   - 跨样本稳健性统计（hub 区室排名是否稳定）
# ============================================================
suppressMessages({
  library(CellChat); library(Seurat); library(ggplot2)
  library(patchwork); library(data.table); library(Matrix)
})
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
LOG <- file.path(ROOT, "logs/P4e_post.log")
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)

pal_ct <- c("CAF/MSC" = "#E31A1C", "Osteoblastic" = "#FF7F00", "Proliferating" = "#1F78B4",
            "Myeloid" = "#33A02C", "Tcell" = "#6A3D9A", "Endothelial" = "#B15928",
            "Bcell" = "#A6CEE3", "NK" = "#FB9A99", "Osteoclast" = "#B2DF8A")

# ---------- 1) 通路 × 样本 热图 ----------
pw_f <- file.path(TM, "P4e_pathway_weight_by_sample.csv")
if (file.exists(pw_f)) {
  pw <- fread(pw_f)
  samples <- grep("^OS_", colnames(pw), value = TRUE)
  logmsg("pathway table:", nrow(pw), "pathways x", length(samples), "samples")

  top <- head(pw[order(-mean_weight)], 30)
  long <- melt(top, id.vars = c("pathway", "mean_weight", "sd_weight"),
               measure.vars = samples, variable.name = "sample", value.name = "weight")
  # 行归一化（每通路 z-score），凸显样本间相对差异
  long[, z := (weight - mean(weight, na.rm = TRUE)) / (sd(weight, na.rm = TRUE) + 1e-9),
       by = pathway]
  long[, pathway := factor(pathway, levels = rev(top$pathway))]

  p <- ggplot(long, aes(sample, pathway, fill = z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "z-score") +
    labs(title = "Top 30 signaling pathways — per-sample relative strength",
         subtitle = "row z-score of summed pathway communication probability",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 0))
  ggsave(file.path(FG, "P4e_pathway_heatmap_by_sample.png"), p,
         width = 8, height = 10, dpi = 300)
  logmsg("saved P4e_pathway_heatmap_by_sample.png")

  # 通路 CV（跨样本变异系数）—— 低 CV = 稳健
  pw[, cv := sd_weight / mean_weight]
  fwrite(pw[, .(pathway, mean_weight, sd_weight, cv, n_detected = rowSums(!is.na(.SD)))],
         file.path(TM, "P4e_pathway_cv_across_samples.csv"))
  logmsg("most stable pathways (lowest CV):")
  print(head(pw[order(cv), .(pathway, mean_weight, cv)], 10))
}

# ---------- 2) hub 区室发送强度一致性 ----------
hub_f <- file.path(TM, "P4e_hub_outgoing_by_sample.csv")
if (file.exists(hub_f)) {
  cons <- fread(hub_f)
  agg <- cons[, .(mean_w = mean(out_weight), sd_w = sd(out_weight),
                  cv = sd(out_weight) / mean(out_weight),
                  min_w = min(out_weight), max_w = max(out_weight)),
              by = celltype][order(-mean_w)]
  fwrite(agg, file.path(TM, "P4e_hub_outgoing_consistency.csv"))
  logmsg("hub compartment outgoing consistency across 6 samples:")
  print(agg)

  p2 <- ggplot(cons, aes(x = sample, y = out_weight, fill = celltype)) +
    geom_col(position = "dodge", width = 0.76, color = "grey25", linewidth = 0.25) +
    scale_fill_manual(values = pal_ct, name = "") +
    labs(title = "Hub-compartment outgoing signaling across 6 samples",
         y = "Outgoing signaling weight", x = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  ggsave(file.path(FG, "P4e_hub_outgoing_by_sample_v2.png"), p2,
         width = 9, height = 5.5, dpi = 300)

  # 各样本中 hub 区室发送强度排名（验证 RUNX2 轴基质区室是否稳定居首）
  rk <- cons[, .(celltype, out_weight)][order(sample, -out_weight)]
  rk <- cons[, rank := frank(-out_weight), by = sample]
  fwrite(rk, file.path(TM, "P4e_hub_rank_by_sample.csv"))
  logmsg("rank of hub compartments in each sample:")
  print(dcast(rk, celltype ~ sample, value.var = "rank"))

  # 基质(CAF/MSC) vs 增殖(Proliferating) 比值
  pv <- dcast(cons, sample ~ celltype, value.var = "out_weight")
  if (all(c("CAF/MSC", "Proliferating") %in% colnames(pv))) {
    pv[, ratio := `CAF/MSC` / Proliferating]
    fwrite(pv, file.path(TM, "P4e_stroma_vs_proliferation_ratio.csv"))
    logmsg("stroma(CAF/MSC)/proliferation outgoing ratio:")
    print(pv[, .(sample, caf = `CAF/MSC`, prolif = Proliferating, ratio)])

    p3 <- ggplot(pv, aes(x = sample, y = ratio)) +
      geom_col(fill = "#E31A1C", width = 0.65, color = "grey25", linewidth = 0.3) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      geom_text(aes(label = round(ratio, 2)), vjust = -0.4, size = 3.6) +
      labs(title = "Stromal (CAF/MSC, RUNX2-high) vs Proliferative (BUB1-high) outgoing ratio",
           subtitle = "ratio > 1 in all samples ⇒ stromal axis dominates signaling consistently",
           y = "CAF/MSC ÷ Proliferating outgoing weight", x = NULL) +
      theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
    ggsave(file.path(FG, "P4e_stroma_vs_proliferation_ratio.png"), p3,
           width = 8, height = 5, dpi = 300)
    logmsg("saved P4e_stroma_vs_proliferation_ratio.png")
  }
}

# ---------- 3) merged 对象的 per-sample 差异通路（CellChat 官方 rankNet） ----------
mm_f <- file.path(ROOT, "data/processed/P4e_cellchat_merged.rds")
if (file.exists(mm_f)) {
  cc <- readRDS(mm_f)
  logmsg("merged object loaded | groups:", length(cc@net))
  # 每个 hub 区室在各样本中接收/发送的信号
  try({
    hub_ct <- c("CAF/MSC", "Osteoblastic", "Proliferating")
    rec <- list()
    for (ct in hub_ct) {
      for (s in names(cc@net)) {
        w <- cc@net[[s]]$weight
        if (!ct %in% rownames(w)) next
        rec[[length(rec) + 1]] <- data.table(
          celltype = ct, sample = s,
          out = sum(w[ct, ], na.rm = TRUE),
          inn = sum(w[, ct], na.rm = TRUE))
      }
    }
    rd <- rbindlist(rec)
    fwrite(rd, file.path(TM, "P4e_hub_inout_by_sample.csv"))
    logmsg("hub in/out by sample (merged net):")
    print(dcast(rd, celltype + sample ~ ., value.var = c("out"), fun.aggregate = sum)[1:6])
  })
}

logmsg("=== P4e-post complete ===")
