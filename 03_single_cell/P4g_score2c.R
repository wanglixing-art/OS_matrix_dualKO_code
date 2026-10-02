# ============================================================
# P4g-score2c: 补第 8 步 —— hub 表达 vs inferCNV 恶性（score2 因热图崩溃未执行）
# ============================================================
suppressMessages({ library(data.table); library(ggplot2); library(Seurat) })
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
LOG <- file.path(ROOT, "logs/P4g_score2c.log")
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)

res <- fread(file.path(TM, "P4g_percell_cnv_score.csv"))
logmsg("res rows:", nrow(res))

seu <- readRDS(file.path(ROOT, "data/processed/P4_seurat_annotated.rds"))
hb <- FetchData(seu, vars = c("BUB1", "RUNX2"), layer = "data")
hb$cell <- rownames(hb)
hbt <- as.data.table(hb)
mg2 <- merge(hbt, res[, .(cell, malignant_95, cnv_score, cnv_z, celltype, group)], by = "cell")
logmsg("merged rows:", nrow(mg2))
fwrite(mg2, file.path(TM, "P4g_hub_expr_vs_infercnv_malignant.csv"))

# Wilcoxon 检验（obs 内，Malignant vs Non-malignant）
obs <- mg2[group == "obs"]
stat <- obs[, .(
  hub = c("BUB1", "RUNX2"),
  med_mal = c(median(BUB1[malignant_95 == "Malignant"]), median(RUNX2[malignant_95 == "Malignant"])),
  med_nonmal = c(median(BUB1[malignant_95 == "Non-malignant"]), median(RUNX2[malignant_95 == "Non-malignant"])),
  p = c(wilcox.test(BUB1 ~ malignant_95, obs)$p.value,
        wilcox.test(RUNX2 ~ malignant_95, obs)$p.value)
)]
fwrite(stat, file.path(TM, "P4g_hub_wilcox_vs_malignant.csv"))
logmsg("hub Wilcoxon (obs):"); print(stat)

long <- melt(mg2, id.vars = c("cell", "malignant_95", "group", "celltype"),
             measure.vars = c("BUB1", "RUNX2"), variable.name = "hub", value.name = "expr")
ph <- ggplot(long[group == "obs"], aes(malignant_95, expr, fill = malignant_95)) +
  geom_violin(scale = "width", alpha = 0.9, linewidth = 0.3, color = "grey25") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", linewidth = 0.3) +
  facet_wrap(~hub, scales = "free_y") +
  scale_fill_manual(values = c("Malignant" = "#D7263D", "Non-malignant" = "#2E86AB"), name = "") +
  labs(title = "Hub expression vs inferCNV malignancy (obs cells)",
       x = NULL, y = "log-normalized expression") +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(FG, "P4g_hub_vs_infercnv_malignant.png"), ph,
       width = 8, height = 4.8, dpi = 300)
logmsg("saved P4g_hub_vs_infercnv_malignant.png")
logmsg("P4g-score2c complete.")
