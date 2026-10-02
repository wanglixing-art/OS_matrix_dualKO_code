# ============================================================
# P7b: ESTIMATE 基质/免疫评分 × 14 基因签名（TARGET-OS 训练集）
#
# 口径：
#   表达 = P1 voom log2CPM 矩阵（13738 基因 × 88 样本）
#   签名 RiskScore = P6 训练系数（P6_signature_genes.csv）
#   样本对齐 = P6 规则（submitter_id "-"→"."），EFS 可评估者 n=85
#   ESTIMATE: filterCommonGenes(id="GeneSymbol") → estimateScore(platform="illumina")
#   统计: 高/低风险组（median split）Wilcoxon + RiskScore 连续 Spearman
# 输出:
#   P7_estimate_scores.csv（每样本 score）
#   P7_estimate_stats.csv（组间比较 + 相关）
#   P7_estimate_boxplots.png（三联箱线图）
# ============================================================
suppressMessages({
  library(data.table); library(estimate); library(ggplot2); library(patchwork)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
LOG <- file.path(PROJ, "logs/P7b_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
TMP <- file.path(PROJ, "data/processed"); dir.create(TMP, showWarnings = FALSE)

# ---------- 1. 表达 + 签名打分（P6 口径） ----------
obj <- readRDS(file.path(PROJ, "data/processed/P1_objects.rds"))
E <- obj$E                       # gene x sample, voom log2CPM
sigdf <- fread(file.path(TM, "P6_signature_genes.csv"))
sig <- sigdf$coef; names(sig) <- sigdf$gene
gvars <- names(sig)[names(sig) %in% rownames(E)]
X <- t(E[gvars, ]); X <- X[, apply(X, 2, sd) > 0]
Xz <- scale(X)
rs <- as.numeric(Xz %*% sig[colnames(X)])
names(rs) <- rownames(X)
logmsg("RiskScore computed on", ncol(X), "genes x", length(rs), "samples")

# ---------- 2. ESTIMATE ----------
# 输入需 GCT： GeneSymbol x sample（estimate 内部取交集common probes）
emat <- data.frame(GeneSymbol = rownames(E), E, check.names = FALSE)
in_f  <- file.path(TMP, "P7_estimate_input.txt")
out_gct <- file.path(TMP, "P7_estimate_input.gct")
out_score <- file.path(TMP, "P7_estimate_score.gct")
fwrite(emat, in_f, sep = "\t", quote = FALSE, row.names = FALSE)
filterCommonGenes(input.f = in_f, output.f = out_gct, id = "GeneSymbol")
estimateScore(input.ds = out_gct, output.ds = out_score, platform = "illumina")
sc <- read.table(out_score, skip = 2, header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)
sc <- as.data.frame(t(sc[, -1]))        # sample x 4 scores
sc$sample <- rownames(sc)
logmsg("ESTIMATE done: ", nrow(sc), "samples; cols:", paste(names(sc), collapse = ","))

# ---------- 3. 合并 + 统计（对齐 P6 训练口径：EFS 可评估 n=85） ----------
clin <- fread(file.path(PROJ, "data/processed/TARGET-OS_clinical.csv"), data.table = FALSE)
clin$id_dot <- gsub("-", ".", clin$submitter_id)
ok85 <- clin$id_dot[!is.na(clin$EFS_time_days) & !is.na(clin$EFS_event) &
                      clin$EFS_time_days > 0]
d <- data.frame(sample = names(rs), RiskScore = rs)
d <- merge(d, sc, by = "sample")
d <- d[d$sample %in% ok85, ]
d$grp <- factor(ifelse(d$RiskScore > median(d$RiskScore), "High", "Low"),
                levels = c("Low", "High"))
logmsg("merged n=", nrow(d), " (EFS-evaluable, P6 training口径)")

svars <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
stat_rows <- data.frame()
for (v in svars) {
  mw <- wilcox.test(as.formula(paste(v, "~ grp")), data = d)
  sp <- suppressWarnings(cor.test(d$RiskScore, d[[v]], method = "spearman"))
  stat_rows <- rbind(stat_rows, data.frame(
    score = v,
    median_low = median(d[[v]][d$grp == "Low"]),
    median_high = median(d[[v]][d$grp == "High"]),
    wilcox_p = signif(mw$p.value, 3),
    spearman_rho = signif(unname(sp$estimate), 3),
    spearman_p = signif(sp$p.value, 3)))
  logmsg(v, ": median(L/H)=", round(median(d[[v]][d$grp == "Low"]), 0), "/",
         round(median(d[[v]][d$grp == "High"]), 0),
         " wilcox p=", signif(mw$p.value, 3),
         " spearman rho=", round(unname(sp$estimate), 3),
         " p=", signif(sp$p.value, 3))
}
fwrite(d, file.path(TM, "P7_estimate_scores.csv"))
fwrite(stat_rows, file.path(TM, "P7_estimate_stats.csv"))

# ---------- 4. 图（三联箱线图，彩色可发表） ----------
pd <- rbind(
  data.frame(grp = d$grp, score = d$StromalScore, panel = "Stromal score"),
  data.frame(grp = d$grp, score = d$ImmuneScore, panel = "Immune score"),
  data.frame(grp = d$grp, score = d$ESTIMATEScore, panel = "ESTIMATE score"))
pd$panel <- factor(pd$panel, levels = c("Stromal score", "Immune score", "ESTIMATE score"))
p <- ggplot(pd, aes(grp, score, fill = grp)) +
  geom_boxplot(width = 0.55, outlier.shape = 21, outlier.size = 1.8,
               color = "grey25", linewidth = 0.45) +
  geom_jitter(width = 0.14, size = 1.5, alpha = 0.5, color = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.8, color = "black") +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#2166AC", "High" = "#D7301F"), name = NULL) +
  labs(title = "Tumor microenvironment scores by LASSO-Cox risk group (TARGET-OS)",
       subtitle = "ESTIMATE on voom-normalized RNA-seq; median-split risk groups",
       x = NULL, y = "ESTIMATE score") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))
ggsave(file.path(FG, "P7_estimate_boxplots.png"), p, width = 8.6, height = 4.4, dpi = 300)
logmsg("=== P7b complete ===")
