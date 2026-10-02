# ============================================================
# P6: LASSO-Cox EFS 签名构建 + GSE21257 外部验证
#
# 设计（CKAP2 模板方法学 + P2 拍板的 ML 目标 A = EFS 事件预测）：
#   训练: TARGET-OS (n=85, EFS 可用) — voom 表达 + 临床
#   候选: P2 候选池 v3 (299 基因, DEG 上调 ∩ P1 显著模块)
#   流程: 单因素 Cox 预筛(p<0.1) → glmnet LASSO-Cox(10 折 CV)
#   口径: 各数据集内基因 z-score 后线性打分（跨平台可迁移）
#   内评: C-index / timeROC AUC(1-4y) / KM median / 多因素独立性
#   外验: GSE21257 (n=53, OS 终点 + 转移维度) — 不重拟合，用训练系数
# ============================================================
suppressMessages({
  library(data.table); library(survival); library(glmnet); library(timeROC)
  library(GEOquery); library(ggplot2); library(survminer); library(patchwork)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
LOG <- file.path(PROJ, "logs/P6_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
set.seed(2026)

# ---------- Part 1: 数据加载 ----------
logmsg("Part 1: load TARGET-OS expression + clinical")
obj <- readRDS(file.path(PROJ, "data/processed/P1_objects.rds"))
E <- obj$E                                   # voom 归一化表达 (gene x sample)
clin <- fread(file.path(PROJ, "data/processed/TARGET-OS_clinical.csv"), data.table = FALSE)
pool <- fread(file.path(TM, "P2_candidate_pool_v3_matrixaxis.csv"), data.table = FALSE)$gene
logmsg("E:", nrow(E), "genes x", ncol(E), "samples | pool:", length(pool))

# 对齐样本（P1 规则：表达列名 = submitter_id 的 "-" 替换为 "."）
clin$id_dot <- gsub("-", ".", clin$submitter_id)
common <- intersect(colnames(E), clin$id_dot)
clin <- clin[match(common, clin$id_dot), ]; E <- E[, common]
ok <- !is.na(clin$EFS_time_days) & !is.na(clin$EFS_event) & clin$EFS_time_days > 0
clin <- clin[ok, ]; E <- E[, ok]
logmsg("aligned samples:", ncol(E), "| EFS events:", sum(clin$EFS_event),
       "| OS events:", sum(clin$OS_event))

genes <- pool[pool %in% rownames(E)]
logmsg("pool genes in matrix:", length(genes))

# ---------- Part 2: 单因素 Cox 预筛 + LASSO ----------
y_t <- clin$EFS_time_days; y_e <- clin$EFS_event
uni <- t(vapply(genes, function(g) {
  d <- data.frame(t = y_t, e = y_e, x = as.numeric(E[g, ]))
  s <- summary(coxph(Surv(t, e) ~ x, d))
  c(HR = s$coef[2], p = s$coef[5])
}, numeric(2)))
uni <- data.frame(gene = rownames(uni), uni, row.names = NULL)
keep <- uni$p < 0.1 & !is.na(uni$p)
logmsg("univariate Cox p<0.1:", sum(keep), "/", nrow(uni))
sel <- uni$gene[keep]

X <- t(E[sel, ])                             # sample x gene
X <- X[, apply(X, 2, sd, na.rm = TRUE) > 0]  # 防 sd=0 基因
Xz <- scale(X)                               # 数据集内 z-score
cv.fit <- cv.glmnet(Xz, Surv(y_t, y_e), family = "cox", alpha = 1,
                    nfolds = 10, type.measure = "deviance")
n_1se  <- sum(coef(cv.fit, s = "lambda.1se")  != 0)
n_min  <- sum(coef(cv.fit, s = "lambda.min")  != 0)
logmsg("LASSO: lambda.1se keeps ", n_1se, " genes | lambda.min keeps ", n_min,
       " genes | lambda.min=", round(cv.fit$lambda.min, 5),
       " lambda.1se=", round(cv.fit$lambda.1se, 5))
lam <- if (n_1se >= 5) "lambda.1se" else "lambda.min"
beta <- as.numeric(coef(cv.fit, s = lam))
names(beta) <- rownames(coef(cv.fit, s = lam))
sig <- beta[beta != 0]
sig <- sig[order(-abs(sig))]
logmsg("signature (", lam, "): ", length(sig), " genes")
print(round(sig, 4))
fwrite(data.frame(gene = names(sig), coef = as.numeric(sig)),
       file.path(TM, "P6_signature_genes.csv"))
fwrite(data.frame(gene = uni$gene, HR = uni$HR, p_uni = uni$p, prescreen = keep),
       file.path(TM, "P6_univariate_screen.csv"))

rs_train <- as.numeric(Xz[, names(sig)] %*% sig)
score_train <- data.frame(submitter_id = colnames(E), RiskScore = rs_train,
                          EFS_time = y_t, EFS_event = y_e,
                          OS_time = clin$OS_time_days, OS_event = clin$OS_event,
                          age = clin$age_at_diagnosis_years, sex = clin$gender)
fwrite(score_train, file.path(TM, "P6_riskscore_TARGET_train.csv"))

# ---------- Part 3: 训练内部评估 ----------
perf_row <- function(time, event, score) {
  s <- summary(coxph(Surv(time, event) ~ score))
  c(C_index = s$concordance["C"], se = s$concordance["sestd.C"],
    HR = s$coef[2], HR_lo = s$conf.int[3], HR_hi = s$conf.int[4], p = s$coef[5])
}
perf_train <- perf_row(y_t, y_e, rs_train)
logmsg("TRAIN EFS: C=", round(perf_train[1], 3), " HR=",
       round(perf_train[3], 2), " [", round(perf_train[4], 2), "-",
       round(perf_train[5], 2), "] p=", signif(perf_train[6], 3))
perf_train_os <- perf_row(clin$OS_time_days, clin$OS_event, rs_train)
logmsg("TRAIN OS : C=", round(perf_train_os[1], 3), " p=", signif(perf_train_os[6], 3))

# KM median split
km_eval <- function(time, event, score, cut = NULL) {
  grp <- if (is.null(cut)) factor(ifelse(score > median(score), "High", "Low"),
                                  levels = c("Low", "High")) else
         factor(ifelse(score > cut, "High", "Low"), levels = c("Low", "High"))
  sd <- survdiff(Surv(time, event) ~ grp)
  p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  list(grp = grp, p_logrank = p)
}
k1 <- km_eval(y_t, y_e, rs_train)
logmsg("TRAIN KM logrank p=", signif(k1$p_logrank, 3))

# 时间依赖 AUC（1/2/3/4 年）
auc_tab <- data.frame()
for (yr in c(1, 2, 3, 4)) {
  tt <- yr * 365.25
  if (sum(y_t >= tt) < 10) next
  r <- timeROC(T = y_t, delta = y_e, marker = rs_train, cause = 1,
               times = tt, ROC = TRUE)
  auc_tab <- rbind(auc_tab, data.frame(year = yr, AUC = r$AUC[2]))
}
logmsg("TRAIN time-dependent AUC:"); print(auc_tab)
fwrite(auc_tab, file.path(TM, "P6_timeAUC_train.csv"))

# 多因素独立性（score + age + sex）
mv <- coxph(Surv(y_t, y_e) ~ rs_train + age + sex, score_train)
mv_tab <- data.frame(term = names(coef(mv)), summary(mv)$conf.int[, c(1, 3, 4)],
                     p = summary(mv)$coef[, 5])
logmsg("multivariable (train EFS):"); print(mv_tab)
fwrite(mv_tab, file.path(TM, "P6_multivariable_train.csv"))
fwrite(data.frame(metric = c(names(perf_train), paste0("AUC_", auc_tab$year, "y_train")),
                  value = c(perf_train, auc_tab$AUC)),
       file.path(TM, "P6_train_performance.csv"))

# ---------- Part 4: GSE21257 外部验证 ----------
logmsg("Part 4: GSE21257 external validation")
gse <- getGEO(filename = file.path(PROJ, "data/raw/GSE21257/series_matrix.txt.gz"),
              getGPL = FALSE)
ex <- exprs(gse); pd <- pData(gse)
plat <- fread(file.path(PROJ, "data/raw/GSE42352/GPL10295_probe_symbol.tsv"),
              data.table = FALSE, select = c("ID", "Symbol"))
plat <- plat[plat$Symbol != "" & !is.na(plat$Symbol), ]
plat <- plat[!duplicated(plat$ID), ]
ex <- ex[rownames(ex) %in% plat$ID, ]; plat <- plat[match(rownames(ex), plat$ID), ]
mmean <- rowMeans(ex); ord <- order(plat$Symbol, -mmean)
ex2 <- ex[ord, ]; pl2 <- plat[ord, ]; first <- !duplicated(pl2$Symbol)
ex2 <- ex2[first, ]; rownames(ex2) <- pl2$Symbol[first]
logmsg("GSE21257 mapped:", nrow(ex2), "genes x", ncol(ex2), "samples")

miss <- setdiff(names(sig), rownames(ex2))
logmsg("signature genes missing in GSE21257:", length(miss),
       if (length(miss) > 0) paste0("(", paste(miss, collapse = ","), ")") else "")

# 解析临床（P3 同款解析器）
chars <- lapply(grep("^characteristics_ch1", colnames(pd), value = TRUE),
                function(cn) pd[[cn]])
getfield <- function(key) {
  out <- rep(NA_character_, nrow(pd))
  for (v in chars) {
    hit <- grepl(paste0("^", key, ":"), v, ignore.case = TRUE)
    out[hit] <- sub(paste0("^", key, ":\\s*"), "", v[hit], ignore.case = TRUE)
  }
  out
}
st <- getfield("status")
os_time <- as.numeric(gsub("[^0-9]", "", sub(".*?(\\d+)\\s*months.*", "\\1", st)))
os_time[!grepl("month", st, ignore.case = TRUE)] <- NA
os_time <- os_time * 30.44                                # 月 -> 天（与 TARGET 同单位）
os_event <- ifelse(grepl("^Deceased", st, ignore.case = TRUE), 1,
             ifelse(grepl("^Alive", st, ignore.case = TRUE), 0, NA))
grp <- getfield("group")
met <- ifelse(grepl("present at diagnosis", grp, ignore.case = TRUE), "Met_at_dx",
        ifelse(grepl("^Metastases", grp, ignore.case = TRUE), "Met_later",
        ifelse(grepl("No metastases|Non-metastatic", grp, ignore.case = TRUE), "No_met", NA)))

gvars <- names(sig)[names(sig) %in% rownames(ex2)]
Xg <- t(ex2[gvars, colnames(ex2), drop = FALSE])
Xg <- Xg[, apply(Xg, 2, sd) > 0]
Xgz <- scale(Xg)
rs_gse <- as.numeric(Xgz[, names(sig)[names(sig) %in% gvars]] %*%
                       sig[names(sig)[names(sig) %in% gvars]])
score_gse <- data.frame(sample = colnames(ex2), RiskScore = rs_gse,
                        OS_time = os_time, OS_event = os_event, met = met)
score_gse <- score_gse[!is.na(score_gse$OS_time) & !is.na(score_gse$OS_event) &
                         score_gse$OS_time > 0, ]
rs_gse2 <- score_gse$RiskScore
logmsg("GSE21257 evaluable: n=", nrow(score_gse), " events=", sum(score_gse$OS_event))

perf_ext <- perf_row(score_gse$OS_time, score_gse$OS_event, rs_gse2)
logmsg("EXT OS: C=", round(perf_ext[1], 3), " HR=", round(perf_ext[3], 2),
       " [", round(perf_ext[4], 2), "-", round(perf_ext[5], 2), "] p=",
       signif(perf_ext[6], 3))
k2 <- km_eval(score_gse$OS_time, score_gse$OS_event, rs_gse2)
logmsg("EXT KM logrank p=", signif(k2$p_logrank, 3))

# 转移维度（Met_later + Met_at_dx 合并 vs No_met）
score_gse$met_any <- ifelse(score_gse$met %in% c("Met_later", "Met_at_dx"), "Met", "No_met")
md <- score_gse[!is.na(score_gse$met_any), ]
mw2 <- wilcox.test(RiskScore ~ met_any, data = md)
logmsg("EXT metastasis (Met vs No_met): median ",
       round(median(md$RiskScore[md$met_any == "Met"]), 3),
       " vs ",
       round(median(md$RiskScore[md$met_any == "No_met"]), 3),
       " | Wilcoxon p=", signif(mw2$p.value, 3))

# 外部时间 AUC（2/3/5 年）
auc_ext <- data.frame()
for (yr in c(2, 3, 5)) {
  tt <- yr * 365.25
  if (sum(score_gse$OS_time >= tt, na.rm = TRUE) < 8) next
  r <- timeROC(T = score_gse$OS_time, delta = score_gse$OS_event,
               marker = rs_gse2, cause = 1, times = tt, ROC = TRUE)
  auc_ext <- rbind(auc_ext, data.frame(year = yr, AUC = r$AUC[2]))
}
logmsg("EXT time AUC:"); print(auc_ext)
fwrite(score_gse, file.path(TM, "P6_riskscore_GSE21257.csv"))
if (nrow(auc_ext) > 0) {
  fwrite(data.frame(metric = c("C_EXT_OS", names(perf_ext)[-1],
                               "logrank_EXT", paste0("AUC_", auc_ext$year, "y_EXT"),
                               "wilcox_met_EXT"),
                    value = c(perf_ext[1], perf_ext[-1], k2$p_logrank, auc_ext$AUC,
                              mw2$p.value)),
         file.path(TM, "P6_external_performance.csv"))
} else {
  fwrite(data.frame(metric = c("C_EXT_OS", names(perf_ext)[-1], "logrank_EXT",
                               "wilcox_met_EXT"),
                    value = c(perf_ext[1], perf_ext[-1], k2$p_logrank, mw2$p.value)),
         file.path(TM, "P6_external_performance.csv"))
}

# ---------- Part 5: KM 图（训练 + 外部）----------
plot_km <- function(time, event, grp, ttl) {
  fit <- survfit(Surv(time, event) ~ grp)
  p <- ggsurvplot(fit, data = data.frame(time, event, grp),
                  pval = paste0("logrank p = ", signif(km_eval(time, event, as.numeric(grp))$p_logrank, 3)),
                  risk.table = TRUE, conf.int = TRUE,
                  palette = c("#2166AC", "#D7301F"),
                  legend.labs = c("Low risk", "High risk"),
                  title = ttl, xlab = "Time (days)",
                  font.title = 12, font.subtitle = 10, surv.median.line = "hv")
  p
}
png(file.path(FG, "P6_KM_train.png"), 1500, 1300, res = 150)
print(plot_km(y_t, y_e, k1$grp,
              paste0("TARGET-OS training (n=", length(y_t), "), EFS endpoint")))
dev.off()
png(file.path(FG, "P6_KM_GSE21257.png"), 1500, 1300, res = 150)
print(plot_km(score_gse$OS_time, score_gse$OS_event, k2$grp,
              paste0("GSE21257 external (n=", nrow(score_gse), "), OS endpoint")))
dev.off()

# 时间 AUC 曲线图（外部 AUC 可能为空，防御处理）
auc_all <- auc_tab
if (nrow(auc_ext) > 0) {
  auc_all <- rbind(cbind(auc_tab, cohort = "Training (TARGET-OS, EFS)"),
                   cbind(auc_ext, cohort = "External (GSE21257, OS)"))
} else {
  auc_all <- cbind(auc_tab, cohort = "Training (TARGET-OS, EFS)")
}
pa <- ggplot(auc_all, aes(factor(year), AUC, fill = cohort, group = cohort)) +
  geom_line(position = position_dodge(0.2), color = "grey55", linewidth = 0.5) +
  geom_point(size = 3.2, shape = 21, color = "grey20", position = position_dodge(0.2)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") +
  scale_fill_manual(values = c("Training (TARGET-OS, EFS)" = "#E31A1C",
                               "External (GSE21257, OS)" = "#1F78B4"), name = "") +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(title = "Time-dependent AUC of the LASSO-Cox risk score",
       x = "Time (years)", y = "AUC") +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(FG, "P6_timeAUC.png"), pa, width = 7, height = 5.2, dpi = 300)

# 签名系数条形图
sg <- data.frame(gene = names(sig), coef = as.numeric(sig))
sg <- sg[order(sg$coef), ]
sg$gene <- factor(sg$gene, levels = sg$gene)
ps <- ggplot(sg, aes(gene, coef, fill = coef > 0)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.25) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#2166AC"),
                    labels = c(`TRUE` = "Higher risk", `FALSE` = "Lower risk"),
                    name = "") +
  labs(title = paste0("LASSO-Cox signature (", nrow(sg), " genes, ", lam, ")"),
       subtitle = "Coefficients from TARGET-OS training (EFS)",
       x = NULL, y = "LASSO coefficient") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(FG, "P6_signature_coefs.png"), ps, width = 7.5, height = 5.8, dpi = 300)

logmsg("=== P6 complete ===")
cat("\nSIGNATURE (", length(sig), " genes):\n"); print(round(sig, 4))
cat("\nTrain C (EFS):", round(perf_train[1], 3), "| Ext C (OS):", round(perf_ext[1], 3), "\n")
