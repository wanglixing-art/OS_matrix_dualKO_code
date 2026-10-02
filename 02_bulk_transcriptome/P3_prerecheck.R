# ============================================================
# P3: 预检索决策点（半天，可回滚）
# 双 hub = BUB1（周期轴） vs RUNX2（成骨轴）
# 内容: TARGET-OS 表达分层 KM(OS+EFS) + Cox 单变量 + 临床关联 + GSE21257 重现性
# ============================================================
suppressMessages({
  library(data.table); library(GEOquery); library(survival); library(survminer)
  library(ggplot2); library(ggpubr)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
HUBS <- c("BUB1", "RUNX2")
set.seed(20260928)

# ---------- Part 1: TARGET-OS 数据准备 ----------
logmsg("Part 1: TARGET-OS data ...")
E <- read.csv(file.path(PROJ, "data/processed/TARGET-OS_voom_E.csv"), row.names = 1, check.names = FALSE)
clin <- fread(file.path(PROJ, "data/processed/TARGET-OS_clinical.csv"), data.table = FALSE)
rownames(clin) <- gsub("-", ".", clin$submitter_id)
common <- intersect(colnames(E), rownames(clin))
E <- E[, common]; clin <- clin[common, ]
logmsg("samples:", ncol(E), "| OS usable:", sum(!is.na(clin$OS_time_days)),
       "| EFS usable:", sum(!is.na(clin$EFS_time_days)))

km_plot <- function(expr, time, event, gene, endpoint, fout) {
  d <- data.frame(time = time, event = event, expr = as.numeric(expr))
  d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$expr) & d$time > 0, ]
  d$group <- factor(ifelse(d$expr > median(d$expr), "High", "Low"), levels = c("Low", "High"))
  n_h <- sum(d$group == "High"); n_l <- sum(d$group == "Low")
  ev_h <- sum(d$event[d$group == "High"] == 1); ev_l <- sum(d$event[d$group == "Low"] == 1)
  fit <- survfit(Surv(time, event) ~ group, data = d)
  sd <- survdiff(Surv(time, event) ~ group, data = d)
  p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  hr <- summary(coxph(Surv(time, event) ~ scale(as.numeric(expr)), data = d))
  png(fout, 1400, 1200, res = 150)
  print(ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                   palette = c("#2166AC", "#D7301F"), legend.labs = c(paste0("Low (n=", n_l, ")"), paste0("High (n=", n_h, ")")),
                   title = paste0(gene, " — ", endpoint, " (TARGET-OS, median split)"),
                   xlab = "Days", ylab = paste0(endpoint, " probability"), risk.table.height = 0.25))
  dev.off()
  data.frame(gene = gene, endpoint = endpoint, n = nrow(d), n_high = n_h, n_low = n_l,
             events_high = ev_h, events_low = ev_l, logrank_p = p,
             HR_per_SD = hr$coef[2], HR_p = hr$conf.int[5], CI_low = hr$conf.int[3], CI_high = hr$conf.int[4])
}

# ---------- Part 2: KM 生存分析（双 hub × 双终点）----------
logmsg("Part 2: KM analysis ...")
km_res <- data.frame()
for (g in HUBS) {
  if (!g %in% rownames(E)) { logmsg("  ", g, " missing in TARGET"); next }
  km_res <- rbind(km_res, km_plot(E[g, ], clin$OS_time_days, clin$OS_event, g, "OS",
                                  file.path(FG, paste0("P3_KM_", g, "_OS.png"))))
  km_res <- rbind(km_res, km_plot(E[g, ], clin$EFS_time_days, clin$EFS_event, g, "EFS",
                                  file.path(FG, paste0("P3_KM_", g, "_EFS.png"))))
}
print(km_res)
fwrite(km_res, file.path(TM, "P3_KM_TARGET.csv"))

# ---------- Part 3: Cox 单变量（hub + 临床协变量）----------
logmsg("Part 3: univariate Cox ...")
cox_one <- function(var, time, event, label, type) {
  d <- data.frame(t = time, e = event, v = as.numeric(var))
  d <- d[!is.na(d$t) & !is.na(d$e) & !is.na(d$v) & d$t > 0, ]
  if (nrow(d) < 10 || length(unique(d$e)) < 2) return(NULL)
  if (type == "cont") d$v <- as.numeric(scale(d$v))
  s <- summary(coxph(Surv(t, e) ~ v, data = d))
  data.frame(variable = label, n = nrow(d), events = sum(d$e == 1),
             HR = s$coef[2], CI_low = s$conf.int[3], CI_high = s$conf.int[4], p = s$coef[5])
}
cox_tab <- data.frame()
for (g in HUBS) {
  if (!g %in% rownames(E)) next
  cox_tab <- rbind(cox_tab, cbind(endpoint = "OS", cox_one(E[g, ], clin$OS_time_days, clin$OS_event, g, "cont")))
  cox_tab <- rbind(cox_tab, cbind(endpoint = "EFS", cox_one(E[g, ], clin$EFS_time_days, clin$EFS_event, g, "cont")))
}
cox_tab <- rbind(cox_tab,
  cbind(endpoint = "OS", cox_one(as.numeric(clin$age_at_diagnosis_years), clin$OS_time_days, clin$OS_event, "Age (per SD)", "cont")),
  cbind(endpoint = "EFS", cox_one(as.numeric(clin$age_at_diagnosis_years), clin$EFS_time_days, clin$EFS_event, "Age (per SD)", "cont")),
  cbind(endpoint = "OS", cox_one(ifelse(clin$gender == "male", 1, 0), clin$OS_time_days, clin$OS_event, "Male vs Female", "bin")),
  cbind(endpoint = "EFS", cox_one(ifelse(clin$gender == "male", 1, 0), clin$EFS_time_days, clin$EFS_event, "Male vs Female", "bin")))
print(cox_tab, row.names = FALSE)
fwrite(cox_tab, file.path(TM, "P3_Cox_univariate_TARGET.csv"))

# ---------- Part 4: 临床关联（Huvos/转移 不可用于 TARGET，改做表达-年龄/性别）----------
logmsg("Part 4: expression vs clinical ...")
cor_age <- sapply(HUBS, function(g) {
  if (!g %in% rownames(E)) return(c(r = NA, p = NA))
  ct <- cor.test(as.numeric(E[g, ]), as.numeric(clin$age_at_diagnosis_years), method = "spearman")
  c(r = ct$estimate, p = ct$p.value)
})
wt <- lapply(HUBS, function(g) {
  if (!g %in% rownames(E)) return(NULL)
  d <- data.frame(expr = as.numeric(E[g, ]), sex = ifelse(clin$gender == "male", "Male", "Female"))
  d <- d[!is.na(d$expr), ]
  data.frame(gene = g, wilcox_p = wilcox.test(expr ~ sex, d)$p.value)
})
fwrite(data.frame(gene = HUBS, age_spearman_r = cor_age[1, ], age_p = cor_age[2, ]),
       file.path(TM, "P3_hub_expr_clinical_TARGET.csv"))
logmsg("age correlation:"); print(cor_age)

# ---------- Part 5: GSE21257 外部重现性（同平台 GPL10295）----------
logmsg("Part 5: GSE21257 external replication ...")
gse <- getGEO(filename = file.path(PROJ, "data/raw/GSE21257/series_matrix.txt.gz"), getGPL = FALSE)
ex <- exprs(gse); pd <- pData(gse)
# 探针 -> symbol
plat <- fread(file.path(PROJ, "data/raw/GSE42352/GPL10295_probe_symbol.tsv"), data.table = FALSE,
              select = c("ID", "Symbol"))
plat <- plat[plat$Symbol != "" & !is.na(plat$Symbol), ]; plat <- plat[!duplicated(plat$ID), ]
ex <- ex[rownames(ex) %in% plat$ID, ]; plat <- plat[match(rownames(ex), plat$ID), ]
mmean <- rowMeans(ex); ord <- order(plat$Symbol, -mmean)
ex2 <- ex[ord, ]; pl2 <- plat[ord, ]; first <- !duplicated(pl2$Symbol)
ex2 <- ex2[first, ]; rownames(ex2) <- pl2$Symbol[first]
logmsg("GSE21257 mapped:", nrow(ex2), "genes x", ncol(ex2), "samples")

# 解析临床：status 字符串 (Deceased at N months / Alive ...)
chars <- lapply(grep("^characteristics_ch1", colnames(pd), value = TRUE), function(cn) pd[[cn]])
names(chars) <- paste0("ch", seq_along(chars))
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
os_event <- ifelse(grepl("^Deceased", st, ignore.case = TRUE), 1,
             ifelse(grepl("^Alive", st, ignore.case = TRUE), 0, NA))
huvos <- as.numeric(gsub("[^0-9]", "", getfield("huvos grade")))
grp <- getfield("group")
met <- ifelse(grepl("present at diagnosis", grp, ignore.case = TRUE), "Met_at_dx",
        ifelse(grepl("^Metastases", grp, ignore.case = TRUE), "Met_later",
        ifelse(grepl("No metastases|Non-metastatic", grp, ignore.case = TRUE), "No_met", NA)))
age_m <- as.numeric(gsub("[^0-9]", "", getfield("age")))
sex <- getfield("gender")
logmsg("GSE21257 clinical: OS time parsable:", sum(!is.na(os_time)),
       "| events:", sum(os_event == 1, na.rm = TRUE),
       "| met group:", paste(table(met), collapse = "/"),
       "| Huvos parsable:", sum(!is.na(huvos)))

# 5.1 hub 在肿瘤中的表达与生存
rep_res <- data.frame()
for (g in HUBS) {
  if (!g %in% rownames(ex2)) { logmsg("  ", g, " missing in GSE21257"); next }
  v <- as.numeric(ex2[g, ])
  d <- data.frame(time = os_time, event = os_event, expr = v)
  d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$expr) & d$time > 0, ]
  d$sexpr <- as.numeric(scale(d$expr))
  s <- summary(coxph(Surv(time, event) ~ sexpr, d))
  d$group <- factor(ifelse(d$expr > median(d$expr), "High", "Low"), levels = c("Low", "High"))
  sd <- survdiff(Surv(time, event) ~ group, d)
  p_lr <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  fout <- file.path(FG, paste0("P3_KM_", g, "_OS_GSE21257.png"))
  fit <- survfit(Surv(time, event) ~ group, d)
  png(fout, 1400, 1200, res = 150)
  print(ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                   palette = c("#2166AC", "#D7301F"), title = paste0(g, " — OS (GSE21257, n=", nrow(d), ")"),
                   xlab = "Months", risk.table.height = 0.25))
  dev.off()
rep_row <- data.frame(gene = g, n = nrow(d), events = sum(d$event == 1),
                      logrank_p = p_lr, HR_per_SD = s$coef[2], HR_p = s$conf.int[5],
                      CI_low = s$conf.int[3], CI_high = s$conf.int[4],
                      met_wilcox_p = NA_real_, met_early_vs_late_p = NA_real_,
                      huvos_cor_p = NA_real_)
  # 转移/化疗反应维度的关联（P3 通过标准的另一维度）
  if (!all(is.na(met))) {
    dm <- data.frame(expr = v, met = met)
    dm <- dm[!is.na(dm$met) & !is.na(dm$expr), ]
    dm$met_bin <- ifelse(dm$met == "No_met", "No", "Yes")   # 有转移(含later) vs 无
    if (length(unique(dm$met_bin)) == 2) {
      w <- wilcox.test(expr ~ met_bin, dm)
      rep_row$met_wilcox_p <- w$p.value
      dm3 <- dm[dm$met != "No_met", ]
      if (length(unique(dm3$met)) == 2) rep_row$met_early_vs_late_p <- wilcox.test(expr ~ met, dm3)$p.value
    }
  }
  if (!all(is.na(huvos))) {
    dh <- data.frame(expr = v, huvos = huvos)
    dh <- dh[!is.na(dh$huvos) & !is.na(dh$expr), ]
    if (length(unique(dh$huvos)) > 1) rep_row$huvos_cor_p <- cor.test(dh$expr, dh$huvos, method = "spearman")$p.value
  }
  rep_res <- rbind(rep_res, rep_row)
}
print(rep_res, row.names = FALSE)
fwrite(rep_res, file.path(TM, "P3_external_replication_GSE21257.csv"))

# 5.2 双 hub 在肿瘤 vs 正常（GSE42352 同平台）的表达重现
logmsg("Part 5.2: hub expression tumor vs normal (GSE42352) ...")
g42352 <- getGEO(filename = file.path(PROJ, "data/raw/GSE42352/GSE42352_series_matrix.txt.gz"), getGPL = FALSE)
exA <- exprs(g42352); pdA <- pData(g42352)
char_cols <- grep("^characteristics_ch1", colnames(pdA), value = TRUE)
char1 <- pdA[[char_cols[1]]]
typeA <- ifelse(grepl("biopsy", char1, ignore.case = TRUE), "Tumor",
          ifelse(grepl("MSC|osteoblast", char1, ignore.case = TRUE), "Normal", "CellLine"))
keepA <- typeA %in% c("Tumor", "Normal")
exA <- exA[, keepA]; typeA <- typeA[keepA]
exA <- exA[rownames(exA) %in% plat$ID, ]; pA <- plat[match(rownames(exA), plat$ID), ]
oA <- order(pA$Symbol, -rowMeans(exA)); exA2 <- exA[oA, ]; pA2 <- pA[oA, ]
fA <- !duplicated(pA2$Symbol); exA2 <- exA2[fA, ]; rownames(exA2) <- pA2$Symbol[fA]
tn <- data.frame()
for (g in HUBS) {
  if (!g %in% rownames(exA2)) next
  d <- data.frame(expr = as.numeric(exA2[g, ]), type = typeA)
  w <- wilcox.test(expr ~ type, d)
  tn <- rbind(tn, data.frame(gene = g, mean_tumor = mean(d$expr[d$type == "Tumor"]),
                             mean_normal = mean(d$expr[d$type == "Normal"]),
                             log2FC = mean(d$expr[d$type == "Tumor"]) - mean(d$expr[d$type == "Normal"]),
                             wilcox_p = w$p.value))
}
print(tn, row.names = FALSE)
fwrite(tn, file.path(TM, "P3_hub_tumor_vs_normal_GSE42352.csv"))

# ---------- Part 6: 决策判定 ----------
logmsg("=== P3 DECISION ===")
for (i in seq_len(nrow(km_res))) {
  cat(sprintf("%s %s: logrank p=%.4f HR/SD=%.3f (p=%.4f)\n",
      km_res$gene[i], km_res$endpoint[i], km_res$logrank_p[i], km_res$HR_per_SD[i], km_res$HR_p[i]))
}
cat("\nGSE21257 replication:\n"); print(rep_res, row.names = FALSE)
cat("\nTumor vs normal:\n"); print(tn, row.names = FALSE)
logmsg("P3 complete.")
