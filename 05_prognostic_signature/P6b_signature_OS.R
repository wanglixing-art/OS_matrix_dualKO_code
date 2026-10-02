# ============================================================
# P6b: 敏感性分析 — OS 端点对齐训练版签名
# 目的: 分离外部验证弱的两种解释
#   H1 端点错配（EFS 训练 vs OS 验证）→ OS-训练版外部应显著改善
#   H2 人群/平台差异 → OS-训练版外部同样弱
# 流程与 P6 完全一致，仅端点换 OS
# ============================================================
suppressMessages({
  library(data.table); library(survival); library(glmnet); library(timeROC)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables")
LOG <- file.path(PROJ, "logs/P6b_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
set.seed(2026)

obj <- readRDS(file.path(PROJ, "data/processed/P1_objects.rds"))
E <- obj$E
clin <- fread(file.path(PROJ, "data/processed/TARGET-OS_clinical.csv"), data.table = FALSE)
pool <- fread(file.path(TM, "P2_candidate_pool_v3_matrixaxis.csv"), data.table = FALSE)$gene

clin$id_dot <- gsub("-", ".", clin$submitter_id)
common <- intersect(colnames(E), clin$id_dot)
clin <- clin[match(common, clin$id_dot), ]; E <- E[, common]
ok <- !is.na(clin$OS_time_days) & !is.na(clin$OS_event) & clin$OS_time_days > 0
clin <- clin[ok, ]; E <- E[, ok]
logmsg("OS-evaluable samples:", ncol(E), "| events:", sum(clin$OS_event))

y_t <- clin$OS_time_days; y_e <- clin$OS_event
genes <- pool[pool %in% rownames(E)]
uni <- t(vapply(genes, function(g) {
  d <- data.frame(t = y_t, e = y_e, x = as.numeric(E[g, ]))
  s <- summary(coxph(Surv(t, e) ~ x, d))
  c(HR = s$coef[2], p = s$coef[5])
}, numeric(2)))
uni <- data.frame(gene = rownames(uni), uni, row.names = NULL)
keep <- uni$p < 0.1 & !is.na(uni$p)
logmsg("univariate OS Cox p<0.1:", sum(keep), "/", nrow(uni))
sel <- uni$gene[keep]

X <- t(E[sel, ]); X <- X[, apply(X, 2, sd) > 0]; Xz <- scale(X)
cv.fit <- cv.glmnet(Xz, Surv(y_t, y_e), family = "cox", alpha = 1,
                    nfolds = 10, type.measure = "deviance")
n_1se <- sum(coef(cv.fit, s = "lambda.1se") != 0)
n_min <- sum(coef(cv.fit, s = "lambda.min") != 0)
lam <- if (n_1se >= 5) "lambda.1se" else "lambda.min"
logmsg("LASSO: 1se keeps", n_1se, "| min keeps", n_min, "-> using", lam)
beta <- as.numeric(coef(cv.fit, s = lam)); names(beta) <- rownames(coef(cv.fit, s = lam))
sig <- sort(beta[beta != 0], decreasing = TRUE)[order(-abs(beta[beta != 0]))]
sig <- beta[beta != 0]; sig <- sig[order(-abs(sig))]
logmsg("OS-signature:", length(sig), " genes: ", paste(names(sig), collapse = ","))
fwrite(data.frame(gene = names(sig), coef = as.numeric(sig)),
       file.path(TM, "P6b_signature_OS_genes.csv"))

rs <- as.numeric(Xz[, names(sig)] %*% sig)
s <- summary(coxph(Surv(y_t, y_e) ~ rs))
logmsg("TRAIN-OS: C=", round(s$concordance["C"], 3),
       " p=", signif(s$coef[5], 3))

# ---------- GSE21257 外部（端点一致: OS）----------
suppressMessages(library(GEOquery))
gse <- getGEO(filename = file.path(PROJ, "data/raw/GSE21257/series_matrix.txt.gz"),
              getGPL = FALSE)
ex <- exprs(gse); pd <- pData(gse)
plat <- fread(file.path(PROJ, "data/raw/GSE42352/GPL10295_probe_symbol.tsv"),
              data.table = FALSE, select = c("ID", "Symbol"))
plat <- plat[plat$Symbol != "" & !is.na(plat$Symbol), ]; plat <- plat[!duplicated(plat$ID), ]
ex <- ex[rownames(ex) %in% plat$ID, ]; plat <- plat[match(rownames(ex), plat$ID), ]
mmean <- rowMeans(ex); ord <- order(plat$Symbol, -mmean)
ex2 <- ex[ord, ]; pl2 <- plat[ord, ]; first <- !duplicated(pl2$Symbol)
ex2 <- ex2[first, ]; rownames(ex2) <- pl2$Symbol[first]

chars <- lapply(grep("^characteristics_ch1", colnames(pd), value = TRUE), function(cn) pd[[cn]])
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
os_time <- os_time * 30.44
os_event <- ifelse(grepl("^Deceased", st, ignore.case = TRUE), 1,
             ifelse(grepl("^Alive", st, ignore.case = TRUE), 0, NA))

miss <- setdiff(names(sig), rownames(ex2))
logmsg("missing in GSE21257:", length(miss),
       if (length(miss) > 0) paste0(" (", paste(miss, collapse = ","), ")") else "")
gvars <- intersect(names(sig), rownames(ex2))
Xg <- t(ex2[gvars, , drop = FALSE]); Xg <- Xg[, apply(Xg, 2, sd) > 0]
Xgz <- scale(Xg)
rs_g <- as.numeric(Xgz[, intersect(names(sig), colnames(Xg))] %*%
                     sig[intersect(names(sig), colnames(Xg))])
dd <- data.frame(rs = rs_g, t = os_time, e = os_event)
dd <- dd[!is.na(dd$t) & !is.na(dd$e) & dd$t > 0, ]
s2 <- summary(coxph(Surv(t, e) ~ rs, dd))
sd0 <- survdiff(Surv(t, e) ~ ifelse(rs > median(rs), "High", "Low"), dd)
p_lr <- 1 - pchisq(sd0$chisq, length(sd0$n) - 1)
auc2 <- data.frame()
for (yr in c(2, 3, 5)) {
  tt <- yr * 365.25
  if (sum(dd$t >= tt, na.rm = TRUE) < 8) next
  r <- timeROC(T = dd$t, delta = dd$e, marker = dd$rs, cause = 1, times = tt, ROC = TRUE)
  auc2 <- rbind(auc2, data.frame(year = yr, AUC = r$AUC[2]))
}
logmsg("EXT-OS (OS-trained): C=", round(s2$concordance["C"], 3),
       " HR=", round(s2$coef[2], 2), " [", round(s2$conf.int[3], 2), "-",
       round(s2$conf.int[4], 2), "] p=", signif(s2$coef[5], 3),
       " | logrank p=", signif(p_lr, 3))
if (nrow(auc2) > 0) logmsg("EXT time AUC: ", paste0(auc2$year, "y=",
       round(auc2$AUC, 3), collapse = " "))
fwrite(data.frame(metric = c("C_EXT_OS", "HR", "HR_lo", "HR_hi", "p_cox",
                             "p_logrank", paste0("AUC_", auc2$year, "y")),
                  value = c(s2$concordance["C"], s2$coef[2], s2$conf.int[3],
                            s2$conf.int[4], s2$coef[5], p_lr, auc2$AUC)),
       file.path(TM, "P6b_external_performance.csv"))
logmsg("=== P6b complete ===")
