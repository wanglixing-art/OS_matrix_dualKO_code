# ============================================================
# P6c: 外部验证扩展 —— GSE39055 / GSE33382 / GSE87624
#
# 背景（用户指令）：GSE21257 外部 OS 端点弱（C=0.547）但转移维度显著
#   (Wilcoxon p=0.0021) → 扩展 2-3 个独立数据集再验证
# 选集依据（series matrix 逐个核实）：
#   1) GSE39055 (n=37, GPL14951 Illumina HT-12 WG-DASL V4.0)
#      - FFPE 活检, recurrence Y/N + time-to-recurrence-or-FU(months)
#      - 端点 = EFS（与签名训练端点完全一致，最强外部验证）
#   2) GSE33382 (n=49 活检, GPL10295 = 与 GSE21257 完全同平台)
#      - LUMC 高级别 OS 化疗前活检, "metastasis within 5yrs: yes/no"
#      - 转移终点验证 + 同平台直接迁移
#   3) GSE87624 (n≈45 肿瘤, RNA-seq HiSeq2000, counts)
#      - tumor type: primary vs metastasis（转移灶组织 vs 原发）
#      - 转移维度 RNA-seq 独立平台验证
# 口径（与 P6 严格一致，不重拟合）：
#   - 数据集内基因 z-score（scale）→ 训练系数线性打分
#   - 缺失基因：从系数向量剔除（贡献=0），并在日志记录缺失数
#   - KM median split / timeROC / Wilcoxon 同 P6
# ============================================================
suppressMessages({
  library(data.table); library(survival); library(timeROC)
  library(GEOquery); library(ggplot2); library(survminer); library(patchwork)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM  <- file.path(PROJ, "results/tables")
FG  <- file.path(PROJ, "results/figures")
EXP <- file.path(PROJ, "data/raw/GEO_expand")
LOG <- file.path(PROJ, "logs/P6c_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
set.seed(2026)

# ---------- Part 0: 训练签名 ----------
sigdf <- fread(file.path(TM, "P6_signature_genes.csv"), data.table = FALSE)
sig <- sigdf$coef; names(sig) <- sigdf$gene
logmsg("P6 signature loaded:", length(sig), "genes")

perf_row <- function(time, event, score) {
  s <- summary(coxph(Surv(time, event) ~ score))
  # concordance 向量按位置取：[1]=C, [2]=se(C)（跨 survival 版本名字不同）
  c(C_index = s$concordance[1], se = s$concordance[2],
    HR = s$coef[2], HR_lo = s$conf.int[3], HR_hi = s$conf.int[4], p = s$coef[5])
}
km_eval <- function(time, event, score) {
  grp <- factor(ifelse(score > median(score), "High", "Low"), levels = c("Low", "High"))
  sd <- survdiff(Surv(time, event) ~ grp)
  p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  list(grp = grp, p_logrank = p)
}
score_and_check <- function(ex2, tag) {
  # ex2: gene(symbol) x sample 表达矩阵 → z-score 打分（P6 口径）
  miss <- setdiff(names(sig), rownames(ex2))
  gvars <- names(sig)[names(sig) %in% rownames(ex2)]
  logmsg(tag, ": signature genes mapped ", length(gvars), "/", length(sig),
         if (length(miss) > 0) paste0(" | missing: ", paste(miss, collapse = ",")) else "")
  Xg <- t(ex2[gvars, , drop = FALSE])
  Xg <- Xg[, apply(Xg, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  gvars2 <- colnames(Xg)
  Xz <- scale(Xg)
  rs <- as.numeric(Xz[, gvars2, drop = FALSE] %*% sig[gvars2])
  list(rs = rs, n_mapped = length(gvars), missing = miss)
}

getfield <- function(pd, key, require_colon = TRUE) {
  chars <- lapply(grep("^characteristics_ch1", colnames(pd), value = TRUE),
                  function(cn) pd[[cn]])
  out <- rep(NA_character_, nrow(pd))
  pat <- paste0("^", key, if (require_colon) ":" else "")
  for (v in chars) {
    hit <- grepl(pat, v, ignore.case = TRUE)
    out[hit] <- sub("^.*:\\s*", "", v[hit])
  }
  out
}

summary_rows <- data.frame()

# ============================================================
# Part 1: GSE39055 —— EFS 端点外部验证（与训练端点一致）
# ============================================================
logmsg("=== Part 1: GSE39055 (n=37, EFS endpoint) ===")
g1 <- getGEO(filename = file.path(EXP, "GSE39055_series_matrix.txt.gz"), getGPL = FALSE)
ex1 <- exprs(g1); pd1 <- pData(g1)
logmsg("matrix:", nrow(ex1), "probes x", ncol(ex1), "samples")

# 平台注释 GPL14951（Illumina HumanHT-12 WG-DASL V4.0）
# 已验证：illuminaHumanv4.db PROBEID 即 ILMN_ 前缀格式；GSE39055 matrix
# 29377 探针（检测子集），14/14 签名基因全覆盖
map_file <- file.path(EXP, "GPL14951_probe_symbol.tsv")
if (!file.exists(map_file)) {
  mp <- NULL
  if (requireNamespace("illuminaHumanv4.db", quietly = TRUE)) {
    library(illuminaHumanv4.db); library(AnnotationDbi)
    ks <- keys(illuminaHumanv4.db, keytype = "PROBEID")
    ann <- AnnotationDbi::select(illuminaHumanv4.db, keys = ks,
                                 columns = "SYMBOL", keytype = "PROBEID")
    mp <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", c("PROBEID", "SYMBOL")]
    names(mp) <- c("ID", "Symbol")
    mp <- mp[!duplicated(mp$ID), ]
    logmsg("illuminaHumanv4.db mapped entries:", nrow(mp))
  }
  if (is.null(mp) || nrow(mp) < 1000) {
    soft <- getGEO("GPL14951", AnnotGPL = FALSE)
    tb <- Table(soft)
    sym_col <- grep("symbol|ilmn_gene|gene_name", colnames(tb), ignore.case = TRUE, value = TRUE)[1]
    id_col  <- grep("^id$", colnames(tb), ignore.case = TRUE, value = TRUE)[1]
    mp <- data.frame(ID = as.character(tb[[id_col]]), Symbol = as.character(tb[[sym_col]]),
                     stringsAsFactors = FALSE)
    mp <- mp[mp$Symbol != "" & !is.na(mp$Symbol) & mp$Symbol != "---", ]
    mp <- mp[!duplicated(mp$ID), ]
  }
  fwrite(mp, map_file)
} else {
  mp <- fread(map_file, data.table = FALSE)
}
mp <- mp[mp$Symbol != "" & !is.na(mp$Symbol), ]
mp <- mp[!duplicated(mp$ID), ]
ex1 <- ex1[rownames(ex1) %in% mp$ID, ]
mp <- mp[match(rownames(ex1), mp$ID), ]
mmean <- rowMeans(ex1); ord <- order(mp$Symbol, -mmean)
ex1 <- ex1[ord, ]; mp1 <- mp[ord, ]; first <- !duplicated(mp1$Symbol)
ex1m <- ex1[first, ]; rownames(ex1m) <- mp1$Symbol[first]
logmsg("GSE39055 mapped:", nrow(ex1m), "genes x", ncol(ex1m), "samples")

# 临床解析：recurrence / death / time (months)
# 注意：时间字段全名 "time until first recurrence or latest follow-up (months)"
# → require_colon = FALSE 前缀匹配（已实测踩坑）
rec  <- getfield(pd1, "recurrence")
tmon <- getfield(pd1, "time until first recurrence", require_colon = FALSE)
ef_t <- as.numeric(tmon) * 30.44
ef_e <- ifelse(rec == "Y", 1, ifelse(rec == "N", 0, NA))
nec  <- getfield(pd1, "percent necrosis")
nec_hi <- as.numeric(grepl(">90|9[0-9]%|90-95|>95|>99", nec))   # Huvos III/IV 近似（>=90% 坏死）
sample1 <- colnames(ex1m)
logmsg("EFS events:", sum(ef_e == 1, na.rm = TRUE), "/", sum(!is.na(ef_e)),
       "| necrosis>=90%:", sum(nec_hi == 1, na.rm = TRUE))

sc1 <- score_and_check(ex1m, "GSE39055")
d1 <- data.frame(sample = sample1, RiskScore = sc1$rs,
                 EFS_time = ef_t, EFS_event = ef_e, nec_hi = nec_hi)
d1 <- d1[!is.na(d1$EFS_time) & !is.na(d1$EFS_event) & d1$EFS_time > 0, ]
logmsg("GSE39055 evaluable: n=", nrow(d1), " events=", sum(d1$EFS_event))
p1 <- perf_row(d1$EFS_time, d1$EFS_event, d1$RiskScore)
k1 <- km_eval(d1$EFS_time, d1$EFS_event, d1$RiskScore)
logmsg("GSE39055 EFS: C=", round(p1[1], 3), " HR=", round(p1[3], 2),
       " [", round(p1[4], 2), "-", round(p1[5], 2), "] p=", signif(p1[6], 3),
       " | KM logrank p=", signif(k1$p_logrank, 3))
auc1 <- data.frame()
for (yr in c(1, 2, 3, 4, 5)) {
  tt <- yr * 365.25
  if (sum(d1$EFS_time >= tt, na.rm = TRUE) < 8) next
  r <- timeROC(T = d1$EFS_time, delta = d1$EFS_event, marker = d1$RiskScore,
               cause = 1, times = tt, ROC = TRUE)
  auc1 <- rbind(auc1, data.frame(year = yr, AUC = r$AUC[2]))
}
logmsg("GSE39055 time AUC:"); print(auc1)
# 化疗反应分层（Huvos>=III 近似）下签名独立性
if (sum(d1$nec_hi == 1, na.rm = TRUE) >= 5 && sum(d1$nec_hi == 0, na.rm = TRUE) >= 5) {
  mv1 <- coxph(Surv(EFS_time, EFS_event) ~ RiskScore + nec_hi, d1)
  logmsg("GSE39055 multivariable (+necrosis>=90%): RiskScore p=",
         signif(summary(mv1)$coef["RiskScore", 5], 3))
  fwrite(data.frame(term = names(coef(mv1)),
                    HR = summary(mv1)$conf.int[, 1],
                    lo = summary(mv1)$conf.int[, 3], hi = summary(mv1)$conf.int[, 4],
                    p = summary(mv1)$coef[, 5]),
         file.path(TM, "P6c_GSE39055_multivariable.csv"))
}
fwrite(d1, file.path(TM, "P6c_riskscore_GSE39055.csv"))
if (nrow(auc1) > 0) fwrite(auc1, file.path(TM, "P6c_timeAUC_GSE39055.csv"))
summary_rows <- rbind(summary_rows, data.frame(
  dataset = "GSE39055", n = nrow(d1), endpoint = "EFS (recurrence)",
  C_index = round(p1[1], 3), C_se = round(p1[2], 3),
  HR = round(p1[3], 2), HR_lo = round(p1[4], 2), HR_hi = round(p1[5], 2),
  cox_p = signif(p1[6], 3), logrank_p = signif(k1$p_logrank, 3),
  genes_mapped = sc1$n_mapped))
png(file.path(FG, "P6c_KM_GSE39055.png"), 1500, 1300, res = 150)
d1$grp <- k1$grp
fit1 <- survfit(Surv(EFS_time, EFS_event) ~ grp, d1)
print(ggsurvplot(fit1, data = d1,
                 pval = paste0("logrank p = ", signif(k1$p_logrank, 3)),
                 risk.table = TRUE, conf.int = TRUE,
                 palette = c("#2166AC", "#D7301F"),
                 legend.labs = c("Low risk", "High risk"),
                 title = paste0("GSE39055 external validation (n=", nrow(d1), "), EFS endpoint"),
                 xlab = "Time (days)", surv.median.line = "hv"))
dev.off()

# ============================================================
# Part 2: GSE33382 —— 同平台 (GPL10295) + 5 年转移终点
# ============================================================
logmsg("=== Part 2: GSE33382 (n=49 biopsy, GPL10295) ===")
g2 <- getGEO(filename = file.path(EXP, "GSE33382_series_matrix.txt.gz"), getGPL = FALSE)
ex2 <- exprs(g2); pd2 <- pData(g2)
logmsg("matrix:", nrow(ex2), "probes x", ncol(ex2), "samples")
mp2 <- fread(file.path(PROJ, "data/raw/GSE42352/GPL10295_probe_symbol.tsv"),
             data.table = FALSE, select = c("ID", "Symbol"))
mp2 <- mp2[mp2$Symbol != "" & !is.na(mp2$Symbol), ]
mp2 <- mp2[!duplicated(mp2$ID), ]
ex2 <- ex2[rownames(ex2) %in% mp2$ID, ]
mp2 <- mp2[match(rownames(ex2), mp2$ID), ]
mmean2 <- rowMeans(ex2); ord2 <- order(mp2$Symbol, -mmean2)
ex2 <- ex2[ord2, ]; mp2o <- mp2[ord2, ]; first2 <- !duplicated(mp2o$Symbol)
ex2m <- ex2[first2, ]; rownames(ex2m) <- mp2o$Symbol[first2]
logmsg("GSE33382 mapped:", nrow(ex2m), "genes x", ncol(ex2m), "samples")

typ <- getfield(pd2, "type")
met5 <- getfield(pd2, "metastasis within 5yrs")
# 列出全部 characteristics 键，确认是否还有生存字段
allkeys <- unique(unlist(lapply(grep("^characteristics_ch1", colnames(pd2), value = TRUE),
                function(cn) sub("^([^:]*):.*", "\\1", pd2[[cn]]))))
logmsg("GSE33382 characteristic keys: ", paste(allkeys, collapse = " | "))
keep2 <- which(typ == "biopsy")
sc2 <- score_and_check(ex2m, "GSE33382")
d2 <- data.frame(sample = colnames(ex2m)[keep2], RiskScore = sc2$rs[keep2],
                 met5 = met5[keep2])
d2 <- d2[!is.na(d2$met5), ]
d2$met5 <- factor(d2$met5, levels = c("no", "yes"))
logmsg("GSE33382 biopsy with met5 info: n=", nrow(d2),
       " (yes=", sum(d2$met5 == "yes"), " no=", sum(d2$met5 == "no"), ")")
if (all(table(d2$met5) >= 5)) {
  mw5 <- wilcox.test(RiskScore ~ met5, data = d2)
  logmsg("GSE33382 met5 (yes vs no): median ",
         round(median(d2$RiskScore[d2$met5 == "yes"]), 3), " vs ",
         round(median(d2$RiskScore[d2$met5 == "no"]), 3),
         " | Wilcoxon p=", signif(mw5$p.value, 3))
  # 逻辑回归方向一致性
  lg5 <- glm(met5 == "yes" ~ RiskScore, data = d2, family = binomial)
  logmsg("GSE33382 logistic OR=", round(exp(coef(lg5)[2]), 2),
         " p=", signif(summary(lg5)$coef[2, 4], 3))
  summary_rows <- rbind(summary_rows, data.frame(
    dataset = "GSE33382", n = nrow(d2), endpoint = "Met within 5 yrs",
    C_index = NA, C_se = NA, HR = exp(coef(lg5)[2]), HR_lo = NA, HR_hi = NA,
    cox_p = signif(summary(lg5)$coef[2, 4], 3), logrank_p = signif(mw5$p.value, 3),
    genes_mapped = sc2$n_mapped))
  fwrite(d2, file.path(TM, "P6c_riskscore_GSE33382.csv"))
}

# ============================================================
# Part 3: GSE87624 —— RNA-seq 转移维度验证
# ============================================================
logmsg("=== Part 3: GSE87624 (RNA-seq, met vs primary) ===")
cnt_file <- file.path(EXP, "GSE87624_Human_masked.txt.gz")
cnt_ok <- FALSE; tmp_head <- NULL
if (file.exists(cnt_file)) {
  hd <- try(readLines(gzfile(cnt_file), n = 1), silent = TRUE)
  closeAllConnections()
  if (!inherits(hd, "try-error") && length(hd) > 0) { cnt_ok <- TRUE; tmp_head <- hd }
}
if (cnt_ok) {
  cn <- tmp_head
  logmsg("GSE87624 counts header: ", substr(cn, 1, 200))
  c3 <- fread(cnt_file, data.table = FALSE)
  logmsg("counts dim:", nrow(c3), "x", ncol(c3))
  # 首列 = 基因标识
  gid <- c3[[1]]; c3[[1]] <- NULL
  M <- as.matrix(c3); rownames(M) <- gid
  # CPM -> log2(CPM+1)
  lib <- colSums(M); M <- log2(t(t(M) / lib) * 1e6 + 1)
  gid2 <- gsub("\\..*$", "", rownames(M))   # 去版本号（若有 Ensembl id）
  # masked 文件可能直接是 symbol 或 Ensembl；用 org.Hs 映射兜底
  is_sym <- mean(grepl("^[A-Z][A-Z0-9-]+$|^C[0-9]orf", rownames(M))) > 0.5
  logmsg("gene id looks like SYMBOL:", is_sym)
  ex3m <- M; rownames(ex3m) <- gid2
  if (!is_sym) {
    # Ensembl -> symbol（org.Hs.eg.db）
    library(org.Hs.eg.db); library(AnnotationDbi)
    mp3 <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(ex3m),
                                 keytype = "ENSEMBL", columns = "SYMBOL")
    mp3 <- mp3[!is.na(mp3$SYMBOL) & !duplicated(mp3$ENSEMBL), ]
    ex3m <- ex3m[rownames(ex3m) %in% mp3$ENSEMBL, ]
    rownames(ex3m) <- mp3$SYMBOL[match(rownames(ex3m), mp3$ENSEMBL)]
  }
  # 聚合重复 symbol（取均值最高）
  gm <- rowMeans(ex3m); o3 <- order(rownames(ex3m), -gm)
  ex3m <- ex3m[o3, , drop = FALSE]
  f3 <- !duplicated(rownames(ex3m)); ex3m <- ex3m[f3, ]
  logmsg("GSE87624 final gene matrix:", nrow(ex3m), "x", ncol(ex3m))

  g3 <- getGEO(filename = file.path(EXP, "GSE87624_series_matrix.txt.gz"), getGPL = FALSE)
  pd3 <- pData(g3)
  tt3 <- getfield(pd3, "tumor type")
  names(tt3) <- pd3$geo_accession
  # sample title 与 counts 列名对齐
  st3 <- getfield(pd3, "title"); if (all(is.na(st3))) st3 <- pd3$title
  names(st3) <- pd3$geo_accession
  logmsg("GSE87624 counts colnames head:", paste(colnames(ex3m)[1:6], collapse = " | "))
  logmsg("GSE87624 sample titles head:", paste(head(unname(st3), 6), collapse = " | "))
  common3 <- intersect(colnames(ex3m), names(tt3))
  if (length(common3) < 5) {  # 列名是 title 而非 GSM
    title_map <- st3; names(title_map) <- gsub(" .*", "", title_map)
    idx <- match(colnames(ex3m), title_map)
    tt3 <- tt3[idx]
    keep3 <- !is.na(tt3) & tt3 %in% c("primary", "metastasis")
  } else {
    keep3 <- common3 %in% names(tt3) & tt3[common3] %in% c("primary", "metastasis")
  }
  sc3 <- score_and_check(ex3m, "GSE87624")
  if (length(common3) < 5) {
    d3 <- data.frame(sample = colnames(ex3m), RiskScore = sc3$rs, tt = tt3)
  } else {
    d3 <- data.frame(sample = common3, RiskScore = sc3$rs[match(common3, colnames(ex3m))],
                     tt = tt3[common3])
  }
  d3 <- d3[!is.na(d3$tt) & d3$tt %in% c("primary", "metastasis"), ]
  logmsg("GSE87624 evaluable: primary=", sum(d3$tt == "primary"),
         " metastasis=", sum(d3$tt == "metastasis"))
  if (all(table(d3$tt) >= 5)) {
    mw3 <- wilcox.test(RiskScore ~ tt, data = d3)
    logmsg("GSE87624 (metastasis vs primary): median ",
           round(median(d3$RiskScore[d3$tt == "metastasis"]), 3), " vs ",
           round(median(d3$RiskScore[d3$tt == "primary"]), 3),
           " | Wilcoxon p=", signif(mw3$p.value, 3))
    summary_rows <- rbind(summary_rows, data.frame(
      dataset = "GSE87624", n = nrow(d3), endpoint = "Met tissue vs primary",
      C_index = NA, C_se = NA, HR = NA, HR_lo = NA, HR_hi = NA,
      cox_p = NA, logrank_p = signif(mw3$p.value, 3), genes_mapped = sc3$n_mapped))
    fwrite(d3, file.path(TM, "P6c_riskscore_GSE87624.csv"))
  }
} else {
  logmsg("GSE87624 counts file missing -> skipped")
}

# ============================================================
# Part 4: 汇总 + 图
# ============================================================
fwrite(summary_rows, file.path(TM, "P6c_external_performance_expand.csv"))
logmsg("=== summary ==="); print(summary_rows)

# 图 A：三数据集 RiskScore 转移/事件对比箱线图
box_list <- list()
if (exists("d2") && nrow(d2) > 0) {
  box_list[["GSE33382\n(met within 5yrs)"]] <- d2[, c("RiskScore", "met5")]
  names(box_list[[length(box_list)]]) <- c("RiskScore", "grp")
}
if (exists("d3") && nrow(d3) > 0) {
  box_list[["GSE87624\n(met tissue)"]] <- d3[, c("RiskScore", "tt")]
  names(box_list[[length(box_list)]]) <- c("RiskScore", "grp")
}
if (length(box_list) > 0) {
  bd <- rbindlist(box_list, idcol = "ds")
  bd$grp <- factor(bd$grp, levels = c("no", "primary", "yes", "metastasis"))
  bd$grp <- factor(bd$grp, levels = levels(droplevels(bd$grp)))
  pb <- ggplot(bd, aes(grp, RiskScore, fill = grp)) +
    geom_boxplot(width = 0.5, outlier.shape = 21, outlier.size = 2,
                 color = "grey25", linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.55, color = "grey25") +
    stat_summary(fun = mean, geom = "point", size = 3, shape = 18, color = "black") +
    scale_fill_manual(values = c("no" = "#2166AC", "primary" = "#2166AC",
                                 "yes" = "#D7301F", "metastasis" = "#D7301F"),
                      guide = "none") +
    facet_wrap(~ds, scales = "free_x") +
    labs(title = "LASSO-Cox risk score by metastasis status (expanded external cohorts)",
         x = NULL, y = "Risk score (z-score based)") +
    theme_bw(base_size = 12) + theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(FG, "P6c_met_boxplots.png"), pb, width = 8, height = 4.6, dpi = 300)
}

# 图 B：森林图（C-index/HR 汇总，训练 + 全部外部）
fr <- fread(file.path(TM, "P6_train_performance.csv"))
train_c <- fr$value[grepl("^C_index", fr$metric)][1]   # 旧文件 metric 名为 "C_index.C"
ext_old <- fread(file.path(TM, "P6_external_performance.csv"))
old_c <- ext_old$value[grepl("^C_EXT_OS", ext_old$metric)][1]
n_eff <- c(TRUE, TRUE, !is.na(summary_rows$C_index))    # GSE33382/GSE87624 无 C（转移终点）
forest_df <- data.frame(
  cohort = c("TARGET-OS (training)", "GSE21257 (OS)", summary_rows$dataset),
  est = c(train_c, old_c, summary_rows$C_index),
  lo = c(NA, NA, summary_rows$C_index - 1.96 * summary_rows$C_se),
  hi = c(NA, NA, summary_rows$C_index + 1.96 * summary_rows$C_se),
  type = c("Training", "External", rep("External", nrow(summary_rows))))
forest_df <- forest_df[n_eff & !is.na(forest_df$est), ]
forest_df$cohort <- factor(forest_df$cohort, levels = rev(forest_df$cohort))
pf <- ggplot(forest_df, aes(est, cohort, color = type)) +
  geom_vline(xintercept = 0.5, linetype = 2, color = "grey55") +
  geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = 0.8,
                  position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("Training" = "#E31A1C", "External" = "#1F78B4"), name = "") +
  coord_cartesian(xlim = c(0.4, 1)) +
  labs(title = "Concordance index across cohorts (14-gene EFS signature)",
       subtitle = "GSE33382 / GSE87624 are metastasis-endpoint cohorts (no C-index; see boxplot figure)",
       x = "C-index", y = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(FG, "P6c_forest_Cindex.png"), pf, width = 8.2, height = 4.4, dpi = 300)

logmsg("=== P6c complete ===")
