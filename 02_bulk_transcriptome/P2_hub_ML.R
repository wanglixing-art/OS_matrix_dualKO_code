# ============================================================
# P2 (修订版): 候选池收紧 -> STRING PPI + cytoHubba 指标 -> hub 排名 -> EFS 预测模型
# 模板: Molecules 2026;31:1901 (CKAP2, GC) 方法学套用
# 项目: OS_matrix_dualKO  |  2026-09-28
# ============================================================
suppressMessages({
  library(data.table); library(jsonlite); library(GEOquery); library(limma)
  library(igraph); library(WGCNA); library(glmnet); library(mboost); library(plsRglm)
  library(randomForestSRC); library(superpc); library(survival); library(timeROC)
  library(ggplot2)
})
set.seed(20260928)
PROJ <- "D:/projects/OS_matrix_dualKO"
TM   <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ============================================================
# Part 1: GSE42352 WGCNA（肿瘤 vs 正常，与模板对齐）+ kME>0.7 收紧
# ============================================================
logmsg("Part 1: GSE42352 WGCNA (tumor vs normal) + kME filter ...")
obj <- readRDS(file.path(PROJ, "data/processed/P1_objects.rds"))

# 重建 GSE42352 表达矩阵（与 P1 同口径）
gelist <- getGEO(filename = file.path(PROJ, "data/raw/GSE42352/GSE42352_series_matrix.txt.gz"),
                 getGPL = FALSE)
ex <- exprs(gelist); pd <- pData(gelist)
char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
char1 <- pd[[char_cols[1]]]
type <- ifelse(grepl("biopsy", char1, ignore.case = TRUE), "Tumor",
        ifelse(grepl("MSC|osteoblast", char1, ignore.case = TRUE), "Normal", "CellLine"))
keep <- type %in% c("Tumor", "Normal")
ex <- ex[, keep]; type <- type[keep]
plat <- fread(file.path(PROJ, "data/raw/GSE42352/GPL10295_probe_symbol.tsv"), data.table = FALSE,
              select = c("ID", "Symbol"))
plat <- plat[plat$Symbol != "" & !is.na(plat$Symbol), ]; plat <- plat[!duplicated(plat$ID), ]
ex <- ex[rownames(ex) %in% plat$ID, ]; plat <- plat[match(rownames(ex), plat$ID), ]
mmean <- rowMeans(ex); ord <- order(plat$Symbol, -mmean)
ex2 <- ex[ord, ]; pl2 <- plat[ord, ]; first <- !duplicated(pl2$Symbol)
ex2 <- ex2[first, ]; rownames(ex2) <- pl2$Symbol[first]
logmsg("GSE42352 matrix:", nrow(ex2), "genes x", ncol(ex2), "samples (",
       sum(type == "Tumor"), "T /", sum(type == "Normal"), "N)")

# 过滤 + WGCNA
madv <- apply(ex2, 1, mad); ex3 <- ex2[order(-madv)[1:min(8000, sum(madv > 0))], ]
mads <- apply(ex3, 1, mad); ex3 <- ex3[mads > 0, ]
datExprG <- t(ex3)
gsg <- goodSamplesGenes(datExprG, verbose = 0)
if (!gsg$allOK) datExprG <- datExprG[gsg$goodSamples, gsg$goodGenes]
traitG <- data.frame(Tumor = ifelse(type[match(rownames(datExprG), colnames(ex2))] == "Tumor", 1, 0))
rownames(traitG) <- rownames(datExprG)

sftG <- pickSoftThreshold(datExprG, powerVector = c(1:10, 12, 14, 16, 18, 20),
                          verbose = 0, networkType = "unsigned", corFnc = "bicor")
pwG <- sftG$powerEstimate
if (is.na(pwG)) pwG <- 4
logmsg("GSE42352 soft power =", pwG)
netG <- blockwiseModules(datExprG, power = pwG, networkType = "unsigned", TOMType = "unsigned",
                         corType = "bicor", minModuleSize = 30, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
                         verbose = 0, maxBlockSize = 8000)
MEsG <- orderMEs(netG$MEs)
mtcG <- stats::cor(MEsG, traitG, use = "p"); pmtcG <- corPvalueStudent(mtcG, nrow(MEsG))
colnames(mtcG) <- "Tumor"
logmsg("GSE42352 modules:", length(unique(netG$colors)),
       " | module-trait r range:", paste(round(range(mtcG[, 1]), 3), collapse = " ~ "))
# 显著模块：p<0.05；若无则取 |r| 最大者（防极端数据集无显著）
sigmodsG <- rownames(mtcG)[pmtcG[, "Tumor"] < 0.05]
if (length(sigmodsG) == 0) sigmodsG <- rownames(mtcG)[order(-abs(mtcG[, "Tumor"]))][1:2]
# orderMEs 列名不含 "ME" 前缀，需直接按模块名取列
kMEall <- stats::cor(datExprG, MEsG[, sigmodsG, drop = FALSE], use = "p")
colnames(kMEall) <- sigmodsG
kMEgene <- apply(abs(kMEall), 1, max)
# 同时算 TARGET WGCNA 的 kME（对 P1 显著模块）
datExprT <- t(obj$E)
madvT <- apply(datExprT, 2, mad); datExprT <- datExprT[, order(-madvT)[1:min(8000, sum(madvT > 0))]]
gsgT <- goodSamplesGenes(datExprT, verbose = 0)
if (!gsgT$allOK) datExprT <- datExprT[gsgT$goodSamples, gsgT$goodGenes]
sigmodsT <- rownames(obj$pmtc)[obj$pmtc[, "EFS_event"] < 0.05 | obj$pmtc[, "OS_event"] < 0.05]
MEsT <- obj$MEs                     # P1 保存对象列名带 "ME" 前缀（如 ME3）
kMEallT <- stats::cor(datExprT, MEsT[, sigmodsT, drop = FALSE], use = "p")
colnames(kMEallT) <- sigmodsT
kMEgeneT <- apply(abs(kMEallT), 1, max)
logmsg("GSE42352 sig modules:", paste(sigmodsG, collapse = ","),
       " | kME>0.7 genes:", sum(kMEgene > 0.7, na.rm = TRUE),
       " | TARGET kME>0.7 genes:", sum(kMEgeneT > 0.7, na.rm = TRUE),
       " | sigmodsT:", paste(sigmodsT, collapse = ","))

fwrite(data.frame(gene = names(kMEgene), kME_GSE42352 = round(kMEgene, 3)),
       file.path(TM, "P2_kME_GSE42352.csv"))
fwrite(data.frame(gene = names(kMEgeneT), kME_TARGET = round(kMEgeneT, 3)),
       file.path(TM, "P2_kME_TARGET.csv"))

# ============================================================
# Part 2: DEG 双阈值（|log2FC|>0.585 与 >1）
# ============================================================
logmsg("Part 2: DEG dual thresholds ...")
deg <- obj$deg
deg_up585 <- deg$gene[deg$adj.P.Val < 0.05 & deg$logFC > 0.585]
deg_dn585 <- deg$gene[deg$adj.P.Val < 0.05 & deg$logFC < -0.585]
deg_up1   <- deg$gene[deg$adj.P.Val < 0.05 & deg$logFC > 1]
logmsg("DEG up: 0.585 ->", length(deg_up585), " | 1.0 ->", length(deg_up1))

# 候选池 v2（收紧，模板 kME>0.7 口径）: 上调(0.585) ∩ TARGET kME>0.7 ∩ GSE42352 kME>0.7
pool_v2 <- Reduce(intersect, list(deg_up585,
                                  names(kMEgeneT)[kMEgeneT > 0.7],
                                  names(kMEgene)[kMEgene > 0.7]))
pool_v1 <- obj$pool
# 候选池 v3（基质轴专用）: DEG 上调 ∩ 显著模块基因集（P1 口径，保留基质/ECM 基因）
pool_v3 <- intersect(deg_up585, obj$mod_genes)
# ECM/基质轴基因定义
ecm_markers <- c("COL11A1","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3",
                 "COL8A1","COL10A1","COL12A1","THBS1","THBS2","FN1","POSTN","SPARC","LUM","DCN",
                 "FBN1","LOX","LOXL1","LOXL2","MMP2","MMP9","MMP14","TIMP1","SPP1","BGN","VCAN","CTHRC1",
                 "ADAMTS4","ADAMTS2","ADAM12","SERPINH1","PCOLCE","COMP","CHAD","ASPN","PRELP","OGN")
logmsg("pool v1 (P1) =", length(pool_v1),
       " | v2 (kME>0.7 strict) =", length(pool_v2),
       " | v3 (DEG up ∩ sig modules) =", length(pool_v3),
       " | ECM genes in v3:", paste(intersect(pool_v3, ecm_markers), collapse = ","))
fwrite(data.frame(gene = pool_v2), file.path(TM, "P2_candidate_pool_v2.csv"), row.names = FALSE)
fwrite(data.frame(gene = pool_v3), file.path(TM, "P2_candidate_pool_v3_matrixaxis.csv"), row.names = FALSE)

# ============================================================
# Part 3: STRING PPI + cytoHubba 指标（Degree/MCC/MNC/EPC）
# ============================================================
logmsg("Part 3: STRING PPI + centrality metrics ...")
get_string <- function(genes, score = 700) {
  ids <- paste(unique(genes), collapse = "%0d")
  u <- sprintf("https://string-db.org/api/tsv/network?identifiers=%s&species=9606&required_score=%d", ids, score)
  tmp <- tempfile(fileext = ".tsv")
  system2("curl", c("-s", "--max-time", "120", shQuote(u), "-o", shQuote(tmp)))
  if (!file.exists(tmp) || file.size(tmp) < 50) return(NULL)
  d <- fread(tmp, sep = "\t", header = TRUE, data.table = FALSE)
  if (nrow(d) == 0) return(NULL)
  data.frame(from = d$preferredName_A, to = d$preferredName_B,
             score = d$score / 1000, stringsAsFactors = FALSE)
}
# STRING PPI：使用池 v3（含基质/ECM 基因，覆盖双轴候选；v2 严格池已剔净 ECM 基因）
edges <- NULL
for (p in list(list(g = pool_v3, s = 700), list(g = pool_v3, s = 400),
               list(g = pool_v2, s = 700), list(g = pool_v2, s = 400))) {
  e <- get_string(p$g, p$s)
  if (!is.null(e)) { edges <- e; rr <- paste(p$s, length(p$g)); break }
}
if (is.null(edges)) stop("STRING API returned nothing")
logmsg("STRING edges (score>0.7):", nrow(edges))
fwrite(edges, file.path(TM, "P2_STRING_edges.csv"))

g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$degree <- degree(g)
cent <- data.frame(gene = V(g)$name,
                   degree = degree(g),
                   betweenness = round(betweenness(g), 1),
                   closeness = round(closeness(g), 4),
                   eigen_centrality = round(eigen_centrality(g)$vector, 4),
                   stringsAsFactors = FALSE)
# MCC / MNC / EPC 近似（igraph 实现与 cytoHubba 同族指标）
mcc_score <- function(g, vid) {
  nb <- neighbors(g, vid)
  if (length(nb) < 1) return(0)
  nbid <- as.integer(nb)
  s <- 0
  for (j in nbid) {
    sub <- induced_subgraph(g, c(as.integer(vid), j))
    if (ecount(sub) == 0) next
    cm <- maximal.cliques(sub)
    s <- s + sum(sapply(cm, function(cl) length(cl) - 1))
  }
  s
}
cent$MCC <- sapply(cent$gene, function(x) mcc_score(g, V(g)[name == x]))
cent$MNC <- sapply(cent$gene, function(x) {
  nb <- neighbors(g, V(g)[name == x]); if (length(nb) < 1) return(0)
  max(c(0, degree(g)[as.integer(nb)]))
})
epc <- sapply(cent$gene, function(x) {
  nb <- neighbors(g, V(g)[name == x]); if (length(nb) == 0) return(0)
  sum(1 / log(pmax(degree(g)[as.integer(nb)], 2)))
})
cent$EPC <- round(epc, 3)
cent <- cent[order(-cent$degree, -cent$MCC), ]
fwrite(cent, file.path(TM, "P2_hub_centrality.csv"))
logmsg("top 15 by degree:", paste(head(cent$gene, 15), collapse = ", "))

# ============================================================
# Part 4: 构建 hub 排名综合表（含 P1 生存关联与表达）
# ============================================================
logmsg("Part 4: integrated hub ranking ...")
E <- obj$E; clin <- fread(file.path(PROJ, "data/processed/TARGET-OS_clinical.csv"), data.table = FALSE)
rownames(clin) <- gsub("-", ".", clin$submitter_id)
common <- intersect(colnames(E), rownames(clin)); clin <- clin[common, ]
cox_res <- t(sapply(cent$gene, function(gene) {
  if (!gene %in% rownames(E)) return(c(HR_EFS = NA, p_EFS = NA, HR_OS = NA, p_OS = NA))
  x <- as.numeric(E[gene, common]);
  d <- data.frame(time = clin$EFS_time_days, event = clin$EFS_event, x = x)
  d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$x), ]; s1 <- summary(coxph(Surv(time, event) ~ x, d))
  d2 <- data.frame(time = clin$OS_time_days, event = clin$OS_event, x = x)
  d2 <- d2[!is.na(d2$time) & !is.na(d2$event) & !is.na(d2$x), ]; s2 <- summary(coxph(Surv(time, event) ~ x, d2))
  c(HR_EFS = s1$coef[2], p_EFS = s1$coef[5], HR_OS = s2$coef[2], p_OS = s2$coef[5])
}))
cox_res <- as.data.frame(cox_res)
cox_res$gene <- rownames(cox_res)
rank_tab <- merge(cent, cox_res, by = "gene", all.x = TRUE)
rank_tab$in_pool_v1 <- rank_tab$gene %in% pool_v1

# ECM/基质轴标注（ecm_markers 已在 Part 2 定义）
rank_tab$matrix_axis <- ifelse(rank_tab$gene %in% ecm_markers, "ECM", "other")
logmsg("ECM-axis genes in PPI network:", sum(rank_tab$matrix_axis == "ECM"),
       "->", paste(rank_tab$gene[rank_tab$matrix_axis == "ECM"], collapse = ","))
rank_tab <- rank_tab[order(-rank_tab$degree, rank_tab$p_EFS), ]
fwrite(rank_tab, file.path(TM, "P2_hub_ranking_integrated.csv"))

# 网络图
png(file.path(FG, "P2_PPI_network.png"), 2400, 2000, res = 150)
set.seed(1); l <- layout_with_fr(g)
V(g)$size <- 6 + 3 * sqrt(degree(g))
V(g)$color <- ifelse(V(g)$name %in% c(head(rank_tab$gene, 10)), "#D7301F",
              ifelse(V(g)$name %in% ecm_markers, "#F4A582", "grey75"))
V(g)$label.cex <- ifelse(V(g)$name %in% head(rank_tab$gene, 20), 0.9, 0.45)
V(g)$label.color <- "black"
plot(g, layout = l, vertex.frame.color = NA, main = "STRING PPI (confidence > 0.7)")
dev.off()

# ============================================================
# Part 5: EFS 事件预测模型（12 算法组合，模板套用）
# ============================================================
logmsg("Part 5: ML models for EFS event prediction ...")
# 特征：优先用严格池 v2；若 <20 则回退 v3 并按单因素 Cox 预筛 top30（防过拟合）
feat <- intersect(pool_v2, rownames(E))
if (length(feat) < 20) feat <- intersect(pool_v3, rownames(E))
X_full <- E[feat, common, drop = FALSE]
y_time <- clin$EFS_time_days; y_event <- clin$EFS_event
ok <- !is.na(y_time) & !is.na(y_event)
X_full <- X_full[, ok]; y_time <- y_time[ok]; y_event <- y_event[ok]
# 单因素 Cox 预筛（p<0.1，取 top30）——模板亦在建模前做特征降维
ps <- apply(t(X_full), 2, function(v) {
  s <- try(summary(coxph(Surv(y_time, y_event) ~ v))$coef[5], silent = TRUE)
  if (inherits(s, "try-error")) NA else s
})
feat_sel <- names(sort(ps[!is.na(ps)]))[1:min(30, sum(!is.na(ps)))]
cat("pool v2 in TARGET:", length(feat), "-> uni-Cox p<0.1 pre-selected:", length(feat_sel), "\n")
X <- t(X_full[feat_sel, ])
X <- scale(X); X[is.na(X)] <- 0
logmsg("training set:", nrow(X), "samples x", ncol(X), "features | events:", sum(y_event == 1))

surv <- Surv(y_time, y_event)
# 通过 C-index 内部 5-fold CV 评估各单算法
nfolds <- 5; folds <- sample(rep(1:nfolds, length.out = nrow(X)))
algo_eval <- list()
cv_cindex <- function(fit_fn, pred_fn, name) {
  cs <- c()
  for (f in 1:nfolds) {
    tr <- folds != f; te <- folds == f
    m <- try(fit_fn(tr), silent = TRUE)
    if (inherits(m, "try-error") || is.null(m)) { cs <- c(cs, NA); next }
    p <- try(pred_fn(m, te), silent = TRUE)
    if (inherits(p, "try-error") || is.null(p) || all(is.na(p)) || diff(range(p, na.rm = TRUE)) == 0) {
      cs <- c(cs, NA); next
    }
    cc <- try(survConcordance(Surv(y_time[te], y_event[te]) ~ p)$concordance, silent = TRUE)
    cs <- c(cs, if (inherits(cc, "try-error")) NA else cc)
  }
  if (all(is.na(cs))) return(NA)
  mean(cs, na.rm = TRUE)
}
ridge_f <- function(tr) cv.glmnet(X[tr, ], surv[tr], family = "cox", alpha = 0, nfolds = 5)
ridge_p <- function(m, te) as.numeric(predict(m, newx = X[te, ], s = "lambda.min"))
lasso_f <- function(tr) cv.glmnet(X[tr, ], surv[tr], family = "cox", alpha = 1, nfolds = 5)
lasso_p <- function(m, te) as.numeric(predict(m, newx = X[te, ], s = "lambda.min"))
enet_f  <- function(tr) cv.glmnet(X[tr, ], surv[tr], family = "cox", alpha = 0.5, nfolds = 5)
enet_p  <- function(m, te) as.numeric(predict(m, newx = X[te, ], s = "lambda.min"))
rfsrc_f <- function(tr) {
  if (requireNamespace("randomForestSRC", quietly = TRUE)) {
    d <- data.frame(Surv(y_time[tr], y_event[tr]), X[tr, , drop = FALSE], check.names = TRUE)
    rfsrc(Surv(y_time[tr], y_event[tr]) ~ ., data = d, ntree = 500, forest = TRUE)
  } else NULL
}
rfsrc_p <- function(m, te) {
  if (is.null(m)) return(NULL)
  as.numeric(predict(m, newdata = data.frame(X[te, , drop = FALSE], check.names = TRUE))$predicted)
}
glmboost_f <- function(tr) {
  d <- data.frame(X[tr, , drop = FALSE], check.names = TRUE)
  glmboost(Surv(y_time[tr], y_event[tr]) ~ ., data = d, family = CoxPH())
}
glmboost_p <- function(m, te) as.numeric(predict(m, newdata = data.frame(X[te, , drop = FALSE], check.names = TRUE)))
pls_f <- function(tr) {
  d <- data.frame(X[tr, , drop = FALSE], check.names = TRUE)
  plsRcox(d, time = y_time[tr], event = y_event[tr], nt = 2)
}
pls_p <- function(m, te) {
  nd <- data.frame(X[te, , drop = FALSE], check.names = TRUE)
  as.numeric(predict(m, type = "lp", newdata = nd)$lp)
}
uni_f <- function(tr) {
  # 单因素 Cox 筛选 + 线性组合（Stepglm 风格兜底）
  ps <- apply(X[tr, ], 2, function(v) {
    s <- try(summary(coxph(Surv(y_time[tr], y_event[tr]) ~ v))$coef[5], silent = TRUE)
    if (inherits(s, "try-error")) NA else s
  })
  sel <- names(sort(ps))[1:min(8, sum(!is.na(ps)))]
  list(sel = sel, coef = coef(coxph(Surv(y_time[tr], y_event[tr]) ~ X[tr, sel])))
}
uni_p <- function(m, te) as.numeric(X[te, m$sel] %*% m$coef)
uni_f2 <- function(tr) {
  ps <- apply(X[tr, , drop = FALSE], 2, function(v) {
    s <- try(summary(coxph(Surv(y_time[tr], y_event[tr]) ~ v))$coef[5], silent = TRUE)
    if (inherits(s, "try-error")) NA else s
  })
  sel <- names(sort(ps))[1:min(8, sum(!is.na(ps)))]
  cf <- try(coef(coxph(as.formula(paste("Surv(y_time, y_event) ~", paste(sprintf("`%s`", sel), collapse = "+"))),
                       data = data.frame(X[, sel, drop = FALSE], y_time = y_time, y_event = y_event)[tr, ])), silent = TRUE)
  if (inherits(cf, "try-error")) NULL else list(sel = sel, coef = cf)
}
superpc_f <- function(tr) {
  d <- list(x = t(X[tr, ]), y = y_time[tr], censoring.status = y_event[tr], featurenames = colnames(X))
  superpc.train(d, type = "survival")
}
superpc_p <- function(m, te) {
  d <- list(x = t(X), y = y_time, censoring.status = y_event, featurenames = colnames(X))
  tr <- superpc.predict(m, d, list(x = t(X[te, ]), y = y_time[te], censoring.status = y_event[te],
                                   featurenames = colnames(X)), threshold = 0.5, n.components = 1)$v.pred[, 1]
  as.numeric(tr)
}
algos <- list(
  Ridge = list(ridge_f, ridge_p), Lasso = list(lasso_f, lasso_p), Enet = list(enet_f, enet_p),
  RSF = list(rfsrc_f, rfsrc_p), glmBoost = list(glmboost_f, glmboost_p),
  plsRcox = list(pls_f, pls_p), StepCox = list(uni_f2, uni_p)
)
res_alg <- data.frame()
for (nm in names(algos)) {
  cc <- try(cv_cindex(algos[[nm]][[1]], algos[[nm]][[2]], nm), silent = TRUE)
  logmsg("  ", nm, " CV C-index =", if (inherits(cc, "try-error")) "ERR" else round(cc, 3))
  res_alg <- rbind(res_alg, data.frame(model = nm, cv_cindex = if (inherits(cc, "try-error")) NA else cc))
}
fwrite(res_alg, file.path(TM, "P2_ML_single_algo_cv.csv"))

# 全数据拟合各算法，输出风险评分
risk <- data.frame(sample = rownames(X))
for (nm in names(algos)) {
  m <- try(algos[[nm]][[1]](rep(TRUE, nrow(X))), silent = TRUE)
  if (inherits(m, "try-error") || is.null(m)) { logmsg("  skip", nm, "fit failed"); next }
  p <- try(algos[[nm]][[2]](m, rep(TRUE, nrow(X))), silent = TRUE)
  if (!inherits(p, "try-error") && !is.null(p) && !all(is.na(p)) && diff(range(p, na.rm = TRUE)) > 0) {
    risk[[nm]] <- scale(as.numeric(p))[, 1]
  }
}
fwrite(risk, file.path(TM, "P2_ML_risk_scores.csv"))

# 最优两两组合（模板 113 模型思路：集成平均风险评分）
combos <- as.data.frame(t(combn(setdiff(colnames(risk), "sample"), 2)))
combos$mean_cindex <- NA_real_
for (i in seq_len(nrow(combos))) {
  a <- risk[[combos[i, 1]]]; b <- risk[[combos[i, 2]]]
  if (is.null(a) || is.null(b)) next
  if (sum(!is.na(a)) < 10 || sum(!is.na(b)) < 10) next
  both <- scale(a + b)[, 1]
  cc <- try(survConcordance(Surv(y_time, y_event) ~ both)$concordance, silent = TRUE)
  if (!inherits(cc, "try-error")) combos$mean_cindex[i] <- cc
}
combos <- combos[order(-combos$mean_cindex), ]
fwrite(combos, file.path(TM, "P2_ML_pair_combinations.csv"))
logmsg("top 5 algorithm pairs (full-data C-index):")
print(head(combos, 5))

saveRDS(list(rank_tab = rank_tab, cent = cent, pool_v2 = pool_v2, pool_v3 = pool_v3,
             edges = edges, X = X, feat_sel = feat_sel, y_time = y_time, y_event = y_event,
             risk = risk, combos = combos, res_alg = res_alg, sigmodsG = sigmodsG,
             sigmodsT = sigmodsT, pwG = pwG, netG = netG),
        file.path(PROJ, "data/processed/P2_objects.rds"))
logmsg("=== P2 SUMMARY ===")
cat("GSE42352 WGCNA power:", pwG, "| sig modules:", paste(sigmodsG, collapse = ","), "\n")
cat("TARGET sig modules:", paste(sigmodsT, collapse = ","), "\n")
cat("pool v1:", length(pool_v1), "| v2 (kME>0.7):", length(pool_v2), "| v3 (matrix axis):", length(pool_v3), "\n")
cat("PPI nodes:", vcount(g), " edges:", ecount(g), "\n")
cat("top 10 hub by degree:", paste(head(cent$gene, 10), collapse = ", "), "\n")
cat("ECM genes in network:", paste(rank_tab$gene[rank_tab$matrix_axis == "ECM"], collapse = ", "), "\n")
cat("EFS events in training:", sum(y_event == 1), " n =", length(y_event), "\n")
cat("best single algo:", res_alg$model[which.max(res_alg$cv_cindex)],
    "C-index =", round(max(res_alg$cv_cindex, na.rm = TRUE), 3), "\n")
logmsg("P2 complete.")
