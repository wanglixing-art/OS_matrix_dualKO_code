# ============================================================
# P1: TARGET-OS 合并 -> 归一化 -> WGCNA + GSE42352 DEG -> 候选池
# 项目: OS_matrix_dualKO  日期: 2026-09-28
# ============================================================
suppressMessages({
  library(jsonlite); library(edgeR); library(limma); library(WGCNA)
  library(clusterProfiler); library(org.Hs.eg.db); library(GEOquery)
  library(data.table); library(ggplot2)
})
options(stringsAsFactors = FALSE)
PROJ <- "D:/projects/OS_matrix_dualKO"
dir.create(file.path(PROJ, "results/tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJ, "results/figures"), recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ---------- Part 1: 合并 88 个 TARGET-OS count TSV ----------
logmsg("Part 1: merging TARGET-OS counts ...")
man <- fromJSON(file.path(PROJ, "data/raw/TARGET-OS/gdc_files_manifest.json"), simplifyVector = FALSE)
hits <- man$data$hits
filemap <- data.frame(file_id  = vapply(hits, function(h) h$file_id, ""),
                      case_id = vapply(hits, function(h) h$cases[[1]]$submitter_id, ""))
filemap$expr_id <- gsub("-", ".", filemap$case_id)

fdir <- file.path(PROJ, "data/raw/TARGET-OS/files")
sym_per_file <- lapply(seq_len(nrow(filemap)), function(i) {
  f <- file.path(fdir, paste0(filemap$file_id[i], ".tsv"))
  d <- read.delim(f, comment.char = "#", check.names = FALSE)
  d <- d[d$gene_type == "protein_coding", ]
  tapply(d$unstranded, d$gene_name, sum)
})
# ---------- Part 1b: 按样本列绑定生成最终 count 矩阵 ----------
all_sym <- sort(unique(unlist(lapply(sym_per_file, names))))
M <- matrix(0L, nrow = length(all_sym), ncol = length(sym_per_file),
            dimnames = list(all_sym, filemap$expr_id))
for (i in seq_along(sym_per_file)) {
  v <- sym_per_file[[i]]
  M[names(v), i] <- as.integer(v)
}
logmsg("count matrix:", nrow(M), "genes x", ncol(M), "samples")

# ---------- Part 2: edgeR 过滤 + TMM + voom ----------
logmsg("Part 2: normalization ...")
clin <- fread(file.path(PROJ, "data/processed/TARGET-OS_clinical.csv"), data.table = FALSE)
rownames(clin) <- gsub("-", ".", clin$submitter_id)
common <- intersect(colnames(M), rownames(clin))
M2 <- M[, common]; clin2 <- clin[common, ]
dge <- DGEList(counts = M2, samples = clin2)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ 1, data = dge$samples)
v <- voom(dge, design, plot = FALSE)
E <- v$E
logmsg("filtered genes:", nrow(dge), " | voom E dim:", dim(E)[1], "x", dim(E)[2])

fwrite(data.frame(gene = rownames(M2), M2, check.names = FALSE),
       file.path(PROJ, "data/processed/TARGET-OS_counts_symbol.csv"))
fwrite(data.frame(gene = rownames(E), E, check.names = FALSE),
       file.path(PROJ, "data/processed/TARGET-OS_voom_E.csv"))

# ---------- Part 3: WGCNA ----------
logmsg("Part 3: WGCNA ...")
allowWGCNAThreads()
datExpr <- t(E)
madv <- apply(datExpr, 2, mad)
datExpr <- datExpr[, order(-madv)[1:min(8000, sum(madv > 0))]]
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
logmsg("WGCNA input:", dim(datExpr)[1], "samples x", dim(datExpr)[2], "genes")

traits <- data.frame(
  EFS_event = as.numeric(clin2[rownames(datExpr), "EFS_event"]),
  OS_event  = as.numeric(clin2[rownames(datExpr), "OS_event"]),
  age       = as.numeric(clin2[rownames(datExpr), "age_at_diagnosis_years"]),
  male      = ifelse(clin2[rownames(datExpr), "gender"] == "male", 1, 0))
rownames(traits) <- rownames(datExpr)

powers <- c(1:10, 12, 14, 16, 18, 20)
# 注: BiocGenerics(经 org.Hs.eg.db 加载) 的 S4 泛型 cor 会遮蔽 WGCNA 内部 cor,
# 故全流程统一用 WGCNA 自带 bicor 规避
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0,
                         networkType = "unsigned", corFnc = "bicor")
pw <- sft$powerEstimate; if (is.na(pw)) pw <- powers[which(sft$fitIndices$SFT.R.sq > 0.8)[1]]
if (is.na(pw)) pw <- 6
logmsg("soft power =", pw, " (SFT R2 =", signif(sft$fitIndices$SFT.R.sq[match(pw, powers)], 3), ")")

png(file.path(PROJ, "results/figures/P1_WGCNA_softthreshold.png"), 1600, 800, res = 150)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold (power)", ylab = "Scale-free Topology R^2", type = "b", pch = 19)
abline(h = 0.8, col = "red", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k., xlab = "power", ylab = "mean connectivity", type = "b", pch = 19)
dev.off()

net <- blockwiseModules(datExpr, power = pw, networkType = "unsigned", TOMType = "unsigned",
                        corType = "bicor", minModuleSize = 30, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, verbose = 0, maxBlockSize = 8000)
modColors <- labels2colors(net$colors)
logmsg("modules found:", length(unique(net$colors)))

MEs <- net$MEs
MEs <- orderMEs(MEs)
mtc <- stats::cor(MEs, traits, use = "p")
pmtc <- corPvalueStudent(mtc, nrow(MEs))
res_mtc <- data.frame(module = rownames(mtc), round(mtc, 3), signif(pmtc, 3), check.names = FALSE)
colnames(res_mtc) <- c("module", paste0("r_", colnames(mtc)), paste0("p_", colnames(mtc)))
fwrite(res_mtc, file.path(PROJ, "results/tables/P1_WGCNA_module_trait_correlations.csv"))
fwrite(data.frame(gene = colnames(datExpr), module_num = net$colors, module_color = modColors),
       file.path(PROJ, "results/tables/P1_WGCNA_module_assignments.csv"))

png(file.path(PROJ, "results/figures/P1_WGCNA_module_trait.png"), 2000, 1400, res = 150)
labeledHeatmap(Matrix = mtc, xLabels = colnames(traits), yLabels = rownames(mtc),
               ySymbols = rownames(mtc), colorLabels = FALSE,
               setStdMargins = FALSE, textMatrix = paste0(signif(mtc, 2), "\n(", signif(pmtc, 1), ")"),
               cex.text = 0.6, colors = blueWhiteRed(50), main = "Module-trait relationships")
dev.off()

# 显著模块: 与 EFS_event 或 OS_event 相关 p<0.05
sig_mods <- unique(sub("^ME", "", rownames(pmtc)[pmtc[, "EFS_event"] < 0.05 | pmtc[, "OS_event"] < 0.05]))
mod_genes <- colnames(datExpr)[net$colors %in% sig_mods]
logmsg("significant modules (EFS/OS event p<0.05):", paste(sig_mods, collapse = ","),
       " | module genes:", length(mod_genes))

# ---------- Part 4: GSE42352 DEG (tumor biopsy vs MSC/osteoblast) ----------
logmsg("Part 4: GSE42352 DEG ...")
gse <- getGEO(filename = file.path(PROJ, "data/raw/GSE42352/GSE42352_series_matrix.txt.gz"),
              getGPL = FALSE, AnnotGPL = FALSE)
ex <- exprs(gse)
pd <- pData(gse)
char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
# 细胞类型在第一列 (biopsy / cell type: MSC|osteoblast / cell line: ...)
char1 <- pd[[char_cols[1]]]
type <- ifelse(grepl("biopsy", char1, ignore.case = TRUE), "Tumor",
        ifelse(grepl("MSC|osteoblast", char1, ignore.case = TRUE), "Normal", "CellLine"))
print(table(type))
keep_s <- type %in% c("Tumor", "Normal")
ex <- ex[, keep_s]; type <- type[keep_s]
# 探针 -> symbol, 同基因取平均表达最高的探针
plat <- fread(file.path(PROJ, "data/raw/GSE42352/GPL10295_probe_symbol.tsv"), data.table = FALSE,
              select = c("ID", "Symbol"))
plat <- plat[plat$Symbol != "" & !is.na(plat$Symbol), ]
plat <- plat[!duplicated(plat$ID), ]
ex <- ex[rownames(ex) %in% plat$ID, ]
plat <- plat[match(rownames(ex), plat$ID), ]
mmean <- rowMeans(ex)
ord <- order(plat$Symbol, -mmean)
ex2 <- ex[ord, ]; pl2 <- plat[ord, ]
first <- !duplicated(pl2$Symbol)
ex2 <- ex2[first, ]; rownames(ex2) <- pl2$Symbol[first]
logmsg("GSE42352 mapped:", nrow(ex2), "genes |", sum(type == "Tumor"), "tumor vs", sum(type == "Normal"), "normal")

design2 <- model.matrix(~ 0 + factor(type)); colnames(design2) <- c("Normal", "Tumor")
fit <- lmFit(ex2, design2)
contrast <- makeContrasts(Tumor - Normal, levels = design2)
fit2 <- eBayes(contrasts.fit(fit, contrast))
deg <- topTable(fit2, number = Inf, adjust.method = "BH")
deg$gene <- rownames(deg)
fwrite(deg, file.path(PROJ, "results/tables/P1_DEG_GSE42352_tumor_vs_normal.csv"))
deg_up <- deg$gene[deg$adj.P.Val < 0.05 & deg$logFC > 1]
deg_dn <- deg$gene[deg$adj.P.Val < 0.05 & deg$logFC < -1]
logmsg("DEG up:", length(deg_up), " down:", length(deg_dn))

png(file.path(PROJ, "results/figures/P1_DEG_volcano.png"), 1400, 1200, res = 150)
deg$sig <- with(deg, ifelse(adj.P.Val < 0.05 & logFC > 1, "Up",
                     ifelse(adj.P.Val < 0.05 & logFC < -1, "Down", "NS")))
ggplot(deg, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c(Up = "#D7301F", Down = "#2166AC", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), lty = 2, col = "grey40") +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "grey40") +
  theme_bw() + labs(title = "GSE42352: OS biopsy vs MSC/Osteoblast", x = "log2FC", y = "-log10(FDR)")
dev.off()

# ---------- Part 5: 候选池 = 肿瘤上调 DEG ∩ 显著模块基因 ----------
pool <- intersect(deg_up, mod_genes)
logmsg("candidate pool size:", length(pool))
fwrite(data.frame(gene = pool), file.path(PROJ, "results/tables/P1_candidate_pool.csv"), row.names = FALSE)

# ---------- Part 6: GO/KEGG 富集 (Q1) ----------
logmsg("Part 6: enrichment ...")
eg <- bitr(pool, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
               readable = TRUE, pvalueCutoff = 0.05, qvalueCutoff = 0.1)
fwrite(as.data.frame(go), file.path(PROJ, "results/tables/P1_GO_BP_enrichment.csv"))
kegg_ok <- TRUE
tryCatch({
  kegg <- enrichKEGG(eg$ENTREZID, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.1)
  kegg <- setReadable(kegg, org.Hs.eg.db, "ENTREZID")
  fwrite(as.data.frame(kegg), file.path(PROJ, "results/tables/P1_KEGG_enrichment.csv"))
}, error = function(e) { kegg_ok <<- FALSE; logmsg("KEGG failed (offline?):", conditionMessage(e)) })
if (nrow(as.data.frame(go)) > 0) {
  png(file.path(PROJ, "results/figures/P1_GO_BP_dotplot.png"), 1800, 1400, res = 150)
  print(dotplot(go, showCategory = 20) + ggtitle("Candidate pool - GO BP") + theme_bw())
  dev.off()
}

# ---------- 汇报 ----------
logmsg("=== P1 SUMMARY ===")
cat("TARGET-OS samples:", ncol(M2), " filtered genes:", nrow(dge), "\n")
cat("WGCNA: power =", pw, " modules =", length(unique(net$colors)),
    " sig modules =", paste(sig_mods, collapse = ","), "\n")
cat("GSE42352 DEG: up =", length(deg_up), " down =", length(deg_dn), "\n")
cat("Candidate pool =", length(pool), "genes\n")
if (nrow(as.data.frame(go)) > 0) {
  cat("Top GO BP terms:\n")
  print(head(as.data.frame(go)[, c("Description", "p.adjust", "Count")], 10))
}
if (length(pool) > 0) {
  cat("Pool contains COL11A1/THBS2/COL1A1/COL10A1? ",
      paste(c("COL11A1","THBS2","COL1A1","COL10A1") %in% pool, collapse = ","), "\n")
}
saveRDS(list(counts = M2, E = E, deg = deg, net = net, MEs = MEs, mtc = mtc, pmtc = pmtc,
             pool = pool, deg_up = deg_up, mod_genes = mod_genes),
        file.path(PROJ, "data/processed/P1_objects.rds"))
logmsg("P1 complete.")
