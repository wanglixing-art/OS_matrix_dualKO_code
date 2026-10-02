# ============================================================
# P4g-score2b: 补跑热图 + 缺失的汇总打印（前序已成功，仅修 melt API）
# ============================================================
suppressMessages({ library(data.table); library(ggplot2) })
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
INFCNV <- file.path(ROOT, "results/infercnv")
LOG <- file.path(ROOT, "logs/P4g_score2b.log")
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)

res <- fread(file.path(TM, "P4g_percell_cnv_score.csv"))
logmsg("res rows:", nrow(res))

# ---------- 打印已算出但前次未显示的汇总 ----------
logmsg("=== crosstab (inferCNV vs proxy, obs cells) ===")
ct_x <- res[group == "obs", .N, by = .(inferCNV = malignant_95, proxy = malignant_proxy)]
print(dcast(ct_x, inferCNV ~ proxy, value.var = "N", fill = 0))

logmsg("=== HMM state by celltype ===")
hf <- file.path(TM, "P4g_hmm_state_by_celltype.csv")
if (file.exists(hf) && file.info(hf)$size > 100) print(fread(hf)) else logmsg("(empty - will recount below)")

logmsg("=== concordance ===")
cf <- file.path(TM, "P4g_evidence_concordance.csv")
if (file.exists(cf) && file.info(cf)$size > 100) print(fread(cf))

# ---------- 热图（确定性解析 + as.table melt）----------
f <- file.path(INFCNV, "expr.infercnv.preliminary.dat")
ann_txt <- fread(file.path(TM, "P4g_infercnv_annotation.txt"), header = FALSE,
                 col.names = c("cell", "group"))
logmsg("reading matrix for heatmap ...")
hdr_line <- readLines(f, n = 1)
cells <- strsplit(hdr_line, "\t", fixed = TRUE)[[1]]
dt <- fread(f, sep = "\t", header = FALSE, skip = 1, showProgress = FALSE)
hm <- as.matrix(dt[, -1, drop = FALSE])
colnames(hm) <- cells
rownames(hm) <- as.character(dt[[1]])
rm(dt); gc()
stopifnot(all(colnames(hm) == ann_txt$cell))
logmsg("matrix:", nrow(hm), "x", ncol(hm))

sc <- setNames(res$cnv_score, res$cell)
ct <- setNames(res$celltype, res$cell)
obs_j <- which(ann_txt$group == "obs")
set.seed(1)
sel <- unlist(tapply(obs_j, ct[ann_txt$cell[obs_j]], function(jj)
  jj[order(-sc[ann_txt$cell[jj]])][seq_len(min(length(jj), 500))]), use.names = FALSE)
sub <- hm[, sel, drop = FALSE]
cell_ct <- ct[ann_txt$cell[sel]]
o <- order(match(cell_ct, c("Osteoblastic", "Proliferating")))
sub <- sub[, o, drop = FALSE]; cell_ct <- cell_ct[o]
logmsg("heatmap:", nrow(sub), "genes x", ncol(sub), "cells")

df <- as.data.table(as.table(sub))
setnames(df, c("gene", "cell", "val"))
df[, gene := as.character(gene)][, cell := as.character(cell)]
df[, gidx := match(gene, rownames(sub))]
# 染色体顺序：gene_order 已排序，rownames(sub) 保持该顺序
vc <- quantile(abs(df$val - 1), 0.99, na.rm = TRUE)
logmsg("val range:", round(min(df$val), 3), "-", round(max(df$val), 3),
       "| clip:", round(vc, 4))

phm <- ggplot(df, aes(cell, gidx, fill = val)) +
  geom_raster() +
  scale_fill_gradient2(low = "#00008B", mid = "white", high = "#8B0000",
                       midpoint = 1, limits = c(1 - vc, 1 + vc),
                       oob = scales::squish, name = "CNV\nresidual") +
  labs(title = "inferCNV heatmap: osteoblastic lineage (obs) cells",
       subtitle = "Blue = copy loss, red = copy gain; genes in chromosome order; cells grouped by type, sorted by CNV score",
       x = "Obs cells (Osteoblastic | Proliferating)", y = "Genes (chromosome order)") +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(FG, "P4g_infercnv_heatmap_obs.png"), phm,
       width = 14, height = 6, dpi = 300)
logmsg("saved P4g_infercnv_heatmap_obs.png")
logmsg("P4g-score2b complete.")
