suppressMessages({ library(Seurat); library(scTenifoldNet); library(data.table) })
seu <- readRDS("data/processed/P4_seurat_annotated.rds")
ct <- as.character(seu$celltype_marker)
idx <- which(ct == "Proliferating"); set.seed(1); idx <- sample(idx, 400)
sub <- subset(seu, cells = colnames(seu)[idx])
cnt <- as.matrix(LayerData(sub, layer = "counts"))
cnt <- cnt[rowSums(cnt > 0) >= 3, , drop = FALSE]
cpm <- t(t(cnt) / colSums(cnt)) * 1e6
v <- apply(log2(cpm + 1), 1, var)
g <- rownames(cnt)[order(v, decreasing = TRUE)]
for (N in c(1000, 2000, 3000)) {
  gg <- unique(c(g[seq_len(N)], "BUB1"))
  m <- as.matrix(cnt[gg, , drop = FALSE])
  set.seed(1); smp <- sample(ncol(m), 200)
  t0 <- Sys.time()
  net <- pcNet(m[, smp], nComp = 3, nCores = 4, q = 0.9, verbose = FALSE)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("N=%-5d net1 of 200 cells: %6.1f s | nnz=%d\n", length(gg), el, length(net@x)))
}
