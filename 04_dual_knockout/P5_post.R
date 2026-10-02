# ============================================================
# P5-post: 双 KO 结果后处理（GO/KEGG 富集 + P4 证据链整合 + 可视化）
#
# 输入（由 P5_final_drg.R 产出）:
#   results/tables/P5_DRG_RUNX2_stroma.csv     (1501 基因, 含 Z/FC/p/DRG 列)
#   results/tables/P5_DRG_BUB1_prolif.csv
#   results/tables/P5_DRG_overlap_summary.csv / shared / RUNX2_specific / BUB1_specific
#   results/tables/P5_DRG_wide_both_KO.csv
#
# DRG 口径：Z > 2（正向，scTenifoldKnk 文献标准），已排除 distance 下溢伪影。
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(clusterProfiler); library(org.Hs.eg.db)
})
ROOT <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(ROOT, "results/tables"); FG <- file.path(ROOT, "results/figures")
DP <- file.path(ROOT, "data/processed")
LOG <- file.path(ROOT, "logs/P5_post.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}
Z_THR <- 2

f1 <- file.path(TM, "P5_DRG_RUNX2_stroma.csv")
f2 <- file.path(TM, "P5_DRG_BUB1_prolif.csv")
stopifnot(file.exists(f1), file.exists(f2))
dr1 <- fread(f1); dr2 <- fread(f2)
logmsg("=== P5-post start ===")
logmsg("RUNX2-KO table:", nrow(dr1), "genes | BUB1-KO table:", nrow(dr2), "genes")

# DRG 判定（与 P5_final_drg.R 同口径，兜底重算以防列缺失）
for (dt in c("dr1", "dr2")) {
  d <- get(dt)
  if (!"underflow" %in% colnames(d)) d[, underflow := distance < 1e-10]
  if (!"DRG" %in% colnames(d)) d[, DRG := (Z > Z_THR) & (!underflow)]
  assign(dt, d)
}
logmsg("RUNX2-KO DRG (Z>2):", sum(dr1$DRG), "| BUB1-KO DRG (Z>2):", sum(dr2$DRG))
logmsg("cols:", paste(colnames(dr1), collapse = ", "))

gw <- file.path(TM, "P5_DRG_wide_both_KO.csv")
if (file.exists(gw)) {
  w <- fread(gw)
} else {
  w <- merge(dr1[, .(gene, Z_RUNX2 = Z, DRG_RUNX2 = DRG)],
             dr2[, .(gene, Z_BUB1 = Z, DRG_BUB1 = DRG)], by = "gene", all = TRUE)
}
s1 <- dr1[DRG == TRUE]$gene; s2 <- dr2[DRG == TRUE]$gene
ov <- intersect(s1, s2); only1 <- setdiff(s1, s2); only2 <- setdiff(s2, s1)
logmsg("shared:", length(ov), "| RUNX2-specific:", length(only1), "| BUB1-specific:", length(only2))

# ============================================================
# 1) GO BP / KEGG 富集（三组：shared / RUNX2-specific / BUB1-specific）
# ============================================================
run_enrich <- function(genes, tag, min_size = 5) {
  if (length(genes) < min_size) { logmsg("  [skip]", tag, "n =", length(genes), "(< 5)"); return(NULL) }
  eg <- tryCatch(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL)
  if (is.null(eg) || nrow(eg) == 0) { logmsg("  [skip]", tag, "symbol->entrez failed"); return(NULL) }
  logmsg("  [", tag, "] mapped", nrow(eg), "/", length(genes), "to ENTREZ")

  go <- tryCatch(
    enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 0.25,
             readable = TRUE),
    error = function(e) { logmsg("  GO error:", conditionMessage(e)); NULL })
  if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
    go <- setReadable(go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    fwrite(as.data.table(as.data.frame(go)), file.path(TM, paste0("P5_GO_", tag, ".csv")))
    logmsg("  GO BP:", nrow(as.data.frame(go)), "terms -> P5_GO_", tag, ".csv")
  } else logmsg("  GO BP: 0 enriched terms for", tag)

  # KEGG 需联网（rest.kegg.jp），代理环境下可能挂起 → 设超时
  kg <- tryCatch({
    old_to <- getOption("timeout"); options(timeout = 60)
    on.exit(options(timeout = old_to), add = TRUE)
    enrichKEGG(gene = eg$ENTREZID, organism = "hsa", pAdjustMethod = "BH",
               pvalueCutoff = 0.1, qvalueCutoff = 0.25)
  }, error = function(e) { logmsg("  KEGG error/skip:", conditionMessage(e)); NULL })
  if (!is.null(kg) && nrow(as.data.frame(kg)) > 0) {
    kg <- setReadable(kg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    fwrite(as.data.table(as.data.frame(kg)), file.path(TM, paste0("P5_KEGG_", tag, ".csv")))
    logmsg("  KEGG:", nrow(as.data.frame(kg)), "pathways -> P5_KEGG_", tag, ".csv")
  } else logmsg("  KEGG: 0 enriched pathways for", tag)
  invisible(go)
}

logmsg("--- enrichment: shared ---");          run_enrich(ov,    "shared")
logmsg("--- enrichment: RUNX2_specific ---");  run_enrich(only1, "RUNX2_specific")
logmsg("--- enrichment: BUB1_specific ---");   run_enrich(only2, "BUB1_specific")
# RUNX2 全部 DRG（不限 specific，样本量大些）
logmsg("--- enrichment: RUNX2_all ---");       run_enrich(s1,    "RUNX2_all")
logmsg("--- enrichment: BUB1_all ---");        run_enrich(s2,    "BUB1_all")

# ============================================================
# 2) 与 P4 证据链整合
# ============================================================
# 2a) hub 基因在两 KO 中的扰动强度
hub <- c("RUNX2", "BUB1")
hb <- w[gene %in% hub]
if (nrow(hb) > 0) {
  fwrite(hb, file.path(TM, "P5_hub_perturbation.csv"))
  logmsg("hub perturbation:"); print(hb)
}

# 2b) DRG 与 P4 恶性细胞 hub 表达基因 / 候选池交叉
cand_f <- file.path(TM, "P2_candidate_pool_v3_matrixaxis.csv")
if (file.exists(cand_f)) {
  cp <- fread(cand_f)
  gcol <- grep("gene|symbol", colnames(cp), ignore.case = TRUE, value = TRUE)[1]
  pool <- unique(cp[[gcol]])
  logmsg("P2 candidate pool size:", length(pool))
  ev <- data.table(
    set = c("RUNX2-KO DRG (Z>2)", "BUB1-KO DRG (Z>2)", "Shared DRG", "RUNX2-specific", "BUB1-specific"),
    n   = c(length(s1), length(s2), length(ov), length(only1), length(only2)),
    in_pool = c(sum(s1 %in% pool), sum(s2 %in% pool), sum(ov %in% pool),
                sum(only1 %in% pool), sum(only2 %in% pool)))
  ev[, pct := round(100 * in_pool / n, 1)]
  fwrite(ev, file.path(TM, "P5_DRG_vs_candidate_pool.csv"))
  logmsg("DRG ∩ P2 candidate pool:"); print(ev)
}

# ============================================================
# 3) 可视化
# ============================================================
FG2 <- FG; dir.create(FG2, showWarnings = FALSE, recursive = TRUE)

# 3a) 双 KO DRG 韦恩（手绘两圆，避免 VennDiagram 依赖）
vt <- file.path(TM, "P5_DRG_overlap_summary.csv")
if (file.exists(vt)) {
  cmp <- fread(vt)
  pv <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white") +
    labs(title = "DRG overlap between the two knockouts",
         subtitle = paste0("RUNX2-KO n=", length(s1), " | BUB1-KO n=", length(s2))) +
    theme_void(base_size = 12) +
    annotate("text", x = 0.38, y = 0.5, label = length(only1), size = 7, fontface = "bold", color = "#E31A1C") +
    annotate("text", x = 0.5,  y = 0.5, label = length(ov),    size = 7, fontface = "bold", color = "#333333") +
    annotate("text", x = 0.62, y = 0.5, label = length(only2), size = 7, fontface = "bold", color = "#1F78B4") +
    annotate("text", x = 0.30, y = 0.30, label = "RUNX2-KO\nspecific", size = 3.6, color = "#E31A1C") +
    annotate("text", x = 0.70, y = 0.30, label = "BUB1-KO\nspecific",  size = 3.6, color = "#1F78B4") +
    annotate("text", x = 0.50, y = 0.72, label = "Shared", size = 3.6, color = "#333333") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.15, 0.85), clip = "off") +
    annotate("path", x = 0.38 + 0.14 * cos(seq(0, 2*pi, length.out = 100)),
                     y = 0.5  + 0.22 * sin(seq(0, 2*pi, length.out = 100)),
             color = "#E31A1C", linewidth = 0.9) +
    annotate("path", x = 0.62 + 0.14 * cos(seq(0, 2*pi, length.out = 100)),
                     y = 0.5  + 0.22 * sin(seq(0, 2*pi, length.out = 100)),
             color = "#1F78B4", linewidth = 0.9)
  ggsave(file.path(FG2, "P5_DRG_venn.png"), pv, width = 5.5, height = 4.5, dpi = 300)
  logmsg("saved P5_DRG_venn.png")
}

# 3b) 双 KO Z 值散点（全部 1501 基因，标出 hub 与 DRG）
if (nrow(w) > 0) {
  w[, cat := "ns"]
  w[gene %in% ov,    cat := "Shared DRG"]
  w[gene %in% only1, cat := "RUNX2-specific"]
  w[gene %in% only2, cat := "BUB1-specific"]
  w[, is_hub := gene %in% hub]
  ord <- c("ns", "Shared DRG", "RUNX2-specific", "BUB1-specific")
  w[, cat := factor(cat, levels = ord)]
  pal <- c("ns" = "grey78", "Shared DRG" = "#6A3D9A",
           "RUNX2-specific" = "#E31A1C", "BUB1-specific" = "#1F78B4")
  pz <- ggplot(w, aes(Z_RUNX2, Z_BUB1)) +
    geom_hline(yintercept = Z_THR, linetype = 2, color = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = Z_THR, linetype = 2, color = "grey45", linewidth = 0.4) +
    geom_point(data = w[cat == "ns"], color = "grey80", size = 1.1, alpha = 0.6) +
    geom_point(data = w[cat != "ns"], aes(color = cat), size = 2.2, alpha = 0.9) +
    ggrepel::geom_text_repel(data = w[cat != "ns" | is_hub == TRUE],
                             aes(label = gene, color = cat),
                             size = 2.6, max.overlaps = 30, show.legend = FALSE,
                             segment.color = "grey60", segment.size = 0.25) +
    scale_color_manual(values = pal, name = "") +
    labs(title = "Perturbation landscape of the two knockouts",
         subtitle = "Each point = one of the 1,501 HVGs; dashed lines = Z = 2",
         x = "Z in RUNX2-KO (stromal compartment)",
         y = "Z in BUB1-KO (proliferating compartment)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  ggsave(file.path(FG2, "P5_dualKO_Z_landscape.png"), pz, width = 7.2, height = 6.4, dpi = 300)
  logmsg("saved P5_dualKO_Z_landscape.png")
}

# 3c) 轴标志基因扰动对比（ECM/基质 vs 细胞周期）
ecm <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL11A1","FN1","SPP1",
         "MMP2","MMP9","MMP14","LUM","ASPN","POSTN","THBS2","SPARC","BGN","DCN",
         "PCOLCE","COL6A3","TIMP1","SERPINE1")
cellcycle <- c("BUB1","BUB1B","CCNB1","CDK1","TOP2A","UBE2C","CDC20","AURKA","AURKB",
               "CENPF","MKI67","PLK1","CCNA2","NDC80","TTK","RRM2","NUSAP1","ASPM")
sel <- c(ecm, cellcycle)
cmp_g <- merge(dr1[gene %in% sel, .(gene, Z_RUNX2 = Z, p_RUNX2 = p.value)],
               dr2[gene %in% sel, .(gene, Z_BUB1 = Z, p_BUB1 = p.value)],
               by = "gene", all = TRUE)
if (nrow(cmp_g) > 0) {
  cmp_g[, group := ifelse(gene %in% ecm, "ECM / stromal program", "Cell-cycle program")]
  cmp_g <- cmp_g[order(group, -pmax(Z_RUNX2, Z_BUB1, na.rm = TRUE))]
  fwrite(cmp_g, file.path(TM, "P5_ECM_cellcycle_perturbation.csv"))
  logmsg("ECM/cell-cycle perturbation table:"); print(cmp_g)

  lg <- melt(cmp_g, id.vars = c("gene", "group"),
             measure.vars = c("Z_RUNX2", "Z_BUB1"),
             variable.name = "KO", value.name = "Z")
  lg[, KO := ifelse(KO == "Z_RUNX2", "RUNX2-KO", "BUB1-KO")]
  lg[, gene := factor(gene, levels = rev(cmp_g$gene))]
  pg <- ggplot(lg, aes(gene, Z, fill = KO)) +
    geom_col(position = "dodge", width = 0.72, color = "grey20", linewidth = 0.25) +
    facet_wrap(~group, scales = "free_y") + coord_flip() +
    scale_fill_manual(values = c("RUNX2-KO" = "#E31A1C", "BUB1-KO" = "#1F78B4"), name = "") +
    labs(title = "Axis-marker genes: perturbation strength in each knockout",
         subtitle = "Bars = Z-score from scTenifoldKnk differential regulation",
         x = NULL, y = "Z-score", fill = "") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  ggsave(file.path(FG2, "P5_axis_markers_perturbation.png"), pg,
         width = 10, height = 7, dpi = 300)
  logmsg("saved P5_axis_markers_perturbation.png")
}

# 3d) 共享 DRG 的 Z 一致性
if (length(ov) > 3) {
  so <- merge(dr1[gene %in% ov, .(gene, Z1 = Z)], dr2[gene %in% ov, .(gene, Z2 = Z)], by = "gene")
  rho <- suppressWarnings(cor(so$Z1, so$Z2, method = "spearman"))
  saveRDS(list(rho = rho, n = nrow(so), data = so), file.path(DP, "P5_shared_concordance.rds"))
  logmsg("shared DRG Z correlation (Spearman):", round(rho, 3), "| n =", nrow(so))
  ps <- ggplot(so, aes(Z1, Z2)) +
    geom_point(color = "#6A3D9A", size = 2.4, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey25", linewidth = 0.7) +
    ggrepel::geom_text_repel(aes(label = gene), size = 2.6, max.overlaps = 20,
                             segment.color = "grey60", segment.size = 0.25) +
    labs(title = "Shared DRGs: perturbation concordance",
         subtitle = paste0("Spearman rho = ", round(rho, 3), " (n = ", nrow(so), ")"),
         x = "Z in RUNX2-KO", y = "Z in BUB1-KO") +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(FG2, "P5_shared_DRG_concordance.png"), ps,
         width = 6.4, height = 5.8, dpi = 300)
  logmsg("saved P5_shared_DRG_concordance.png")
} else logmsg("shared DRG too few (n =", length(ov), ") for concordance plot")

# 3e) GO 富集条形图
for (nm in c("shared", "RUNX2_specific", "BUB1_specific", "RUNX2_all", "BUB1_all")) {
  ef <- file.path(TM, paste0("P5_GO_", nm, ".csv"))
  if (!file.exists(ef)) next
  eg <- fread(ef)
  if (nrow(eg) == 0) next
  eg <- head(eg[order(p.adjust)], 15)
  eg[, Description := factor(Description, levels = rev(Description))]
  pe <- ggplot(eg, aes(Description, -log10(p.adjust), fill = Count)) +
    geom_col(width = 0.72, color = "grey20", linewidth = 0.25) +
    coord_flip() +
    scale_fill_gradient(low = "#FDD49E", high = "#B30000", name = "Count") +
    labs(title = paste0("GO BP enrichment: ", gsub("_", " ", nm), " DRGs"),
         subtitle = paste0("top ", nrow(eg), " terms by FDR"),
         x = NULL, y = "-log10(FDR)") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(FG2, paste0("P5_GO_", nm, ".png")), pe,
         width = 9.5, height = 5.5, dpi = 300)
  logmsg("saved P5_GO_", nm, ".png")
}

# 3f) 证据链汇总表（P4 hub 定位 × P4g CNV × P5 DRG）
ev_files <- c(
  hub_expr   = file.path(TM, "P4_hub_expr_by_celltype.csv"),
  cnv_malig  = file.path(TM, "P4g_hub_expr_vs_infercnv_malignant.csv"),
  cellchat   = file.path(TM, "P4e_hub_outgoing_by_sample.csv"),
  drg        = file.path(TM, "P5_hub_perturbation.csv"))
for (i in seq_along(ev_files)) {
  if (file.exists(ev_files[i])) logmsg("evidence chain [",
      names(ev_files)[i], "] OK:", basename(ev_files[i]))
  else logmsg("evidence chain [", names(ev_files)[i], "] MISSING")
}

logmsg("=== P5-post complete ===")
