# ============================================================
# P7c: TARGET-OS TMB 计算 × 14 基因签名
#
# 数据：GDC masked MAF（167 文件 / 142 case；WXS；open tier）
# 口径：
#   non-silent 突变 = Variant_Classification ∈ 9 类非同义（排除 Silent/IGR 等）
#   一 case 多 aliquot → 突变按 (Chr, Start, Ref, Alt) 去重后合并
#   TMB = non-silent counts / 38 Mb（WXS 靶向区域惯例）
#   与 RiskScore（n=88 全样本）Spearman + 高低风险（median split）Wilcoxon
# 输出：P7_tmb_scores.csv / P7_tmb_stats.csv / P7_tmb_scatter.png
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
MAFDIR <- file.path(PROJ, "data/raw/TARGET-OS_maf")
LOG <- file.path(PROJ, "logs/P7c_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}

NONSYN <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
            "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
            "Translation_Start_Site", "Nonstop_Mutation")

# case 映射
cmap <- fread(file.path(PROJ, "data/raw/TARGET-OS_maf_casemap.tsv"))
files <- list.files(MAFDIR, pattern = "\\.maf\\.gz$", full.names = TRUE)
logmsg("MAF files found:", length(files))
stopifnot(length(files) > 100)

# 逐文件解析
all_mut <- list()
for (i in seq_along(files)) {
  fid <- tools::file_path_sans_ext(basename(files[i]))
  fid <- sub("\\.maf$", "", fid)
  cs <- cmap$case_submitter_id[cmap$file_id == fid]
  if (length(cs) == 0) next
  m <- tryCatch(fread(files[i], sep = "\t", header = TRUE, skip = "Hugo_Symbol",
                      select = c("Chromosome", "Start_Position", "Reference_Allele",
                                 "Tumor_Seq_Allele2", "Variant_Classification"),
                      showProgress = FALSE), error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0) next
  m <- m[Variant_Classification %in% NONSYN]
  m$case <- cs[1]
  all_mut[[length(all_mut) + 1]] <- m
}
mut <- rbindlist(all_mut)
logmsg("total non-silent rows (pre-dedup):", nrow(mut))

# case 内去重 + 计数
mut_u <- unique(mut, by = c("case", "Chromosome", "Start_Position",
                            "Reference_Allele", "Tumor_Seq_Allele2"))
tmb <- mut_u[, .(n_nonsyn = .N), by = case]
tmb[, TMB := n_nonsyn / 38]     # muts/Mb
logmsg("cases with TMB:", nrow(tmb), "| median TMB:", round(median(tmb$TMB), 2),
       "muts/Mb | range:", round(min(tmb$TMB), 2), "-", round(max(tmb$TMB), 2))

# RiskScore（88 样本全口径，TMB 与生存无关）
sigdf <- fread(file.path(TM, "P6_signature_genes.csv"))
obj <- readRDS(file.path(PROJ, "data/processed/P1_objects.rds"))
E <- obj$E
sig <- sigdf$coef; names(sig) <- sigdf$gene
gvars <- names(sig)[names(sig) %in% rownames(E)]
X <- t(E[gvars, ]); X <- X[, apply(X, 2, sd) > 0]
rs <- as.numeric(scale(X) %*% sig[colnames(X)])
names(rs) <- rownames(X)
# submitter_id 映射（表达列名 TARGET.40.XXXX 点分隔 → case TARGET-40-XXX）
d <- data.frame(sample = names(rs), RiskScore = rs)
d$case <- gsub("\\.", "-", d$sample)
d <- merge(d, tmb, by = "case")
logmsg("merged TMB x RiskScore: n=", nrow(d))

sp <- suppressWarnings(cor.test(d$RiskScore, d$TMB, method = "spearman"))
d$grp <- factor(ifelse(d$RiskScore > median(d$RiskScore), "High", "Low"),
                levels = c("Low", "High"))
mw <- wilcox.test(TMB ~ grp, data = d)
logmsg("TMB~RiskScore spearman rho=", round(unname(sp$estimate), 3),
       " p=", signif(sp$p.value, 3),
       " | High/Low median TMB=", round(median(d$TMB[d$grp == "High"]), 2), "/",
       round(median(d$TMB[d$grp == "Low"]), 2),
       " wilcox p=", signif(mw$p.value, 3))

fwrite(d, file.path(TM, "P7_tmb_scores.csv"))
fwrite(data.frame(metric = c("spearman_rho", "spearman_p", "wilcox_p",
                             "median_TMB_high", "median_TMB_low", "n"),
                  value = c(unname(sp$estimate), sp$p.value, mw$p.value,
                            median(d$TMB[d$grp == "High"]),
                            median(d$TMB[d$grp == "Low"]), nrow(d))),
       file.path(TM, "P7_tmb_stats.csv"))

p <- ggplot(d, aes(RiskScore, TMB, color = grp)) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", fill = "grey75",
              linewidth = 0.6, alpha = 0.3) +
  scale_color_manual(values = c("Low" = "#2166AC", "High" = "#D7301F"), name = "Risk group") +
  labs(title = paste0("Tumor mutational burden vs risk score (TARGET-OS, n=", nrow(d), ")"),
       subtitle = paste0("Spearman rho=", round(unname(sp$estimate), 2),
                         ", p=", signif(sp$p.value, 2),
                         " | non-silent muts / 38 Mb (GDC masked MAF)"),
       x = "LASSO-Cox risk score", y = "TMB (mutations / Mb)") +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(FG, "P7_tmb_scatter.png"), p, width = 6.6, height = 5.2, dpi = 300)
logmsg("=== P7c complete ===")
