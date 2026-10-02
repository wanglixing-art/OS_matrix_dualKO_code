# ============================================================
# P0_build_probe_map.R
# 从 GEO 平台注释生成 探针 -> gene symbol 映射表（P1 与 P6c 的依赖文件）
#
# 产物:
#   data/raw/GSE42352/GPL10295_probe_symbol.tsv   完整平台注释表(制表符分隔, 含 ID/Symbol 列)
#   data/raw/GEO_expand/GPL14951_probe_symbol.tsv 两列映射表(ID,Symbol; 逗号分隔)
#
# 用法: Rscript P0_build_probe_map.R [项目根目录]
# 前置: 已运行 P0_download_geo_expand.sh（GPL10295 family.soft.gz）
# 说明: GPL14951 的映射优先走 illuminaHumanv4.db（无需下载 2.4 GB 平台文件）；
#       仅当该注释包不可用时才回退到 GEOquery 在线获取平台文件。
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "D:/projects/OS_matrix_dualKO"

GSE42352 <- file.path(ROOT, "data/raw/GSE42352")
EXPAND   <- file.path(ROOT, "data/raw/GEO_expand")
dir.create(EXPAND, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ---------- 1) GPL10295 (Illumina HumanWG-6 v3.0) ----------
# P1_analysis.R 用 fread(select = c("ID","Symbol")) 读取，故必须是含表头的制表符文本，
# 且保留完整平台表（48,701 行）——不要裁剪成两列。
out1 <- file.path(GSE42352, "GPL10295_probe_symbol.tsv")
soft1 <- file.path(GSE42352, "GPL10295_family.soft.gz")

if (!file.exists(out1)) {
  if (!file.exists(soft1)) {
    stop("缺少 ", soft1, "\n请先运行 P0_download_geo_expand.sh")
  }
  suppressMessages({ library(GEOquery); library(data.table) })
  msg("读取 GPL10295 平台注释 (约 320 MB soft 文件, 需要几分钟) ...")
  gpl <- getGEO(filename = soft1)
  tb  <- Table(gpl)
  msg("平台表:", nrow(tb), "行 x", ncol(tb), "列")
  if (!all(c("ID", "Symbol") %in% colnames(tb))) {
    stop("GPL10295 平台表中缺少 ID 或 Symbol 列, 实际列名: ",
         paste(colnames(tb), collapse = ", "))
  }
  fwrite(tb, out1, sep = "\t")
  msg("写出", out1, "|", nrow(tb), "行")
} else {
  msg("已存在, 跳过:", basename(out1))
}

# ---------- 2) GPL14951 (Illumina HumanHT-12 WG-DASL V4.0) ----------
# P6c_external_expand.R 优先用 illuminaHumanv4.db 的 PROBEID -> SYMBOL 映射；
# 该注释包已覆盖 GSE39055 矩阵中 14/14 个签名基因，故默认走这条路，不必下平台文件。
out2 <- file.path(EXPAND, "GPL14951_probe_symbol.tsv")

if (!file.exists(out2)) {
  mp <- NULL
  if (requireNamespace("illuminaHumanv4.db", quietly = TRUE)) {
    suppressMessages({ library(illuminaHumanv4.db); library(AnnotationDbi) })
    ks  <- keys(illuminaHumanv4.db, keytype = "PROBEID")
    ann <- AnnotationDbi::select(illuminaHumanv4.db, keys = ks,
                                 columns = "SYMBOL", keytype = "PROBEID")
    mp  <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", c("PROBEID", "SYMBOL")]
    names(mp) <- c("ID", "Symbol")
    mp <- mp[!duplicated(mp$ID), ]
    msg("illuminaHumanv4.db 映射条目:", nrow(mp))
  }
  if (is.null(mp) || nrow(mp) < 1000) {
    # 回退：在线获取平台文件（约 2.4 GB，慢）
    msg("回退到 GEOquery 在线获取 GPL14951 ...")
    suppressMessages(library(GEOquery))
    tb2 <- Table(getGEO("GPL14951", AnnotGPL = FALSE))
    sym_col <- grep("symbol|ilmn_gene|gene_name", colnames(tb2), ignore.case = TRUE, value = TRUE)[1]
    mp <- data.frame(ID = as.character(tb2[["ID"]]),
                     Symbol = as.character(tb2[[sym_col]]), stringsAsFactors = FALSE)
    mp <- mp[!is.na(mp$Symbol) & mp$Symbol != "" & mp$Symbol != "---", ]
    mp <- mp[!duplicated(mp$ID), ]
  }
  write.csv(mp, out2, row.names = FALSE, quote = FALSE)
  msg("写出", out2, "|", nrow(mp), "行")
} else {
  msg("已存在, 跳过:", basename(out2))
}

msg("全部完成。")
