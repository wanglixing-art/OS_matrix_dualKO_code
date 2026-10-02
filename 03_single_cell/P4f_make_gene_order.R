# ============================================================
# 生成 inferCNV 所需的 gene order 文件（gene / chr / start / end）
# 方案：直接用 loadDb 读 TxDb sqlite（绕过 TxDb 包 lazy-load 安装问题）
# 输出：data/raw/infercnv_ref/gene_order_hg38.txt（按染色体+位置排序）
# ============================================================
suppressMessages({ library(AnnotationDbi); library(GenomicFeatures); library(data.table) })
PROJ <- "D:/projects/OS_matrix_dualKO"
OUT <- file.path(PROJ, "data/raw/infercnv_ref/gene_order_hg38.txt")
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

sq <- "D:/rtools45/tmp/pkgsrc/TxDb.Hsapiens.UCSC.hg38.knownGene/inst/extdata/TxDb.Hsapiens.UCSC.hg38.knownGene.sqlite"
logmsg("Loading TxDb from sqlite (", round(file.size(sq)/1e6,1), "MB ) ...")
txdb <- loadDb(sq)

logmsg("Extracting gene coordinates ...")
g <- GenomicFeatures::genes(txdb)      # GRanges, 每个 gene_id 一个区间
logmsg("total gene ranges:", length(g))

dt <- data.table(
  gene  = names(g),
  chr   = as.character(GenomeInfoDb::seqnames(g)),
  start = as.integer(GenomicRanges::start(g)),
  end   = as.integer(GenomicRanges::end(g))
)

# 只保留标准染色体，按 cytoBand 排序
keep_chr <- paste0("chr", c(1:22, "X", "Y"))
dt <- dt[chr %in% keep_chr]
dt[, chr_num := match(chr, keep_chr)]
setorder(dt, chr_num, start, end)
dt[, chr_num := NULL]

logmsg("genes on standard chromosomes:", nrow(dt))
logmsg("per-chromosome counts:")
print(dt[, .N, by = chr][order(match(chr, keep_chr))])

# gene_id 是 Entrez，转成 gene symbol（inferCNV 要与表达矩阵的 symbol 匹配）
logmsg("Mapping Entrez -> gene symbol ...")
suppressMessages(library(org.Hs.eg.db))
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
        keys = unique(dt$gene), columns = "SYMBOL", keytype = "ENTREZID"))
map <- as.data.table(map)
map <- map[!is.na(SYMBOL)]
map <- map[!duplicated(ENTREZID)]
dt <- merge(dt, map, by.x = "gene", by.y = "ENTREZID", all.x = FALSE)
setnames(dt, "SYMBOL", "symbol")
logmsg("with symbol:", nrow(dt))

# inferCNV 格式：symbol chr start end（tab 分隔，无表头）
out_dt <- dt[, .(symbol, chr, start, end)]
out_dt <- out_dt[!duplicated(symbol)]
out_dt[, chr_num := match(chr, keep_chr)]
setorder(out_dt, chr_num, start)
out_dt[, chr_num := NULL]
fwrite(out_dt, OUT, sep = "\t", col.names = FALSE, quote = FALSE)
logmsg("written:", OUT)
logmsg("lines:", nrow(out_dt))
logmsg("head:")
print(head(out_dt, 5))
