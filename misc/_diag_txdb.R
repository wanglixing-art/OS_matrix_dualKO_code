# 诊断：手动重现 TxDb 的 lazy load 失败原因
pkgdir <- "D:/rtools45/tmp/pkgsrc/TxDb.Hsapiens.UCSC.hg38.knownGene"
cat("=== DESCRIPTION Depends/Imports ===\n")
d <- read.dcf(file.path(pkgdir, "DESCRIPTION"))
print(d[, intersect(c("Depends","Imports","LinkingTo"), colnames(d)), drop = FALSE])

cat("\n=== 依赖包可用性 ===\n")
deps <- c("GenomicFeatures","AnnotationDbi","S4Vectors","IRanges","GenomicRanges",
          "BiocGenerics","GenomeInfoDb","RSQLite","DBI","methods","utils","stats")
for (p in deps) {
  ok <- requireNamespace(p, quietly = TRUE)
  v <- if (ok) as.character(packageVersion(p)) else "-"
  cat(sprintf("  %-16s %-6s %s\n", p, ifelse(ok,"OK","MISS"), v))
}

cat("\n=== 手动加载 R 代码（捕获真实错误）===\n")
r_files <- list.files(file.path(pkgdir, "R"), full.names = TRUE)
cat("R files:", length(r_files), "\n")
if (length(r_files)) {
  for (f in r_files) {
    r <- try(sys.source(f, envir = new.env()), silent = TRUE)
    if (inherits(r, "try-error")) {
      cat("!! FAILED in:", basename(f), "\n")
      cat(as.character(r), "\n")
    }
  }
}
cat("=== done ===\n")
