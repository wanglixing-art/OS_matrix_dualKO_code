cat("=== 手动测试 loadDb ===\n")
suppressMessages({ library(AnnotationDbi); library(GenomicFeatures) })
sq <- "D:/rtools45/tmp/pkgsrc/TxDb.Hsapiens.UCSC.hg38.knownGene/inst/extdata/TxDb.Hsapiens.UCSC.hg38.knownGene.sqlite"
cat("sqlite exists:", file.exists(sq), "| size MB:", round(file.size(sq)/1e6, 1), "\n")
db <- try(loadDb(sq), silent = TRUE)
if (inherits(db, "try-error")) {
  cat("!! loadDb FAILED:\n"); cat(as.character(db), "\n")
} else {
  cat("loadDb OK | class:", class(db)[1], "\n")
  cat("seqlevels:", length(GenomeInfoDb::seqlevels(db)), "\n")
}
