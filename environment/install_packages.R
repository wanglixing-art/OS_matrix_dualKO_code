# Install the R dependency set for this project.
#
# Tested with R 4.6.1 on Windows 11 (ucrt). For exact versions used in the
# published analysis, see R_sessionInfo.txt in this folder.
#
# Windows users: a C++ toolchain is required for scTenifoldKnk / CellChat /
# infercnv. Install Rtools 45 (matches R 4.5.x and R 4.6.x) and make sure
# pkgbuild::has_build_tools(debug = TRUE) returns TRUE before running this.

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(install.packages.check.source = "no")   # prefer binaries on Windows

# ---- CRAN ----
cran <- c(
  "data.table", "ggplot2", "patchwork", "survival", "survminer", "timeROC",
  "glmnet", "Matrix", "jsonlite", "igraph", "tidygraph", "ggraph", "ggpubr",
  "dplyr", "NMF", "superpc", "randomForestSRC", "plsRglm", "mboost",
  "GEOquery", "WGCNA", "estimate", "BiocManager", "remotes"
)
install.packages(setdiff(cran, rownames(installed.packages())))

# ---- Bioconductor ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c(
  "edgeR", "limma", "AnnotationDbi", "org.Hs.eg.db", "GenomicFeatures",
  "clusterProfiler", "infercnv", "SingleR", "illuminaHumanv4.db"
)
BiocManager::install(setdiff(bioc, rownames(installed.packages())), ask = FALSE)

# ---- GitHub ----
remotes::install_github("jinworks/CellChat")            # CellChat 2.2.0.9001
remotes::install_github("cailab-tamu/scTenifoldKnk")    # 1.1
remotes::install_github("cailab-tamu/scTenifoldNet")    # 1.4
# harmony is on CRAN in recent versions; install from GitHub if unavailable:
if (!requireNamespace("harmony", quietly = TRUE))
  remotes::install_github("immunogenomics/harmony")

# ---- verify ----
needed <- c(cran, bioc, "harmony", "CellChat", "scTenifoldKnk", "scTenifoldNet")
missing <- needed[!vapply(needed, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
if (length(missing)) {
  message("STILL MISSING: ", paste(missing, collapse = ", "))
} else {
  message("All dependencies installed.")
}
