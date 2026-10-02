# P8_pool_enrichment_test.R — 补存 DRG 与候选池重叠的超几何检验结果（可追溯性）
# 背景基因集 = scTenifoldKnk 网络中保留的 1,501 个基因
suppressMessages(library(data.table))
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables")

pool <- fread(file.path(TM, "P2_candidate_pool_v3_matrixaxis.csv"))
POOL_N <- length(unique(pool[[1]]))
BG <- 1501L                      # scTenifoldKnk 输入/保留基因数
obs <- fread(file.path(TM, "P5_DRG_vs_candidate_pool.csv"))

res <- data.table(
  set             = obs$set,
  observed_in_pool= obs$in_pool,
  n_DGR           = obs$n,
  pool_size       = POOL_N,
  background      = BG,
  expected        = round(obs$n * POOL_N / BG, 2)
)
res[, p_value := mapply(function(k, n) {
  if (is.na(n) || n == 0) return(NA_real_)
  phyper(k - 1, POOL_N, BG - POOL_N, n, lower.tail = FALSE)
}, observed_in_pool, n_DGR)]

fwrite(res, file.path(TM, "P8_pool_enrichment_test.csv"))
print(res)
cat("\n[write]", file.path(TM, "P8_pool_enrichment_test.csv"), "\n")
