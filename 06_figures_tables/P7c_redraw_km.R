# 重画 TARGET 双枢纽 4 张 KM 图（修复 em dash 乱码 -> ASCII "-"）
suppressMessages({library(data.table); library(survival); library(survminer); library(ggplot2)})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")

E <- read.csv(file.path(PROJ, "data/processed/TARGET-OS_voom_E.csv"), row.names = 1, check.names = FALSE)
clin <- fread(file.path(PROJ, "data/processed/TARGET-OS_clinical.csv"), data.table = FALSE)
rownames(clin) <- gsub("-", ".", clin$submitter_id)
common <- intersect(colnames(E), rownames(clin))
E <- E[, common]; clin <- clin[common, ]
cat("samples:", ncol(E), "\n")

km_plot <- function(expr, time, event, gene, endpoint, fout) {
  d <- data.frame(time = time, event = event, expr = as.numeric(expr))
  d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$expr) & d$time > 0, ]
  d$group <- factor(ifelse(d$expr > median(d$expr), "High", "Low"), levels = c("Low", "High"))
  n_h <- sum(d$group == "High"); n_l <- sum(d$group == "Low")
  fit <- survfit(Surv(time, event) ~ group, data = d)
  png(fout, 1400, 1200, res = 150)
  print(ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                   palette = c("#2166AC", "#D7301F"),
                   legend.labs = c(paste0("Low (n=", n_l, ")"), paste0("High (n=", n_h, ")")),
                   title = paste0(gene, " - ", endpoint, " (TARGET-OS, median split)"),
                   xlab = "Days", ylab = paste0(endpoint, " probability"),
                   risk.table.height = 0.25))
  dev.off()
}
for (g in c("BUB1", "RUNX2")) {
  km_plot(E[g, ], clin$OS_time_days, clin$OS_event, g, "OS",
          file.path(FG, paste0("P3_KM_", g, "_OS.png")))
  km_plot(E[g, ], clin$EFS_time_days, clin$EFS_event, g, "EFS",
          file.path(FG, paste0("P3_KM_", g, "_EFS.png")))
  cat("done", g, "\n")
}
