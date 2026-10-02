# ============================================================
# P7a2: 逐张独立的 ROC / AUC 曲线图
#   用户指令："请分别画出 auc 曲线图"
#   理解：把每个队列 x 每个时间点的 ROC 曲线拆成独立单图
#        输出目录 results/figures/roc_individual/
#   打分：全部沿用 P6/P6c 口径（数据集内 z-score + 训练系数），不重拟合
#   特色：单图内嵌 AUC + 95% CI（timeROC iid=TRUE 渐近方差）
# ============================================================
suppressMessages({
  library(data.table); library(survival); library(timeROC); library(ggplot2)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables")
FG <- file.path(PROJ, "results/figures")
OUTD <- file.path(FG, "roc_individual")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(PROJ, "logs/P7a2_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}

COL <- c("TARGET-OS" = "#D7301F", "GSE21257" = "#1F78B4", "GSE39055" = "#33A02C")

# ---------- 单张 ROC 曲线 ----------
plot_single_roc <- function(time, event, marker, horizon, cohort, endpoint,
                            col, fname) {
  tt <- horizon * 365.25
  keep <- !is.na(time) & !is.na(event) & time > 0
  time <- time[keep]; event <- event[keep]; marker <- marker[keep]
  r <- timeROC(T = time, delta = event, marker = marker, cause = 1,
               times = tt, ROC = TRUE, iid = TRUE)
  # 关键：timeROC 在「单时间点」时会自动插入 t=0 列（列名 "t=0"），
  # 必须按列名 t=... 匹配目标时间点，不能用列号硬取（否则取到 t=0 退化线）。
  cn <- colnames(r$FP)
  if (is.null(cn)) cn <- paste0("t=", r$times)
  idx <- which.min(abs(as.numeric(sub("^t=", "", cn)) - tt))
  FP <- if (is.matrix(r$FP)) r$FP[, idx] else r$FP
  TP <- if (is.matrix(r$TP)) r$TP[, idx] else r$TP
  auc <- as.numeric(r$AUC[idx])
  # confint 返回 list（$CI_AUC：行=时间点、列=2.5%/97.5%、值域 0-100）
  ci <- tryCatch({
    ci_all <- confint(r, level = 0.95)$CI_AUC
    j <- which.min(abs(as.numeric(sub("^t=", "", rownames(ci_all))) - tt))
    as.numeric(ci_all[j, ]) / 100
  }, error = function(e) c(NA, NA))
  df <- data.frame(FP = FP, TP = TP)
  lbl <- if (all(!is.na(ci))) sprintf("AUC = %.3f\n(95%% CI %.3f-%.3f)", auc, ci[1], ci[2]) else
    sprintf("AUC = %.3f", auc)
  ti <- sprintf("%s - %dy %s", cohort, horizon, endpoint)
  st <- sprintf("n = %d, events = %d", length(time), sum(event))
  p <- ggplot(df, aes(FP, TP)) +
    geom_ribbon(aes(ymin = 0, ymax = TP), fill = col, alpha = 0.13) +
    geom_path(color = col, linewidth = 1.15) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    annotate("label", x = 0.56, y = 0.075, label = lbl, size = 4.1,
             hjust = 0, vjust = 0, lineheight = 0.95, colour = col,
             fontface = "bold", fill = "white", label.r = unit(0.12, "lines")) +
    labs(title = ti, subtitle = st,
         x = "1 - Specificity (false positive rate)",
         y = "Sensitivity (true positive rate)") +
    theme_bw(base_size = 12.5) +
    theme(plot.title = element_text(face = "bold", size = 13.5, colour = col),
          plot.subtitle = element_text(colour = "grey30", size = 11),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUTD, fname), p, width = 5.4, height = 5.2, dpi = 300)
  logmsg(sprintf("%-46s AUC=%.3f  CI[%s]", fname, auc,
                 paste(sprintf("%.3f", ci), collapse = ", ")))
  invisible(data.frame(cohort = cohort, horizon = horizon, endpoint = endpoint,
                       n = length(time), events = sum(event),
                       AUC = auc, CI_low = ci[1], CI_high = ci[2]))
}

# ---------- 单队列 AUC-time 折线 ----------
plot_auc_curve <- function(auc_df, cohort, col, fname) {
  p <- ggplot(auc_df, aes(horizon, AUC)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "grey55") +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = col, alpha = 0.15) +
    geom_line(color = col, linewidth = 1.0) +
    geom_point(color = col, size = 3.0) +
    geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -1.1, size = 3.5,
              colour = col, fontface = "bold") +
    scale_x_continuous(breaks = auc_df$horizon) +
    coord_cartesian(ylim = c(0.45, 1.0)) +
    labs(title = paste0(cohort, " - time-dependent AUC"),
         x = "Time (years)", y = "AUC") +
    theme_bw(base_size = 12.5) +
    theme(plot.title = element_text(face = "bold", size = 13.5, colour = col),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUTD, fname), p, width = 5.6, height = 4.6, dpi = 300)
  logmsg(paste("AUC curve written:", fname))
}

all_rows <- data.frame()

# ---------- 队列 1：TARGET-OS 训练集（EFS 1/2/3/4y） ----------
d1 <- fread(file.path(TM, "P6_riskscore_TARGET_train.csv"))
h1 <- c(1, 2, 3, 4)
for (h in h1) {
  all_rows <- rbind(all_rows, plot_single_roc(
    d1$EFS_time, d1$EFS_event, d1$RiskScore, h,
    "TARGET-OS training cohort", "EFS", COL["TARGET-OS"],
    sprintf("P7_ROC_TARGET_train_%dy.png", h)))
}
a1 <- all_rows[all_rows$cohort == "TARGET-OS training cohort", ]
plot_auc_curve(a1, "TARGET-OS training cohort (n=85, 41 EFS events)", COL["TARGET-OS"],
               "P7_AUCcurve_TARGET_train.png")

# ---------- 队列 2：GSE21257（OS 2/3/5y） ----------
d2 <- fread(file.path(TM, "P6_riskscore_GSE21257.csv"))
d2 <- d2[!is.na(OS_time) & !is.na(OS_event) & OS_time > 0]
h2 <- c(2, 3, 5)
n2 <- nrow(all_rows)
for (h in h2) {
  all_rows <- rbind(all_rows, plot_single_roc(
    d2$OS_time, d2$OS_event, d2$RiskScore, h,
    "GSE21257 external cohort", "OS", COL["GSE21257"],
    sprintf("P7_ROC_GSE21257_%dy.png", h)))
}
a2 <- all_rows[(n2 + 1):nrow(all_rows), ]
plot_auc_curve(a2, "GSE21257 external cohort (n=53, 23 OS events)", COL["GSE21257"],
               "P7_AUCcurve_GSE21257.png")

# ---------- 队列 3：GSE39055（EFS 1/2/3/5y） ----------
d3 <- fread(file.path(TM, "P6c_riskscore_GSE39055.csv"))
d3 <- d3[!is.na(EFS_time) & !is.na(EFS_event) & EFS_time > 0]
h3 <- c(1, 2, 3, 5)   # 4y 与 3y AUC 相同（数据特性）
n3 <- nrow(all_rows)
for (h in h3) {
  all_rows <- rbind(all_rows, plot_single_roc(
    d3$EFS_time, d3$EFS_event, d3$RiskScore, h,
    "GSE39055 external cohort", "EFS", COL["GSE39055"],
    sprintf("P7_ROC_GSE39055_%dy.png", h)))
}
a3 <- all_rows[(n3 + 1):nrow(all_rows), ]
plot_auc_curve(a3, "GSE39055 external cohort (n=36, 18 EFS events)", COL["GSE39055"],
               "P7_AUCcurve_GSE39055.png")

fwrite(all_rows, file.path(TM, "P7_roc_individual_summary.csv"))
logmsg("=== P7a2 complete: ", nrow(all_rows), " individual ROC panels ===")
