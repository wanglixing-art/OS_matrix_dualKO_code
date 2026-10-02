# ============================================================
# P7a: 三队列 time-ROC 曲线（用户指令：训练集与测试集都出 ROC）
#
# 队列与打分（全部沿用 P6/P6c 口径，不重拟合）：
#   1) TARGET-OS 训练集（n=85, EFS）— P6_riskscore_TARGET_train.csv
#   2) GSE21257（n=53, OS）— P6_riskscore_GSE21257.csv
#   3) GSE39055（n=36, EFS）— P6c_riskscore_GSE39055.csv
# 画法：timeROC(ROC=TRUE) 多时间点叠画一张图/队列，AUC 图例内嵌；
#       每队列独立配色（红-橙-蓝绿渐变，符合可发表彩色风格）
# ============================================================
suppressMessages({
  library(data.table); library(survival); library(timeROC)
  library(ggplot2); library(patchwork)
})
PROJ <- "D:/projects/OS_matrix_dualKO"
TM <- file.path(PROJ, "results/tables"); FG <- file.path(PROJ, "results/figures")
LOG <- file.path(PROJ, "logs/P7a_run.log")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = " "))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE)
}

# ---------- 通用绘制函数 ----------
plot_timeROC <- function(time, event, marker, horizons, title, unit_days = TRUE) {
  tt <- horizons * if (unit_days) 365.25 else 365.25
  keep <- !is.na(time) & !is.na(event) & time > 0
  time <- time[keep]; event <- event[keep]; marker <- marker[keep]
  r <- timeROC(T = time, delta = event, marker = marker, cause = 1,
               times = tt, ROC = TRUE)
  # r$FP / r$TP 是 (阈值点 × 时间) 矩阵，按列取（已实测：dim 86×4）
  # color = 时间点（含 AUC），保证每个 horizon 独立颜色（AUC 相同也不合并）
  df <- do.call(rbind, lapply(seq_along(tt), function(i) {
    data.frame(FP = r$FP[, i], TP = r$TP[, i],
               horizon = paste0(horizons[i], "y (AUC ",
                                sprintf("%.3f", r$AUC[i]), ")"))
  }))
  lv <- unique(df$horizon)
  pal <- setNames(c("#D7301F", "#F16913", "#1F78B4", "#33A02C", "#6A3D9A"), lv)
  p <- ggplot(df, aes(FP, TP, color = horizon)) +
    geom_path(linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
    scale_color_manual(values = pal, name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "1 - Specificity (false positive rate)",
         y = "Sensitivity (true positive rate)") +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.72, 0.22),
          legend.background = element_rect(fill = "white", color = "grey75",
                                           linewidth = 0.3),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = 13))
  p
}

# ---------- 队列 1：TARGET-OS 训练集（EFS，1/2/3/4 年） ----------
d1 <- fread(file.path(TM, "P6_riskscore_TARGET_train.csv"))
p1 <- plot_timeROC(d1$EFS_time, d1$EFS_event, d1$RiskScore, c(1, 2, 3, 4),
                   paste0("TARGET-OS training cohort (n=", nrow(d1),
                          ", ", sum(d1$EFS_event), " events), EFS endpoint"))
ggsave(file.path(FG, "P7_timeROC_TARGET_train.png"), p1, width = 6.4, height = 6.0, dpi = 300)
auc1v <- timeROC(T = d1$EFS_time, delta = d1$EFS_event, marker = d1$RiskScore,
                 cause = 1, times = c(1, 2, 3, 4) * 365.25)$AUC
logmsg("TARGET-OS timeROC done: AUC 1-4y = ",
       paste(round(auc1v, 3), collapse = ", "))

# ---------- 队列 2：GSE21257（OS，2/3/5 年） ----------
d2 <- fread(file.path(TM, "P6_riskscore_GSE21257.csv"))
d2 <- d2[!is.na(OS_time) & !is.na(OS_event) & OS_time > 0]
p2 <- plot_timeROC(d2$OS_time, d2$OS_event, d2$RiskScore, c(2, 3, 5),
                   paste0("GSE21257 external cohort (n=", nrow(d2),
                          ", ", sum(d2$OS_event), " events), OS endpoint"))
ggsave(file.path(FG, "P7_timeROC_GSE21257.png"), p2, width = 6.4, height = 6.0, dpi = 300)
logmsg("GSE21257 timeROC done")

# ---------- 队列 3：GSE39055（EFS，1-5 年） ----------
d3 <- fread(file.path(TM, "P6c_riskscore_GSE39055.csv"))
d3 <- d3[!is.na(EFS_time) & !is.na(EFS_event) & EFS_time > 0]
h3 <- c(1, 2, 3, 5)   # 4y 与 3y AUC 相同（数据特性），取 1/2/3/5
p3 <- plot_timeROC(d3$EFS_time, d3$EFS_event, d3$RiskScore, h3,
                   paste0("GSE39055 external cohort (n=", nrow(d3),
                          ", ", sum(d3$EFS_event), " events), EFS endpoint"))
ggsave(file.path(FG, "P7_timeROC_GSE39055.png"), p3, width = 6.4, height = 6.0, dpi = 300)
logmsg("GSE39055 timeROC done")

# ---------- 汇总图（三联版，供手稿复合图用） ----------
pal_all <- c("TARGET-OS (EFS)" = "#D7301F", "GSE21257 (OS)" = "#1F78B4",
             "GSE39055 (EFS)" = "#33A02C")
auc_rows <- data.frame()
mk_auc <- function(time, event, marker, horizons, cohort) {
  for (yr in horizons) {
    r <- timeROC(T = time, delta = event, marker = marker, cause = 1,
                 times = yr * 365.25)
    auc_rows <<- rbind(auc_rows, data.frame(cohort = cohort, year = yr,
                                            AUC = r$AUC[2]))
  }
}
mk_auc(d1$EFS_time, d1$EFS_event, d1$RiskScore, c(1,2,3,4), "TARGET-OS (EFS)")
mk_auc(d2$OS_time, d2$OS_event, d2$RiskScore, c(2,3,5), "GSE21257 (OS)")
mk_auc(d3$EFS_time, d3$EFS_event, d3$RiskScore, c(1,2,3,5), "GSE39055 (EFS)")
fwrite(auc_rows, file.path(TM, "P7_timeAUC_all_cohorts.csv"))

pa <- ggplot(auc_rows, aes(factor(year), AUC, group = cohort, color = cohort)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey55") +
  geom_line(linewidth = 0.55, alpha = 0.8,
            position = position_dodge(0.15)) +
  geom_point(size = 2.6, position = position_dodge(0.15)) +
  scale_color_manual(values = pal_all, name = NULL) +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(title = "Time-dependent AUC of the 14-gene LASSO-Cox risk score",
       subtitle = "Training (TARGET-OS, EFS) and two external cohorts (GSE21257 OS; GSE39055 EFS)",
       x = "Time (years)", y = "AUC") +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(FG, "P7_timeAUC_all_cohorts.png"), pa, width = 7.6, height = 5.0, dpi = 300)
logmsg("=== P7a complete ===")
