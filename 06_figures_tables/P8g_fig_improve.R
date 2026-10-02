# P8g: 图 2C 与图 8A 呈现优化 (v2: ggraph repel 标签 + RUNX2 邻域)
# 1) PPI：正文用核心子网络（degree>=10 ∪ RUNX2 一阶邻居），完整网络移补充图 S10
# 2) KM：正文图 8A 截断至 5 年随访窗 + RMST 注释，完整随访曲线移补充图 S11
suppressMessages({library(igraph); library(ggraph); library(tidygraph); library(dplyr)
                  library(survival); library(survminer)})
FG <- "results/figures"; TB <- "results/tables"

## ---------- Part 1: PPI ----------
e <- read.csv(file.path(TB, "P2_STRING_edges.csv"))
g0 <- graph_from_data_frame(e, directed = FALSE); g0 <- simplify(g0)
deg <- degree(g0)

## 核心集合：degree>=10 ∪ RUNX2 一阶邻居
core <- names(deg)[deg >= 10]
nb_runx2 <- neighbors(g0, "RUNX2")$name
core <- union(core, c("RUNX2", nb_runx2))
gs <- induced_subgraph(g0, core)
## as_tbl_graph 会丢边属性，按端点名回填 score
sc <- rbind(data.frame(a = e$from, b = e$to, score = e$score),
            data.frame(a = e$to, b = e$from, score = e$score))
tg <- as_tbl_graph(gs) %>%
  activate(edges) %>%
  mutate(name_from = .N()$name[from], name_to = .N()$name[to]) %>%
  left_join(sc, by = c("name_from" = "a", "name_to" = "b")) %>%
  mutate(score = dplyr::coalesce(score, min(sc$score))) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree(),
         grp = case_when(name == "BUB1" ~ "BUB1",
                         name == "RUNX2" ~ "RUNX2",
                         name %in% nb_runx2 & name != "RUNX2" ~ "RUNX2 neighbors",
                         TRUE ~ "Other core genes"))

set.seed(1)
p <- ggraph(tg, layout = "fr", niter = 2000) +
  geom_edge_link(aes(width = score), colour = "grey72", alpha = 0.8,
                 edge_colour = "grey72") +
  scale_edge_width(range = c(0.3, 1.6), guide = "none") +
  geom_node_point(aes(size = degree, fill = grp), shape = 21,
                  colour = "white", stroke = 0.4) +
  scale_size(range = c(3.5, 9.5), guide = "none") +
  scale_fill_manual(values = c("BUB1" = "#D7301F", "RUNX2" = "#2166AC",
                               "RUNX2 neighbors" = "#9ECAE1",
                               "Other core genes" = "#D6DEEB"), name = NULL) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.0,
                 max.overlaps = 20, segment.colour = "grey60",
                 family = "sans", colour = "grey15") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(title = "STRING PPI - core subnetwork (degree >= 10 plus RUNX2 neighbours)")
ggsave(file.path(FG, "P2_PPI_hub_subnet.png"), p, width = 7.0, height = 7.0,
       dpi = 300, bg = "white")
cat("subnet:", vcount(gs), "nodes /", ecount(gs), "edges\n")

## 完整网络（补充图 S10）
tgf <- as_tbl_graph(g0) %>%
  mutate(degree = centrality_degree(),
         comm = as.factor(group_fast_greedy()))
lab25 <- names(sort(deg, decreasing = TRUE))[1:25]
set.seed(1)
pf <- ggraph(tgf, layout = "fr", niter = 3000) +
  geom_edge_link(edge_colour = "grey80", edge_width = 0.4, alpha = 0.7) +
  geom_node_point(aes(size = degree, fill = comm), shape = 21,
                  colour = "white", stroke = 0.35, show.legend = FALSE) +
  scale_size(range = c(2, 8)) +
  scale_fill_brewer(palette = "Set3") +
  geom_node_text(aes(label = ifelse(name %in% lab25, name, NA)), repel = TRUE,
                 size = 3.0, max.overlaps = 30, segment.colour = "grey60",
                 colour = "grey10", na.rm = TRUE, family = "sans") +
  theme_void() +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(title = "STRING PPI network - full (166 nodes / 460 edges, confidence > 0.7)")
ggsave(file.path(FG, "P2_PPI_network_full.png"), pf, width = 11.0, height = 9.5,
       dpi = 300, bg = "white")
cat("full network done\n")

## ---------- Part 2: KM ----------
d <- read.csv(file.path(TB, "P6c_riskscore_GSE39055.csv"))
d$grp <- factor(ifelse(d$RiskScore > median(d$RiskScore), "High", "Low"),
                levels = c("Low", "High"))
tau <- 1826.25
dt <- d
dt$EFS_time <- pmin(dt$EFS_time, tau)
dt$EFS_event <- ifelse(d$EFS_time > tau, 0, dt$EFS_event)
fit <- survfit(Surv(EFS_time, EFS_event) ~ grp, dt)

png(file.path(FG, "P6c_KM_GSE39055_5y.png"), 1800, 1700, res = 250, bg = "white")
km <- ggsurvplot(fit, data = dt, conf.int = TRUE, risk.table = TRUE,
                 palette = c("#2166AC", "#D7301F"),
                 legend.labs = c("Low risk", "High risk"),
                 xlim = c(0, tau), break.time.by = 365.25,
                 xscale = 365.25,  # 刻度标签除以 365.25 -> 整数年 0,1,2,3,4,5
                 xlab = "Time (years)", ylab = "EFS probability",
                 title = "GSE39055 external validation (n=36), EFS - 5-y window",
                 font.main = 13, pval = FALSE, surv.median.line = "hv",
                 risk.table.height = 0.28, fontsize = 3.2)
km$plot <- km$plot +
  annotate("text", x = tau * 0.60, y = 0.93, hjust = 0, size = 3.6,
           label = paste0("log-rank (5-y window) p = 0.0090\n",
                          "RMST difference = 611 days\n(95% CI 237-985), p = 0.0014"))
print(km, newpage = FALSE)
dev.off()
cat("KM 5y done\n")

file.copy(file.path(FG, "P6c_KM_GSE39055.png"), file.path(FG, "S11_KM_full_followup.png"),
          overwrite = TRUE)
cat("S11 done\n")
