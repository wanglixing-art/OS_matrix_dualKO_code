#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""P7c_split_figures.py — 图表分级重构：显著/美观图留正文，弱图/技术图移补充材料"""
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

MD = r"D:\projects\OS_matrix_dualKO\manuscript\OS_manuscript_v1_cn.md"
t = open(MD, encoding="utf-8").read()
orig = t

# ---------------------------------------------------------------
# 1) 图注整行替换（图 2/3/4/6/9/10）
NEW_CAPTIONS = {
"图 2｜": "**图 2｜候选池构建与双枢纽基因鉴定。** (A) GSE42352 肿瘤 vs 正常差异表达火山图（556 上调 / 484 下调，|log2FC| > 0.585）；(B) 299 基因候选池的基因本体（GO）生物学过程富集点图；(C) STRING 蛋白互作网络（置信度 > 0.7，166 节点 / 460 边，含细胞周期主聚类、胶原/ECM 亚聚类与免疫聚类）；(D) RUNX2 按中位表达分组的 Kaplan-Meier EFS 曲线（TARGET-OS，n = 85，log-rank p = 0.014）。WGCNA 软阈值筛选与模块-性状相关性分析见补充图 S1；BUB1 的单基因生存分析（中位分组 log-rank p > 0.1 而 per-SD Cox p = 0.0385）见补充图 S2。",
"图 3｜": "**图 3｜单细胞图谱与双枢纽基因的区室定位。** (A) UMAP 细胞类型注释（6 个 GSE162454 样本，质控后 42,784 个细胞；9 种类型：髓系 33.6%、CAF/MSC 19.1%、增殖 13.5%、T 细胞 12.3%、成骨 8.5%、破骨 5.9%、B 细胞 3.6%、内皮 2.1%、NK 1.4%）；(B) 谱系标志基因在各聚类中的表达点图；(C) BUB1 与 RUNX2 的 FeaturePlot 表达梯度；(D) 双枢纽在各细胞类型中的表达小提琴图（RUNX2 在 CAF/MSC 均值 0.921、成骨细胞 0.804；BUB1 在增殖细胞 0.216、恶性细胞 0.209）。质控前后对比、UMAP 聚类与样本来源分布见补充图 S3。",
"图 4｜": "**图 4｜inferCNV 恶性判定与细胞通讯架构。** (A) 各细胞类型恶性比例汇总（成骨细胞 84.2%、增殖细胞 49.2%、参照区室 5.0%）；(B) 双枢纽表达在 CNV 定义恶性 vs 非恶性细胞中的比较（RUNX2 恶性中位 0.778 vs 非恶性 0，p ≈ 0；BUB1 p = 5.8e-35）；(C) CellChat 通讯网络环形图（相互作用权重；CAF/MSC 发出权重 8.34 为全网络最高）；(D) 信号通路权重热图（COLLAGEN 平均权重 11.31，为第二名通路的 2.3 倍）；(E) 逐样本的枢纽区室发出信号量（CAF/MSC 在 6/6 样本中排名第一，平均权重 11.23，CV = 0.385）；(F) 基质至增殖细胞的发送比值（6/6 样本 > 1，范围 1.18-3.33）。inferCNV 观测热图（参照区室 95 分位阈值，thr = 0.050）见补充图 S4。",
"图 6｜": "**图 6｜双基因敲除特异性下游程序的 GO 富集对比。** (A) RUNX2-KO 特异性差异调控基因（DRG，20 个）的 GO 生物学过程富集（首位 TGF-β1 产生，3/20 基因，FDR = 2.15e-4）；(B) BUB1-KO 特异性 DRG（17 个）的 GO 富集（首位白细胞介导的免疫，8/17 基因，FDR = 9.1e-7）。全部差异调控基因的 GO 富集结果见补充图 S5。KEGG 通路水平的结果在正文中报告：RUNX2-KO 首位为 TGF-β 信号通路（z = 7.7），BUB1-KO 首位为 Fcγ 受体介导的吞噬作用（z = 9.8）。",
"图 9｜": "**图 9｜风险评分与转移表型的关联。** 按转移状态分组的风险评分箱线图（z-score 标度）：左，GSE33382 五年内转移 vs 未转移（中位 0.191 vs -0.547，p = 0.0017）；右，GSE87624 转移灶 vs 原发（中位 0.405 vs -0.203，p = 0.349，n = 9 转移灶，方向一致但效能不足）。GSE21257 中按转移状态分组的风险评分差异同样显著（中位 0.174 vs -0.582，p = 0.0021，见正文）；该队列 OS 终点的 Kaplan-Meier 分析（log-rank p = 0.613）见补充图 S6。",
"图 10｜": "**图 10｜风险评分与肿瘤微环境的关联（TARGET-OS）。** 按风险评分中位分组的 ESTIMATE 三评分箱线图（基质：中位 -77.5 vs 358.7，p = 1.95e-4，rho = -0.55；免疫：中位 -454.6 vs -173.9，p = 0.0089，rho = -0.328；综合：中位 -531.2 vs 341.5，p = 2.92e-4，rho = -0.492；n = 85）。风险评分与肿瘤突变负荷的关联（Spearman rho = 0.31，p = 0.0084；高危组中位 TMB 0.61 vs 低危组 0.47 突变/Mb，Wilcoxon p = 0.0255）见补充图 S7。",
}
lines = t.splitlines()
n_repl = 0
for i, L in enumerate(lines):
    for key, new in NEW_CAPTIONS.items():
        if L.startswith("**" + key):
            lines[i] = new
            n_repl += 1
t = "\n".join(lines)
print("图注替换:", n_repl, "/ 6")

# ---------------------------------------------------------------
# 2) 图片 alt 文本同步
ALT = [
    ("![图 6｜双 KO 下游程序富集对比]", "![图 6｜双 KO 特异性下游程序 GO 富集对比]"),
    ("![图 10｜微环境与突变负荷关联]", "![图 10｜微环境 ESTIMATE 关联]"),
]
for old, new in ALT:
    if old in t:
        t = t.replace(old, new)
        print("alt 更新:", new[:30])

# ---------------------------------------------------------------
# 3) 正文面板级引用更新
CITES = [
    # (旧, 新, 说明)
    ("（图 2B-C）", "（补充图 S1A-B）"),                       # WGCNA 两面板 -> 补充
    ("（图 2D-E）", "（图 2B-C）"),                            # GO+PPI 新字母
    ("p = 0.0176；图 2F-G）", "p = 0.0176；图 2D；BUB1 见补充图 S2）"),
    ("（图 4A-C）", "（图 4A-B；inferCNV 观测热图见补充图 S4）"),
    ("（图 4D-E）", "（图 4C-D）"),
    ("（图 4F-G）", "（图 4E-F）"),
    ("（图 5A-C）", "（图 5A-C）"),                            # 不变，占位
    ("（图 6A-B）", "（图 6A）"),
    ("（图 6C-D）", "（图 6B；全部 DRG 见补充图 S5）"),
    ("（图 7A）", "（图 7A）"),                                # 不变，占位
    ("（图 7B-C）", "（图 7B-C）"),                            # 不变，占位
    ("（图 10A）", "（图 10A）"),                              # 不变，占位
    ("（图 10B）", "（补充图 S7）"),
]
for old, new in CITES:
    cnt = t.count(old)
    if old == new:
        print("占位不变:", old, cnt)
        continue
    if cnt == 0:
        print("!! 未命中:", old)
    else:
        t = t.replace(old, new)
        print("引用更新:", old, "->", new[:40], f"(x{cnt})")

# ---------------------------------------------------------------
# 4) 图表清单节整体替换
m = re.search(r"## 图表清单\n.*?(?=\n## 参考文献)", t, flags=re.S)
assert m, "图表清单节未找到"
NEW_LIST = """## 图表清单

### 正文图（Figure 1-10）

| 编号 | 内容 | 源文件（results/figures/） |
|---|---|---|
| 图 1 | 研究设计与分析流程总览 | F1.pdf（源码 code/P7b_fig1_design.R） |
| 图 2 | 候选池构建与双枢纽鉴定（4 面板） | F2.pdf |
| 图 3 | 单细胞图谱与双枢纽区室定位（4 面板） | F3.pdf |
| 图 4 | 恶性判定与 CellChat 通讯架构（6 面板） | F4.pdf |
| 图 5 | 双基因平行虚拟敲除扰动景观（3 面板） | F5.pdf |
| 图 6 | 双 KO 特异性下游程序 GO 富集对比（2 面板） | F6.pdf |
| 图 7 | 签名系数 + 训练集生存 + 时间依赖 ROC（3 面板） | F7.pdf |
| 图 8 | GSE39055 外部验证（3 面板） | F8.pdf |
| 图 9 | 风险评分与转移表型关联（1 面板） | F9.pdf |
| 图 10 | 微环境 ESTIMATE 评分关联（1 面板） | F10.pdf |

### 补充图（Figure S1-S9）

| 编号 | 内容 | 源文件（results/figures/） |
|---|---|---|
| 补充图 S1 | WGCNA 软阈值与模块-性状分析（2 面板） | S1.pdf |
| 补充图 S2 | 双枢纽单基因生存分析（3 面板） | S2.pdf |
| 补充图 S3 | 单细胞质控与聚类概览（4 面板） | S3.pdf |
| 补充图 S4 | inferCNV 观测热图 | S4.pdf |
| 补充图 S5 | 双 KO 全部 DRG 的 GO 富集（2 面板） | S5.pdf |
| 补充图 S6 | GSE21257 队列签名 Kaplan-Meier OS | S6.pdf |
| 补充图 S7 | 风险评分与 TMB 关联 | S7.pdf |
| 补充图 S8 | 逐时间点独立 ROC 曲线含 95% CI（11 面板） | S8.pdf |
| 补充图 S9 | 各队列 AUC-时间曲线（3 面板） | S9.pdf |

| 表 | 内容 |
|---|---|
| 表 1 | 研究队列与临床特征 |
| 表 2 | 14 基因 LASSO-Cox EFS 签名（系数/方向/功能归属） |
| 表 3 | 风险评分的多因素独立性分析 |
| 表 4 | 训练与四个独立队列的表现汇总 |"""
t = t[:m.start()] + NEW_LIST + t[m.end():]
print("图表清单已更新")

# ---------------------------------------------------------------
# 5) 插入补充材料节（参考文献之前）
SUPPL = """## 补充材料图注（Supplementary Figures）

**补充图 S1｜WGCNA 共表达网络构建细节。** (A) 软阈值筛选（幂 = 4，无标度拓扑 R² = 0.929）；(B) 模块-性状相关性热图（21 个模块，其中 7 个与 EFS/OS 事件显著相关）。

![补充图 S1｜WGCNA 构建细节](../results/figures/S1.png)

**补充图 S2｜双枢纽基因的单基因生存分析（TARGET-OS）。** (A) BUB1 与 EFS（中位分组 log-rank p = 0.126；per-SD HR = 1.39，95% CI 1.02-1.89，p = 0.039）；(B) BUB1 与 OS（log-rank p = 0.739）；(C) RUNX2 与 OS（log-rank p = 0.056；per-SD HR = 1.49，95% CI 0.98-2.26，p = 0.060）。RUNX2 与 EFS 的显著关联见正文图 2D。

![补充图 S2｜双枢纽单基因生存分析](../results/figures/S2.png)

**补充图 S3｜单细胞质控与聚类概览（GSE162454）。** (A) 质控前小提琴图；(B) 质控后小提琴图（42,784 个细胞）；(C) UMAP 聚类（38 簇）；(D) UMAP 样本来源（6 个样本）。

![补充图 S3｜单细胞质控与聚类概览](../results/figures/S3.png)

**补充图 S4｜inferCNV 观测热图。** 参照区室 95 分位阈值（thr = 0.050）下的 CNV 观测热图；基于 CNV 的恶性判定共判别 10,612 个恶性细胞（各细胞类型恶性比例见正文图 4A）。

![补充图 S4｜inferCNV 观测热图](../results/figures/S4.png)

**补充图 S5｜双 KO 全部差异调控基因的 GO 富集。** (A) RUNX2-KO 全部差异调控基因的 GO 生物学过程富集（top 10 按 FDR）；(B) BUB1-KO 全部差异调控基因的 GO 富集（top 10 按 FDR）。特异性 DRG 的富集结果见正文图 6。

![补充图 S5｜双 KO 全部 DRG 的 GO 富集](../results/figures/S5.png)

**补充图 S6｜GSE21257 队列中按风险评分分组的 Kaplan-Meier OS 曲线。** n = 55（低危 27 / 高危 28），log-rank p = 0.613；该队列 OS 终点整体不显著，但按转移状态分组的风险评分差异显著（中位 0.174 vs -0.582，p = 0.0021，见正文）。

![补充图 S6｜GSE21257 Kaplan-Meier OS](../results/figures/S6.png)

**补充图 S7｜风险评分与肿瘤突变负荷的关联（TARGET-OS）。** n = 72（非同义突变，掩码 MAF ≥ 38 Mb），Spearman rho = 0.31，p = 0.0084；高危组中位 TMB 0.61 vs 低危组 0.47 突变/Mb，Wilcoxon p = 0.0255。

![补充图 S7｜TMB 关联](../results/figures/S7.png)

**补充图 S8｜各队列逐时间点的独立时间依赖 ROC 曲线（含 95% CI）。** (A-D) TARGET-OS 训练集 1/2/3/4 年；(E-H) GSE39055 1/2/3/5 年；(I-K) GSE21257 2/3/5 年。汇总 AUC 与 95% CI 见正文图 7-8。

![补充图 S8｜逐时间点独立 ROC](../results/figures/S8.png)

**补充图 S9｜各队列 AUC-时间曲线。** (A) TARGET-OS 训练集；(B) GSE39055；(C) GSE21257。

![补充图 S9｜AUC-时间曲线](../results/figures/S9.png)

---

"""
t = t.replace("## 参考文献", SUPPL + "## 参考文献", 1)

open(MD, "w", encoding="utf-8").write(t)
print("写入完成, len:", len(orig), "->", len(t))
