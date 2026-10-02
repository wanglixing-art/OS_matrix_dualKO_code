#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""P8b_rewrite_methods.py — 重写材料与方法：补全软件版本、真实参数、超几何检验背景"""
import re, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
MD = r"D:\projects\OS_matrix_dualKO\manuscript\OS_manuscript_v1_cn.md"
t = open(MD, encoding="utf-8").read()

NEW = """## 材料与方法

### 数据获取与预处理

骨肉瘤批量 RNA-seq 数据通过 NCI 基因组数据共享平台（GDC）API 自 TARGET-OS 项目获取。下载 88 个单样本 STAR 计数文件（GENCODE v36 注释，非链特异性）并合并为基因×样本计数矩阵；临床实体表（383 例，含随访扩展）取自同一来源。经蛋白编码基因聚合与低表达过滤（edgeR 4.10.5 的 filterByExpr）后保留 13,738 个基因，采用 TMM（trimmed mean of M-values）标准化并经 limma-voom（limma 3.68.5）转换。无事件生存（EFS）与总生存（OS）的时间与事件指示变量取自 GDC 官方临床实体；88 个表达样本中 85 例具可评估 EFS（41 个事件）、86 例具可评估 OS（29 个事件），以 EFS 为主要终点。

公共微阵列与 RNA-seq 队列用于外部验证与差异表达分析：GSE42352（84 例肿瘤活检、15 例正常间充质干细胞/成骨细胞，用于差异表达分析）、GSE21257（n = 53，化疗前活检，与 GSE33382 同平台）、GSE39055（n = 37，EFS 终点，含复发状态、化疗方案与坏死率）、GSE33382（n = 53 例活检，5 年内转移）与 GSE87624（n = 52，原发灶与转移灶 RNA-seq）。系列矩阵文件经 GEOquery 2.80.0 直接解析，临床特征自 characteristics_ch1 字段提取，而非依赖文献描述。探针注释对 GSE21257/GSE33382 使用 GPL10295 平台表，对 GSE39055 使用 illuminaHumanv4.db 包（完整映射 14 个签名基因）。单细胞 RNA-seq 数据取自 GSE162454（6 个 10x Genomics 样本）。

### 差异表达与加权基因共表达网络分析

使用 limma 评估肿瘤活检与正常对照之间的差异表达。基因以 |log2 倍数变化| > 0.585（对应倍数变化 1.5）定义为上调，共获得 556 个上调基因与 484 个下调基因。加权基因共表达网络分析（WGCNA 1.74）基于中位绝对偏差排名前 8,000 的基因，采用双权重中相关（bicor）与无符号网络类型，软阈值幂经 pickSoftThreshold 确定为 4（无标度拓扑拟合 R² = 0.929）；blockwiseModules 的参数为 minModuleSize = 30、mergeCutHeight = 0.25、maxBlockSize = 8,000，共产生 21 个模块。模块与 EFS/OS 事件状态的性状相关性分析识别出 7 个与结局显著相关的模块（p < 0.05，共 2,553 个基因），其中 ME3 与 EFS 事件呈最强正相关（r = 0.396，p = 1.34e-4），ME2 呈最强负相关（r = -0.251，p = 0.0185）。肿瘤上调差异表达基因与结局相关模块基因的交集定义了包含 299 个基因的候选池。

### 蛋白质-蛋白质相互作用与枢纽基因筛选

将候选池映射至 STRING 数据库（置信度阈值 0.7），得到含 166 个节点与 460 条边的网络。使用 cytoHubba 的度（degree）、MCC、MNC、EPC 与中介中心性（betweenness）量化网络中心性，并与针对 EFS 的单因素 Cox 回归整合。在 85 个训练样本上比较多种机器学习算法的表现（glmBoost、Lasso、弹性网络、Ridge、StepCox），预选 30 个 Cox p < 0.1 的特征用于事件预测性能比较。最终选择代表不同功能区室的两个枢纽基因。

### 生存分析与外部验证

在 TARGET-OS 训练集中开展 Kaplan-Meier 曲线（按中位表达分组）、对数秩检验、单因素 Cox 回归（每标准差，per-SD）以及校正年龄与性别的多因素 Cox 回归（survival 3.8.6、survminer 0.5.2）。在 GSE21257 中检验单基因与转移状态的关联，并在 GSE42352 中重现肿瘤与正常组织的表达差异。

### 单细胞 RNA-seq 处理与细胞类型注释

将 GSE162454 的 6 个样本读入为 10x 矩阵（50,780 个细胞），经 nFeature_RNA 200-6,000 与线粒体含量 < 20% 过滤后保留 42,784 个细胞（Seurat 5.5.1，CreateSeuratObject 参数 min.cells = 3、min.features = 200）。标准化后选取 3,000 个高变基因，缩放并进行主成分分析（40 个主成分），随后以 Harmony 2.0.5 进行跨样本整合。基于 harmony 嵌入的前 30 个主成分构建近邻图，在分辨率 1.0 下聚类得到 38 个簇，并经标志物分配与基于两个独立参考集的 SingleR 2.14.2 注释为 9 种细胞类型。使用 inferCNV 1.28.0（HMM i6 模型，analysis_mode = "subclusters"，cutoff = 0.1，denoise = TRUE）判定恶性细胞：参考细胞（T 细胞、NK 细胞、B 细胞、髓系细胞、内皮细胞与破骨细胞，降采样至 4,612 个）与推定恶性区室（成骨细胞样与增殖细胞，6,000 个）在 8,066 个基因上进行比较，CNV 评分超过参考区室第 95 百分位（阈值 0.050）的细胞被判为恶性；参考区室自身产生 5.0% 的恶性率，用于确认阈值校准。

### 细胞-细胞通讯分析

将 CellChat 2.2.0.9001 应用于注释数据集，使用 CellChatDB.human 在细胞类型水平分析，以 100 次置换识别显著信号通路与配体-受体对，并以 filterCommunication（min.cells = 10）过滤低丰度细胞群。为评估稳健性，另在 6 个样本中各自独立推断通讯，随后经 liftCellChat 协调与 mergeCellChat 聚合（128 条通路的并集）。

### scTenifoldKnk 虚拟敲除

使用 scTenifoldKnk 1.1 开展虚拟敲除。高变基因按 log2(CPM+1) 的方差定义，保留前 1,500 个并强制纳入两个枢纽基因（最终 1,501 个基因）。每次敲除从目标区室中采样 1,000 个细胞，按每 300 个细胞构建主成分网络（nc_nComp = 3）；张量经 CP 分解（td_K = 3，最大迭代 1,000 次，误差 1e-5）分解，基因调控网络经流形对齐（d = 2）比对。差异调控由 dRegulation 与经验零分布计算。由于默认卡方零分布在该网络规模下过于严格、且 251 个基因的距离值下溢至约 1e-17（浮点伪影），失调基因（DRG）定义为单向 Z > 2 且距离 ≥ 1e-10。DRG 与候选池的重叠采用超几何检验评估，背景为虚拟敲除网络中的 1,501 个基因。

### 特征构建与验证

首先对候选池应用单因素 Cox 回归，299 个基因中有 141 个 p < 0.1 者进入 LASSO 惩罚 Cox 模型（glmnet 5.0，family = "cox"，alpha = 1，10 折交叉验证，type.measure = "deviance"，种子 2026）。lambda.min 解保留 14 个基因（lambda.1se 解保留基因数不足 5 个，故采用 lambda.min）。风险评分计算为各基因标准化表达（数据集内 z 分数）与其训练系数的加权和；外部验证中系数从不重新拟合。区分度以 Harrell 一致性指数（C-index）评估，时间依赖 AUC 由 timeROC 0.4.1 在 1、2、3、4、5 年计算，独立性以多因素 Cox 回归评估。敏感性分析在 OS 终点上重新训练模型，以检验终点不匹配的贡献。

### 微环境与肿瘤突变负荷

在 illumina 平台上对 85 个训练样本计算基质、免疫与综合 ESTIMATE 评分（estimate 1.0.13），并与风险评分进行相关性分析。TARGET-OS 体细胞突变数据自 GDC 获取，为 167 个掩码体细胞突变 MAF 文件（142 例）。非沉默变异限定为 9 类（错义、无义、移码缺失/插入、框内缺失/插入、剪接位点、翻译起始位点与非终止突变），在同一病例内按染色体位置、参考等位基因与替代等位基因去重，按 38 Mb 外显子组足迹估算肿瘤突变负荷（TMB）。

### 统计分析

全部分析在 R 4.6.1 中进行。连续变量以 Wilcoxon 秩和检验与 Spearman 相关性比较，生存差异以对数秩检验评估，富集检验采用超几何分布。双侧 p < 0.05 视为具有统计学意义。随机过程均固定随机种子（主种子 20260928；LASSO 交叉验证种子 2026），全部中间结果表均保留以确保可重复性。
"""

m = re.search(r"## 材料与方法\n.*?(?=\n!\[图 1｜)", t, flags=re.S)
assert m, "材料与方法节未找到"
t = t[:m.start()] + NEW + "\n" + t[m.end():]
open(MD, "w", encoding="utf-8").write(t)

cn = len(re.findall(r"[\u4e00-\u9fff]", NEW))
subs = len(re.findall(r"^### ", NEW, flags=re.M))
print(f"材料与方法已重写: {cn} 中文汉字, {subs} 小节")
print("版本号抽查:", [x for x in ["edgeR 4.10.5","limma 3.68.5","WGCNA 1.74","Seurat 5.5.1","Harmony 2.0.5","SingleR 2.14.2","inferCNV 1.28.0","CellChat 2.2.0.9001","scTenifoldKnk 1.1","glmnet 5.0","timeROC 0.4.1","survival 3.8.6","estimate 1.0.13"] if x in NEW])
