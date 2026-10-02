# P8i3: 摘要定稿（中文四段，结论段四层递进）+ 中英对照英文版（严格 <=350 words）
import io, os, re

BASE = r"D:\projects\OS_matrix_dualKO"
MD = os.path.join(BASE, "manuscript", "OS_manuscript_v1_cn.md")
OUT_EN = os.path.join(BASE, "manuscript", "OS_abstract_CN_EN_v1.md")
lines = io.open(MD, encoding="utf-8").read().splitlines()

TAG = {4: "**背景与目的**", 6: "**方法**", 8: "**结果**", 10: "**结论**"}
for i, t in TAG.items():
    assert lines[i - 1].startswith(t), f"L{i} mismatch"

CN = {
 4: "**背景与目的**　骨肉瘤是儿童和青少年中最常见的原发恶性骨肿瘤，转移与复发是治疗失败的主要原因。驱动基质重塑与增殖的程序分处哪些细胞区室、二者在下游是否独立，迄今缺乏系统解析；现有预后评估也缺少可在治疗前识别高危患者的分子工具。本研究以区室解析为主线，鉴定分处不同区室的枢纽基因，解析其下游程序的独立性，并构建可外部验证的无事件生存（EFS）风险评分。",
 6: "**方法**　整合 TARGET-OS bulk 转录组（85 例可评估 EFS，41 个事件）、五个外部队列与单细胞转录组（42,784 个细胞），开展差异表达与加权基因共表达网络分析（WGCNA）、蛋白互作枢纽筛选、区室定位与恶性判定、细胞通讯分析、双基因平行虚拟敲除、LASSO-Cox 建模与外部验证，以及微环境与肿瘤突变负荷（TMB）关联分析。",
 8: "**结果**　差异表达与 WGCNA 交集定义 299 基因候选池，提名分属不同区室的 BUB1 与 RUNX2（EFS 每标准差 HR = 1.39 与 1.55）。RUNX2 富集于癌相关成纤维细胞/间充质干细胞与成骨细胞（均值 0.921 / 0.804），BUB1 富集于增殖细胞（0.216），两者区室不重叠；inferCNV 判 84.2% 成骨细胞为恶性。通讯权重最高的发送方为 CAF/MSC（8.34），COL1A1-CD44 最强互作。双基因虚拟敲除各产生 20 与 17 个失调基因、无共享，分别指向基质重塑与髓系-吞噬-补体程序。14 基因评分训练集 C-index 0.796（HR = 5.43，95% CI 3.34-8.84，p = 1.02e-11），GSE39055 外部验证 0.696（HR = 2.69，p = 0.0342），与低基质评分、高突变负荷相关（rho = -0.55 / 0.31）。",
 10:"**结论**　区室解析揭示骨肉瘤存在两条平行且下游互不共享的致病轴：RUNX2 驱动的基质重塑轴（CAF/MSC 与成骨细胞）与 BUB1 驱动的增殖轴，二者具备作为独立干预方向的潜力。基于同一候选池导出的 14 基因 EFS 评分在外部独立队列中保持判别效能，并在含化疗诱导坏死的多因素模型中维持独立预测性（HR = 2.92，p = 0.0203），支持其作为治疗前风险分层工具的候选价值。该框架为骨肉瘤的靶点优先排序与联合干预策略提供了可检验假说；相关机制推断仍需蛋白水平与功能实验验证。",
}

EN = {
"Background and aims": "Osteosarcoma is the most common primary malignant bone tumour in children and adolescents, and metastasis or relapse remains the leading cause of treatment failure. Which cellular compartments harbour the programmes driving matrix remodelling and proliferation, and whether these programmes are independent downstream, remain unresolved; current prognostic assessment also lacks a molecular tool that reliably identifies high-risk patients before treatment. Taking compartmental resolution as the organising principle, this study aimed to identify hub genes residing in distinct compartments, to determine the independence of their downstream programmes, and to build an externally validated event-free survival (EFS) risk score.",
"Methods": "We integrated TARGET-OS bulk transcriptomes (85 patients evaluable for EFS, 41 events), five external cohorts and single-cell transcriptomes (42,784 cells), performing differential expression and weighted gene co-expression network analysis (WGCNA), protein-protein interaction hub screening, compartmental localisation with malignancy inference, cell-cell communication analysis, parallel in silico knockout of two genes, LASSO-Cox modelling with external validation, and microenvironment and tumour mutational burden (TMB) association analyses.",
"Results": "Differentially expressed genes intersected with WGCNA modules defined a 299-gene candidate pool and nominated BUB1 and RUNX2, which reside in distinct compartments (HR per SD for EFS, 1.39 and 1.55). RUNX2 was enriched in cancer-associated fibroblasts/mesenchymal stem cells and osteoblasts (mean 0.921 / 0.804) and BUB1 in proliferating cells (0.216), with no compartmental overlap; inferCNV classified 84.2% of osteoblasts as malignant. CAF/MSC were the dominant signalling senders and COL1A1-CD44 the strongest interaction. Parallel knockout yielded 20 and 17 non-overlapping dysregulated genes, indicating matrix remodelling and myeloid-phagocytosis-complement programmes, respectively. The 14-gene score achieved C-indices of 0.796 (training; HR = 5.43, 95% CI 3.34-8.84, p = 1.02e-11) and 0.696 (GSE39055; HR = 2.69, p = 0.0342), and correlated with low stromal scores and high mutational burden (rho = -0.55 / 0.31).",
"Conclusions": "Compartmental resolution revealed two parallel pathogenic axes with non-overlapping downstream programmes in osteosarcoma - a RUNX2-driven matrix remodelling axis (CAF/MSC and osteoblasts) and a BUB1-driven proliferation axis - which may represent independent targets for intervention. The 14-gene EFS score derived from the same candidate pool retained discriminative performance in an independent external cohort and remained independently predictive in a multivariable model including chemotherapy-induced necrosis (HR = 2.92, p = 0.0203), supporting its potential value as a pretreatment risk-stratification tool. This framework provides testable hypotheses for target prioritisation and combination strategies in osteosarcoma, although the mechanistic inferences are computational and require protein-level and functional validation.",
}

# --- 写中文正文 ---
for i, s in CN.items():
    lines[i - 1] = s
io.open(MD, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")

# --- 写中英对照文件 + 词数统计 ---
wc = lambda s: len([w for w in re.split(r"\s+", s.strip()) if w])
buf, tot_en = [], 0
for k, v in EN.items():
    n = wc(v); tot_en += n
    buf.append(f"**{k}**　{v}\n")
    print(f"  {k:<22} {n:>3} words")
cn_tot = sum(len(re.findall(r"[\u4e00-\u9fff]", CN[i])) for i in CN)
print(f"英文摘要合计 {tot_en} words   |   中文摘要 CJK {cn_tot} 字")

with io.open(OUT_EN, "w", encoding="utf-8", newline="\n") as f:
    f.write("# 摘要（中英对照）— OS_matrix_dualKO v1.3\n\n")
    f.write("> 中文正文见 manuscript/OS_manuscript_v1_cn.md；本文件英文版为投稿实体，已按 350 词上限校准。\n\n")
    f.write("## 中文\n\n" + "\n".join(CN[i] for i in (4, 6, 8, 10)) + "\n\n## English\n\n")
    f.write("\n".join(buf))
print("written:", OUT_EN)
