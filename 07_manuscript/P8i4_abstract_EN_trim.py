# P8i4: 英文摘要压缩至 <=350 词（保留全部关键数值），重写中英对照文件的 English 段
import io, os, re

BASE = r"D:\projects\OS_matrix_dualKO"
OUT_EN = os.path.join(BASE, "manuscript", "OS_abstract_CN_EN_v1.md")

EN = {
"Background and aims": "Osteosarcoma is the most common primary malignant bone tumour in children and adolescents, and metastasis or relapse remains the leading cause of treatment failure. Which compartments harbour the matrix-remodelling and proliferation programmes, and whether these programmes are downstream-independent, remain unresolved; no molecular tool reliably identifies high-risk patients before treatment. This study therefore aimed to identify hub genes in distinct compartments, assess the downstream independence of their programmes, and build an externally validated event-free survival (EFS) risk score.",
"Methods": "We integrated TARGET-OS bulk transcriptomes (85 patients evaluable for EFS, 41 events), five external cohorts and single-cell transcriptomes (42,784 cells), performing differential expression and weighted gene co-expression network analysis (WGCNA), protein-protein interaction hub screening, compartmental localisation with malignancy inference, cell-cell communication analysis, parallel in silico knockout of two genes, LASSO-Cox modelling with external validation, and microenvironment and tumour mutational burden (TMB) association analyses.",
"Results": "Intersecting differentially expressed genes with WGCNA modules defined a 299-gene candidate pool and nominated BUB1 and RUNX2, which reside in distinct compartments (HR per SD for EFS: 1.39 and 1.55). RUNX2 was enriched in cancer-associated fibroblasts/mesenchymal stem cells and osteoblasts (0.921 / 0.804) and BUB1 in proliferating cells (0.216), with no overlap; inferCNV classified 84.2% of osteoblasts as malignant. CAF/MSC were the dominant signal senders and COL1A1-CD44 the strongest interaction. Parallel knockout yielded 20 and 17 non-overlapping dysregulated genes, indicating matrix remodelling and myeloid-phagocytosis-complement programmes. The 14-gene score achieved C-indices of 0.796 in training (HR = 5.43, 95% CI 3.34-8.84, p = 1.02e-11) and 0.696 in GSE39055 (HR = 2.69, p = 0.0342), correlating with low stromal scores and high mutational burden (rho = -0.55 / 0.31).",
"Conclusions": "Compartmental resolution revealed two parallel pathogenic axes with non-overlapping downstream programmes in osteosarcoma: a RUNX2-driven matrix remodelling axis (CAF/MSC, osteoblasts) and a BUB1-driven proliferation axis, potentially representing independent targets for intervention. The 14-gene EFS score derived from the same candidate pool retained discriminative performance externally and remained independently predictive after adjustment for chemotherapy-induced necrosis (HR = 2.92, p = 0.0203), supporting its potential as a pretreatment risk-stratification tool. This framework provides testable hypotheses for target prioritisation and combination strategies, although the mechanistic inferences require protein-level and functional validation.",
}

wc = lambda s: len([w for w in re.split(r"\s+", s.strip()) if w])
tot = 0
for k, v in EN.items():
    n = wc(v); tot += n
    print(f"  {k:<22} {n:>3} words")
print(f"英文摘要合计 {tot} words  (目标 <=350, 余量 {350-tot})")

txt = io.open(OUT_EN, encoding="utf-8").read()
head = txt.split("## English")[0]
body = "\n".join(f"**{k}**　{v}\n" for k, v in EN.items())
io.open(OUT_EN, "w", encoding="utf-8", newline="\n").write(head + "## English\n\n" + body)
print("updated:", OUT_EN)
