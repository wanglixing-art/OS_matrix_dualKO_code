# -*- coding: utf-8 -*-
"""
P9_en_trim.py  --  英文投稿稿合规精简（JBO: 摘要 <=300 词 / 正文 2000-5000 词）

原则（硬约束）：
  1. 只删除冗余措辞，不删除、不修改任何数值、参数、基因名、统计量
  2. 每一处替换都必须唯一匹配（count == 1），否则记为失败并整体报告
  3. 替换前后做数值集合比对：任何数值消失即报警
  4. 引用标记 [n] 数量前后必须一致

用法: python P9_en_trim.py [--dry]
"""
import re
import sys
import shutil
import datetime

BASE = "manuscript"
F_MAIN = f"{BASE}/OS_manuscript_v1_EN.md"
F_ABS = f"{BASE}/OS_abstract_CN_EN_v1.md"

# ---------------------------------------------------------------- 正文精简
BODY_EDITS = [
    # ---- Introduction, para 1 ----
    ("coinciding with the window of rapid bone growth; tumours arise predominantly",
     "coinciding with rapid bone growth; tumours arise predominantly"),
    ("and molecular tools that can reliably identify high-risk patients before treatment are lacking",
     "and molecular tools to identify high-risk patients before treatment are lacking"),
    ("Which compartments drive extracellular matrix remodelling as opposed to proliferation, and whether these programmes are transcriptionally independent downstream, have not been systematically resolved",
     "Which compartments drive extracellular matrix remodelling rather than proliferation, and whether these programmes are transcriptionally independent downstream, remain unresolved"),
    # ---- Introduction, para 2 ----
    ("scRNA-seq has allowed malignant and non-malignant cells to be distinguished at cellular resolution, programme scores compared across cell types, and prognostic genes localised to compartments such as the osteoblastic lineage and macrophages [9-13]; differential expression and co-expression network analyses have also nominated candidates",
     "scRNA-seq has distinguished malignant from non-malignant cells at cellular resolution, compared programme scores across cell types, and localised prognostic genes to compartments such as the osteoblastic lineage and macrophages [9-13]; differential expression and co-expression analyses have also nominated candidates"),
    ("First, most prognostic models are limited by single-cohort training without external independent validation, or",
     "First, most prognostic models are limited by single-cohort training without external validation, or"),
    ("Second, candidate gene screening usually stops at the level of gene lists, without defining the cellular compartments to which key genes belong",
     "Second, candidate screening usually stops at gene lists, without defining the compartments to which key genes belong"),
    ("and computationally inferred hub genes lack dissection of their downstream programmes",
     "and computationally inferred hubs lack dissection of their downstream programmes"),
    ("Most importantly, no study has yet used parallel in silico knockout at single-cell resolution to compare the downstream regulatory programmes of two hub genes located in different compartments.",
     "Most importantly, no study has used parallel in silico knockout at single-cell resolution to compare the downstream programmes of two hub genes in different compartments."),
    # ---- Introduction, para 3 ----
    ("Following the established progression from candidate screening through compartmental localisation and network inference to external validation, we introduced in silico knockout to supply the missing mechanistic step, asking whether osteosarcoma contains two hub genes occupying distinct compartments and independent downstream.",
     "Introducing in silico knockout to supply the missing mechanistic step, we asked whether osteosarcoma contains two hub genes occupying distinct compartments with independent downstream programmes."),
    ("Three principal findings emerged. First, we identified BUB1 and RUNX2 as hubs of the proliferation and stromal compartments",
     "Three principal findings emerged. First, BUB1 and RUNX2 were hubs of the proliferation and stromal compartments"),
    ("Second, parallel in silico knockout of the two genes showed that their downstream programmes differ in nature and share no dysregulated genes, supporting their orthogonality at the level of transcriptional regulation.",
     "Second, parallel in silico knockout showed that their downstream programmes differ in nature and share no dysregulated genes, supporting orthogonality at the transcriptional level."),
    ("We first describe candidate pool construction and hub gene identification, then single-cell compartmental localisation and the communication architecture, then the dual-gene in silico knockout, and finally the risk score, its external validation and its association with the microenvironment and mutational burden.",
     "We describe candidate pool construction and hub identification, single-cell localisation and communication architecture, the dual-gene knockout, and finally the risk score with its external validation and microenvironment associations."),
    # ---- Methods ----
    ("Bulk RNA-seq data for osteosarcoma were obtained from the TARGET-OS project through the NCI Genomic Data Commons (GDC) API. Eighty-eight single-sample STAR count files (GENCODE v36 annotation, non-strand-specific) were downloaded and merged into a gene",
     "Bulk RNA-seq data were obtained from the TARGET-OS project through the NCI Genomic Data Commons (GDC) API. Eighty-eight single-sample STAR count files (GENCODE v36 annotation, non-strand-specific) were merged into a gene"),
    ("After aggregation of protein-coding genes and low-expression filtering",
     "After aggregating protein-coding genes and low-expression filtering"),
    ("Public microarray and RNA-seq cohorts were used for external validation and differential expression analysis:",
     "Public cohorts were used for external validation and differential expression analysis:"),
    ("Series matrix files were parsed directly with GEOquery 2.80.0, and clinical features were extracted from characteristics_ch1 fields rather than from published descriptions.",
     "Series matrix files were parsed with GEOquery 2.80.0, and clinical features were extracted from characteristics_ch1 fields rather than published descriptions."),
    ("and the illuminaHumanv4.db package for GSE39055 (complete mapping of the 14 signature genes). Single-cell RNA-seq data were obtained from GSE162454",
     "and illuminaHumanv4.db for GSE39055 (complete mapping of the 14 signature genes). Single-cell RNA-seq data came from GSE162454"),
    ("Differential expression between tumour biopsies and normal controls was assessed with limma.",
     "Differential expression between tumour biopsies and normal controls used limma."),
    ("biweight midcorrelation (bicor) and an unsigned network type; the soft-thresholding power determined by pickSoftThreshold was 4",
     "biweight midcorrelation (bicor) and an unsigned network; the soft-thresholding power from pickSoftThreshold was 4"),
    ("identified seven modules significantly associated with outcome (p < 0.05, 2,553 genes)",
     "identified seven outcome-associated modules (p < 0.05, 2,553 genes)"),
    ("cytoHubba degree, MCC, MNC, EPC and betweenness centrality were used to quantify network centrality and were integrated with univariable Cox regression for EFS. Several machine-learning algorithms (glmBoost, Lasso, elastic net, Ridge, StepCox) were compared on the 85 training samples, with 30 features preselected at Cox p < 0.1 for comparison of event-prediction performance.",
     "cytoHubba degree, MCC, MNC, EPC and betweenness centrality quantified network centrality and were integrated with univariable Cox regression for EFS. Several machine-learning algorithms (glmBoost, Lasso, elastic net, Ridge, StepCox) were compared on the 85 training samples, with 30 features preselected at Cox p < 0.1 to compare event-prediction performance."),
    ("and the difference in restricted mean survival time (RMST) between groups, together with its 95% confidence interval, was estimated using a Greenwood variance estimator.",
     "and the between-group difference in restricted mean survival time (RMST) with its 95% confidence interval was estimated using a Greenwood variance estimator."),
    ("scaled and subjected to principal component analysis (40 principal components), followed by cross-sample integration with Harmony 2.0.5.",
     "scaled and subjected to principal component analysis (40 components), then integrated across samples with Harmony 2.0.5."),
    ("yielding 38 clusters, which were annotated to nine cell types by marker assignment and by SingleR 2.14.2",
     "yielding 38 clusters annotated to nine cell types by marker assignment and SingleR 2.14.2"),
    ("CellChat 2.2.0.9001 was applied to the annotated dataset at cell type level using CellChatDB.human, with 100 permutations",
     "CellChat 2.2.0.9001 was applied at cell type level using CellChatDB.human, with 100 permutations"),
    ("and correlated with the risk score. TARGET-OS somatic mutation data were obtained from the GDC as 167 masked somatic mutation MAF files (142 patients).",
     "and correlated with the risk score. Somatic mutation data were obtained from the GDC as 167 masked MAF files (142 patients)."),
    ("splice site, translation start site and nonstop mutation), deduplicated within a patient",
     "splice site, translation start site and nonstop), deduplicated within a patient"),
    ("survival differences were assessed by log-rank test, and enrichment was assessed by hypergeometric distribution. Two-sided p < 0.05 was considered statistically significant.",
     "survival differences by log-rank test, and enrichment by hypergeometric distribution. Two-sided p < 0.05 was significant."),
    ("and all intermediate result tables were retained to ensure reproducibility.",
     "and all intermediate result tables were retained for reproducibility."),
    # ---- Results ----
    ("produced 21 modules, seven significantly associated with event-free survival (EFS) or overall survival (OS) event status.",
     "produced 21 modules, seven associated with event-free survival (EFS) or overall survival (OS) event status."),
    ("whose main cluster consisted of cell-cycle genes, with a separate collagen/extracellular matrix subcluster and a third immune cluster.",
     "whose main cluster comprised cell-cycle genes, with separate collagen/extracellular matrix and immune subclusters."),
    ("The first was BUB1, a spindle assembly checkpoint kinase, which had the joint highest degree (22) and was significantly associated with EFS",
     "The first was BUB1, a spindle assembly checkpoint kinase with the joint highest degree (22), significantly associated with EFS"),
    ("and both were significantly overexpressed in tumours relative to normal controls",
     "and both were overexpressed in tumours relative to normal controls"),
    ("In GSE21257, both genes were directionally consistent with poorer prognosis (HR 1.27 and 1.14) and both were significantly associated with metastatic status",
     "In GSE21257, both were directionally consistent with poorer prognosis (HR 1.27 and 1.14) and significantly associated with metastatic status"),
    ("After quality control, 42,784 cells from the six GSE162454 samples clustered into 38 clusters and were annotated as nine cell types",
     "After quality control, 42,784 cells from the six GSE162454 samples formed 38 clusters annotated as nine cell types"),
    ("showed only 48.4% concordance, indicating that approach missed 2,982 cells with CNV evidence of malignancy",
     "showed only 48.4% concordance, that approach missing 2,982 cells with CNV evidence of malignancy"),
    ("CAF/MSC were the dominant signal senders, with an outgoing interaction weight of 8.34",
     "CAF/MSC were the dominant senders, with an outgoing interaction weight of 8.34"),
    ("The COLLAGEN pathway had the highest mean weight (11.31), 2.3-fold that of the second-ranked pathway (CV 0.34).",
     "The COLLAGEN pathway had the highest mean weight (11.31), 2.3-fold the second-ranked pathway (CV 0.34)."),
    ("the AP-1 transcription complex (FOS, FOSB, JUN, JUNB), the small leucine-rich proteoglycan (SLRP) matrix core (DCN Z = 2.76, BGN Z = 2.14, LUM Z = 2.10, PCOLCE Z = 2.06, together with COL3A1 Z = 2.35)",
     "the AP-1 complex (FOS, FOSB, JUN, JUNB), the small leucine-rich proteoglycan (SLRP) matrix core (DCN Z = 2.76, BGN Z = 2.14, LUM Z = 2.10, PCOLCE Z = 2.06, with COL3A1 Z = 2.35)"),
    ("but with a different downstream profile concentrated on myeloid, phagocytic and complement programmes",
     "but with a different profile concentrated on myeloid, phagocytic and complement programmes"),
    ("RUNX2 was unperturbed in the BUB1 knockout network (Z = 0.40) and BUB1 unperturbed in the RUNX2 knockout network (Z = -0.75).",
     "RUNX2 was unperturbed in the BUB1 knockout network (Z = 0.40) and BUB1 in the RUNX2 network (Z = -0.75)."),
    ("BUB1 knockout DRGs were significantly enriched in the candidate pool",
     "BUB1 knockout DRGs were enriched in the candidate pool"),
    ("The coefficients, directions and functional annotations of the 14 genes are given in Table 2. Ten genes contributed risk and four were protective (individual coefficients in Table 2), placing collagen genes",
     "Coefficients, directions and functional annotations are given in Table 2: ten genes contributed risk and four were protective, placing collagen genes"),
    ("(beyond which fewer than five patients remained at risk per group and estimates became unstable), the analysis was also evaluated",
     "(beyond which fewer than five patients remained at risk per group), the analysis was also evaluated"),
    ("RMST does not rely on the proportional hazards assumption and is the recommended summary measure when curves cross; the full follow-up curves are shown in Supplementary Fig. S11.",
     "RMST does not rely on proportional hazards and is the recommended summary measure when curves cross; full follow-up curves are in Supplementary Fig. S11."),
    ("the median-split log-rank test was not significant (p = 0.126) while the per-SD Cox regression was (HR = 1.39, p = 0.039)",
     "the median-split log-rank test was not significant (p = 0.126) while per-SD Cox regression was (HR = 1.39, p = 0.039)"),
    ("Sensitivity analysis showed that retraining on the OS endpoint in GSE21257 produced an 18-gene model (including RUNX2) whose external performance (C-index 0.562, p = 0.186) was no better than the frozen-coefficient model, suggesting the weak performance arises mainly from population, platform and endpoint differences rather than overfitting or endpoint mismatch.",
     "Retraining on the OS endpoint in GSE21257 produced an 18-gene model (including RUNX2) whose external performance (C-index 0.562, p = 0.186) was no better than the frozen-coefficient model, suggesting the weak performance arises mainly from population, platform and endpoint differences rather than overfitting."),
    ("; BUB1 was also a prioritised candidate target in genome-scale CRISPR-Cas9 dependency screens [29].",
     "; BUB1 was also prioritised in genome-scale CRISPR-Cas9 dependency screens [29]."),
    # ---- Discussion ----
    ("provides a cellular-origin explanation for this evidence; BUB1 represents",
     "provides a cellular-origin explanation; BUB1 represents"),
    ("At single-cell resolution, RUNX2 was mainly localised to CAF/MSC and osteoblastic cells and BUB1 to proliferating and malignant cells, with almost no overlap between the two compartments. inferCNV malignancy calling agreed with marker-based proxy classification in only 48.4% of cases",
     "At single-cell resolution, RUNX2 localised mainly to CAF/MSC and osteoblastic cells and BUB1 to proliferating and malignant cells, with almost no overlap. inferCNV malignancy calling agreed with marker-based classification in only 48.4% of cases"),
    ("the COLLAGEN pathway 2.3-fold heavier than the second-ranked pathway, and CAF/MSC ranked first",
     "the COLLAGEN pathway 2.3-fold heavier than the next, and CAF/MSC ranked first"),
    ("BUB1 knockout did not significantly perturb cell-cycle genes (highest NUSAP1 Z = 1.75), and results should be interpreted as a regulatory neighbourhood rather than evidence of functional necessity.",
     "BUB1 knockout did not strongly perturb cell-cycle genes (highest NUSAP1 Z = 1.75), and results represent a regulatory neighbourhood rather than functional necessity."),
    ("but cross-platform and cross-endpoint error remain the main sources of external performance variability;",
     "but cross-platform and cross-endpoint error remain the main sources of external variability;"),
    ("consistent with reports of immune exclusion accompanied by high mutational burden and suppressive compartment enrichment [41-44].",
     "consistent with reports of immune exclusion with high mutational burden and suppressive compartment enrichment [41-44]."),
    ("Second, this is a purely computational analysis: the mechanistic inferences rely on external anchoring to published pharmacological and chromatin evidence and lack wet-laboratory validation here, so causal roles are not established.",
     "Second, this is a purely computational analysis: the mechanistic inferences rely on anchoring to published pharmacological and chromatin evidence and lack wet-laboratory validation, so causal roles are not established."),
    ("and the virtual knockout significance threshold (Z > 2 with distance",
     "and the virtual knockout threshold (Z > 2 with distance"),
    # ---- Conclusions ----
    ("this study built a progressive framework spanning candidate pool screening, hub gene identification, single-cell localisation",
     "this study built a progressive framework spanning candidate screening, hub identification, single-cell localisation"),
    ("although the mechanistic inferences require further protein-level and functional validation.",
     "although the mechanistic inferences require protein-level and functional validation."),
]

# ------------------------------------------------- 摘要精简（主稿 + 中英对照稿）
ABS_EDITS = [
    ("adolescents, and metastasis or relapse remains the leading cause",
     "adolescents; metastasis or relapse remains the leading cause"),
    ("before treatment all remain unresolved",
     "before treatment remain unresolved"),
    ("whether these programmes are downstream-independent",
     "whether these are downstream-independent"),
    ("parallel in silico knockout of two genes, LASSO-Cox modelling",
     "parallel in silico knockout, LASSO-Cox modelling"),
    ("and microenvironment and tumour mutational burden (TMB) analyses.",
     "and microenvironment and mutational burden analyses."),
    ("and BUB1 in proliferating cells (0.216), with no overlap",
     ", BUB1 in proliferating cells (0.216), with no overlap"),
    ("for target prioritisation and combination strategies, pending functional validation.",
     "for target prioritisation and combination, pending functional validation."),
]


# ------------------------------------------------------------------ utilities
def wc(s):
    """Word 口径：去 markdown 修饰后按空白分词"""
    s = re.sub(r"[*_`]", "", s)
    return len(re.findall(r"\S+", s))


def nums(s):
    return re.findall(r"\d+(?:\.\d+)?(?:e-?\d+)?", s)


def cites(s):
    return re.findall(r"\[[\d,\s\-]+\]", s)


def apply_edits(text, edits, label, failures):
    for i, (old, new) in enumerate(edits, 1):
        n = text.count(old)
        if n != 1:
            failures.append((label, i, n, old[:70]))
            continue
        text = text.replace(old, new)
    return text


def body_span(lines):
    a = next(i for i, l in enumerate(lines) if re.match(r"^#+\s*1\.\s*Introduction", l))
    b = next(i for i, l in enumerate(lines) if re.match(r"^#+\s*CRediT", l))
    return a, b


def main():
    dry = "--dry" in sys.argv
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    failures = []

    # ---- 主稿 ----
    src = open(F_MAIN, encoding="utf-8").read()
    lines = src.split("\n")
    a, b = body_span(lines)
    body_before = "\n".join(lines[a:b])
    nums_before = set(nums(body_before))
    cites_before = len(cites(body_before))

    new = apply_edits(src, BODY_EDITS, "body", failures)
    new = apply_edits(new, ABS_EDITS, "abstract", failures)

    nlines = new.split("\n")
    a2, b2 = body_span(nlines)
    body_after = "\n".join(nlines[a2:b2])
    nums_after = set(nums(body_after))
    cites_after = len(cites(body_after))

    lost_nums = sorted(nums_before - nums_after)
    gained_nums = sorted(nums_after - nums_before)

    print("=" * 72)
    print(f"正文（含小标题，Word 口径）: {wc(body_before)} -> {wc(body_after)}")
    print(f"正文（纯叙述）          : {wc(chr(10).join(l for l in nlines[a2:b2] if not l.lstrip().startswith('#')))}"
          f"  (原 {wc(chr(10).join(l for l in lines[a:b] if not l.lstrip().startswith('#')))})")
    print(f"引用标记数               : {cites_before} -> {cites_after}")
    print(f"数值丢失                 : {lost_nums if lost_nums else '无'}")
    print(f"数值新增                 : {gained_nums if gained_nums else '无'}")

    # 摘要
    def abs_words(text, tag):
        ls = text.split("\n")
        i0 = next(i for i, l in enumerate(ls) if re.match(r"^#+\s*Abstract", l))
        i1 = next(i for i, l in enumerate(ls) if re.match(r"^#+\s*Highlights", l))
        seg = "\n".join(ls[i0 + 1:i1])
        seg = re.sub(r"^\*\*Keywords:.*$", "", seg, flags=re.M)
        return wc(seg)

    print(f"摘要（含结构标签）        : {abs_words(src,'o')} -> {abs_words(new,'n')}")
    # 逐段
    print("\n逐段英文摘要词数（含标签）:")
    ls = new.split("\n")
    i0 = next(i for i, l in enumerate(ls) if re.match(r"^#+\s*Abstract", l))
    i1 = next(i for i, l in enumerate(ls) if re.match(r"^#+\s*Highlights", l))
    tot = 0
    for l in ls[i0 + 1:i1]:
        if l.strip() and not l.strip().startswith("**Keywords"):
            lab = re.match(r"\*\*(.+?)\*\*", l)
            print(f"   {lab.group(1) if lab else '?':<20} {wc(l)}")
            tot += wc(l)
    print(f"   {'TOTAL (incl. labels)':<20} {tot}")

    if failures:
        print("\n!!! 未唯一匹配的替换（未应用）:")
        for lab, i, n, s in failures:
            print(f"   [{lab} #{i}] count={n} :: {s}")
        print("\n=> 有失败项，未写盘。请修正后重跑。")
        return 1

    if dry:
        print("\n[dry-run] 未写盘。")
        return 0

    shutil.copy(F_MAIN, f"{BASE}/_pre_p9_{stamp}_OS_manuscript_v1_EN.md")

    # 摘要文件（中英对照，仅改英文段）
    absrc = open(F_ABS, encoding="utf-8").read()
    abnew = apply_edits(absrc, ABS_EDITS, "abstract_file", failures)
    if failures:
        print("摘要文件替换失败:", failures)
        return 1
    shutil.copy(F_ABS, f"{BASE}/_pre_p9_{stamp}_OS_abstract_CN_EN_v1.md")

    open(F_MAIN, "w", encoding="utf-8", newline="\n").write(new)
    open(F_ABS, "w", encoding="utf-8", newline="\n").write(abnew)
    print(f"\n已写盘。快照: _pre_p9_{stamp}_*.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
