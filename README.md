# Compartmental resolution of matrix-remodelling and proliferation programmes in osteosarcoma

Analysis code accompanying the manuscript:

> **Compartmental resolution of matrix-remodelling and proliferation programmes in osteosarcoma:
> parallel dual-gene in silico knockout and a 14-gene event-free survival score**
> Sun R, Chen S, Cheng J, Gong Y, Liu Y, Xiong H, Pi W, Yin C, Wang Q, Cheng J, Wang L\* (*corresponding author)
> Submitted to *Journal of Bone Oncology* (Elsevier)

This repository contains the full analysis pipeline, from raw public data download to the
final figures and tables. It is a flat, script-based pipeline organised into eight stages.

---

## 1. What the study does

The study asks whether osteosarcoma contains **two hub genes occupying distinct cellular
compartments with transcriptionally independent downstream programmes**, and whether a
signature derived from the shared candidate pool can predict event-free survival (EFS).

Three linked claims are made:

| # | Claim | Evidence in the code | Manuscript |
|---|---|---|---|
| 1 | Two hubs occupy non-overlapping compartments | `P4a-c` single-cell localisation; `P4g` inferCNV malignancy calling | Figs. 3, 4; Suppl. Figs. S3, S4 |
| 2 | Their downstream programmes are orthogonal | `P5` parallel dual-gene in silico knockout (shared DRGs = 0) | Figs. 5, 6; Suppl. Fig. S5 |
| 3 | A 14-gene EFS score derived from the same pool validates externally | `P6` LASSO-Cox with frozen coefficients; `P6c` external expansion | Figs. 7, 8, 10; Tables 2–4 |

The two hub genes are **RUNX2** (matrix-remodelling axis; CAF/MSC + osteoblastic
compartments) and **BUB1** (proliferation axis).

---

## 2. Data sources (all public — no new data were generated)

| Dataset | Role in the study | Source |
|---|---|---|
| **TARGET-OS** bulk RNA-seq (STAR counts, GENCODE v36) + masked somatic mutation MAF + clinical entity | Training set for the signature; TMB; ESTIMATE | NCI Genomic Data Commons — <https://portal.gdc.cancer.gov/projects/TARGET-OS> |
| **GSE42352** | Differential expression (84 tumour biopsies vs 15 normal MSC/osteoblast samples) | GEO |
| **GSE21257** | External validation (n = 53, pre-chemotherapy biopsies; OS endpoint) | GEO |
| **GSE39055** | Primary external validation (n = 37, evaluable EFS n = 36, 18 events) | GEO |
| **GSE33382** | Metastasis-within-5-years comparison (n = 53 after filtering) | GEO |
| **GSE87624** | Primary vs metastatic RNA-seq (n = 52) | GEO |
| **GSE162454** | Single-cell RNA-seq (6 samples, 10x Genomics; 42,784 cells after QC) | GEO |

Download scripts: `01_data_download/`. Clinical variables were parsed from the GEO
`characteristics_ch1` fields rather than from published descriptions, and from the GDC
clinical entity export.

---

## 3. Repository layout

```
01_data_download/         Fetch raw data (GDC API, GEO, GDC MAF slice)
02_bulk_transcriptome/    P1 limma DEG + WGCNA -> 299-gene pool; P2 STRING/cytoHubba hubs
03_single_cell/           P4a-c QC/integration/annotation; P4d CellChat; P4e-e per-sample
                          CellChat; P4g inferCNV + hub scoring (run in this order)
04_dual_knockout/         P5 scTenifoldKnk parallel dual-gene knockout (core innovation)
05_prognostic_signature/  P6 LASSO-Cox + training-cohort modelling; P6c external cohorts
06_figures_tables/        P7a time-ROC; P7b ESTIMATE; P7c TMB/KM redraws; figure assembly
07_manuscript/            Markdown -> DOCX conversion and manuscript text utilities
environment/              R session info (real package versions) and setup notes
misc/                     Environment installers and one-off diagnostic scripts
```

Each stage folder is numbered in execution order. Scripts are prefixed `P<n>` matching
the analysis phase described in the manuscript Methods.

---

## 4. Software environment

Analysis was performed under **R 4.6.1 (2026-06-24 ucrt, x86_64-w64-mingw32)** on Windows 11.

Exact package versions are recorded in **`environment/R_sessionInfo.txt`** (generated from a
live `installed.packages()` call, not hand-written). Key versions:

| Package | Version | Package | Version |
|---|---|---|---|
| edgeR | 4.10.5 | Seurat | 5.5.1 |
| limma | 3.68.5 | harmony | 2.0.5 |
| WGCNA | 1.74 | SingleR | 2.14.2 |
| clusterProfiler | 4.20.0 | infercnv | 1.28.0 |
| GEOquery | 2.80.0 | CellChat | 2.2.0.9001 |
| glmnet | 5.0 | scTenifoldKnk | 1.1 |
| timeROC | 0.4.1 | scTenifoldNet | 1.4 |
| survminer | 0.5.2 | estimate | 1.0.13 |
| org.Hs.eg.db | 3.23.1 | illuminaHumanv4.db | 1.26.0 |

Additional runtimes used by individual scripts:

- **Node.js 22** — only for `P0_build_clinical.js` (GDC cases JSON → clinical table)
- **Python 3.13** — only for the DOCX/figure assembly scripts in `06_` and `07_`
  (`python-docx`, `pandas`, `matplotlib`)

`environment/install_packages.R` installs the R dependency set. Note that `scTenifoldKnk`
and `CellChat` are GitHub/Bioconductor packages and require a C++ toolchain
(Rtools 45 on Windows; see `misc/install_*.sh`).

---

## 5. How to run

1. **Set the project root.** Every script begins with an absolute `ROOT` (or `root`)
   constant. Change it to your own checkout path before running:

   ```r
   ROOT <- "D:/projects/OS_matrix_dualKO"   # <- edit this
   ```

2. **Create the data tree** expected by the scripts:

   ```
   data/raw/       # downloaded source files (GDC + GEO)
   data/processed/ # intermediate matrices
   results/        # figures, tables, intermediate RDS objects
   ```

3. **Fetch the data** — run these in order:

   ```bash
   bash 01_data_download/P0_download_gdc.sh          # TARGET-OS (GDC API)
   bash 01_data_download/P0_download_geo.sh          # GSE21257, GSE162454, + inferCNV reference
   bash 01_data_download/P0_download_geo_expand.sh   # GSE42352, GSE39055, GSE33382, GSE87624
   Rscript 01_data_download/P0_build_probe_map.R     # probe -> symbol maps (P1 / P6c need these)
   node 01_data_download/P0_build_clinical.js        # TARGET-OS clinical table
   ```

   Large files (raw RNA-seq, GEO series matrices) are **not** redistributed here; everything
   is re-fetched from the sources in Section 2. **See `01_data_download/README.md` for an
   important note about how four of the cohorts were originally obtained.**

4. **Run the stages in numerical order.** Long-running steps: `P4_scRNA.R`,
   `P4g_infercnv.R`, `P5_dualKO.R` (scTenifoldKnk network construction is the slowest step).

5. **Build the figures** with the `06_figures_tables/` scripts, then the manuscript with
   `07_manuscript/P7b_md_to_docx.py` (usage: `script.py <in.md> <out.docx> [full|main|supp]`).

---

## 6. Script → output map

| Manuscript item | Produced by |
|---|---|
| Fig. 1 (study design) | `06_figures_tables/P7b_fig1_design.R` |
| Fig. 2A (volcano) | `02_bulk_transcriptome/P1_analysis.R` |
| Fig. 2B (GO dot plot, 299-gene pool) | `02_bulk_transcriptome/P8_pool_enrichment_test.R` |
| Fig. 2C (STRING core subnetwork) | `02_bulk_transcriptome/P2_hub_ML.R` |
| Fig. 2D (RUNX2 KM, EFS) | `02_bulk_transcriptome/P3_prerecheck.R` |
| Fig. 3 (single-cell atlas, hub localisation) | `03_single_cell/P4_scRNA.R` |
| Fig. 4A-B (inferCNV malignancy, hub comparisons) | `03_single_cell/P4g_infercnv.R` + `P4g_post_infercnv.R` |
| Fig. 4C-D (CellChat network, pathway heatmap) | `03_single_cell/P4d_cellchat.R` |
| Fig. 4E-F (per-sample signalling) | `03_single_cell/P4e_cellchat_percellchat.R` + `P4e_post.R` |
| Fig. 5 (dual-knockout perturbation landscape) | `04_dual_knockout/P5_dualKO.R` + `P5_post.R` |
| Fig. 6 (GO comparison of the two knockout programmes) | `04_dual_knockout/P5_final_drg.R` |
| Fig. 7 (training-cohort performance) | `05_prognostic_signature/P6_signature.R` |
| Fig. 8 (GSE39055 external validation) | `05_prognostic_signature/P6c_external_expand.R` + `06_figures_tables/P7a_roc_curves.R` |
| Fig. 9 (metastatic phenotype, GSE33382/GSE87624) | `05_prognostic_signature/P6c_external_expand.R` |
| Fig. 10 (ESTIMATE scores vs risk score) | `06_figures_tables/P7b_estimate.R` |
| Table 1 (cohorts and clinical characteristics) | `01_data_download/P0_build_clinical.js` + `02_bulk_transcriptome/P1_analysis.R` |
| Table 2 (14-gene signature coefficients) | `05_prognostic_signature/P6_signature.R` |
| Table 3 (multivariable independence) | `05_prognostic_signature/P6c_external_expand.R` |
| Table 4 (performance across cohorts) | `06_figures_tables/P7a_roc_curves.R` |
| Suppl. Figs. S1–S2, S10 | `02_bulk_transcriptome/P1_analysis.R`, `P2_hub_ML.R`, `P3_prerecheck.R` |
| Suppl. Figs. S3–S5 | `03_single_cell/*`, `04_dual_knockout/*` |
| Suppl. Figs. S6–S7, S11 | `05_prognostic_signature/P6c_external_expand.R`, `06_figures_tables/P7c_tmb.R`, `P7c_redraw_km.R` |
| Figure/table/SI assembly into the manuscript | `06_figures_tables/P7b_make_composites.sh`, `P7c_split_figures.py`, `07_manuscript/P7b_embed_figures_tables.py` |

---

## 7. Methodological notes and caveats

- **Frozen signature coefficients.** The 14-gene LASSO-Cox model is fitted once on
  TARGET-OS and then applied unchanged to every external cohort (`P6c`). No re-fitting on
  validation data. This is deliberate: it makes the training-to-external performance
  gradient (C-index 0.796 → 0.696 → 0.547) interpretable rather than optimistic.
- **inferCNV reference.** Malignancy calling uses normal/reference compartments as the
  baseline with a 95th-percentile threshold (`thr = 0.050`). Osteoblastic cells were called
  malignant in 84.2% of cases; the reference compartments served as the negative control.
- **Virtual knockout is computational, not experimental.** `scTenifoldKnk` (`P5`) models
  the effect of removing one gene from a single-cell gene regulatory network. The reported
  downstream programmes are **hypotheses requiring functional validation**, and the
  manuscript states this explicitly.
- **Endpoint heterogeneity across validation cohorts.** GSE21257 and several other cohorts
  report overall survival rather than EFS; the figures label endpoints per cohort rather
  than pooling them.
- **Follow-up truncation.** In GSE39055 the KM curve is truncated at 5 years because the
  number at risk fell to ≤ 4 thereafter (Suppl. Fig. S11 shows the untruncated curve).
- **No new data.** Every dataset is public; the study is entirely computational.

---

## 8. Citation and licence

See `CITATION.cff`. A citable archived snapshot of this code is intended to be deposited
alongside the paper; if you use these scripts, please cite the manuscript above.

Licence: MIT (see `LICENSE`).

## 9. Contact

Corresponding author: **Lixing Wang**, Department of Gastrointestinal Surgery (Ward II),
The Second Affiliated Hospital of Kunming Medical University, Kunming 650101, Yunnan, China.
E-mail: <2300240803@qq.com>
