# Stage 01 — Data download

All data are public. **No raw data are redistributed in this repository** (raw + intermediate
footprint is ~12 GB). Re-download from the sources below.

## Expected directory tree

Scripts assume this layout relative to the project root (edit the `ROOT` constant at the
top of each script):

```
data/
├── raw/
│   ├── TARGET-OS/
│   │   ├── gdc_cases_expanded.json        # clinical entity export (GDC API)
│   │   └── ...                            # STAR count files (GENCODE v36)
│   ├── TARGET-OS_maf/                     # masked somatic mutation MAF slice
│   ├── GSE42352/                          # DEG cohort (84 tumour vs 15 normal)
│   │   ├── GSE42352_series_matrix.txt.gz
│   │   ├── GPL10295_family.soft.gz        # platform annotation (320 MB)
│   │   └── GPL10295_probe_symbol.tsv      # <- built by P0_build_probe_map.R
│   ├── GSE21257/                          # external validation (OS endpoint)
│   ├── GSE162454/                         # single-cell (10x)
│   └── GEO_expand/                        # external validation cohorts for P6c
│       ├── GSE39055_series_matrix.txt.gz
│       ├── GSE33382_series_matrix.txt.gz
│       ├── GSE87624_series_matrix.txt.gz
│       ├── GSE87624_Human_masked.txt.gz
│       └── GPL14951_probe_symbol.tsv      # <- built by P0_build_probe_map.R
├── processed/                             # intermediate matrices written by the pipeline
└── (results/ at project root)             # figures, tables, RDS objects
```

Create the empty tree before the first run.

## Scripts and order

| Order | Script | What it fetches | Notes |
|---|---|---|---|
| 1 | `P0_download_gdc.sh` | TARGET-OS STAR counts + clinical entity via the GDC API | 88 single-sample count files (GENCODE v36, non-strand-specific) |
| 2 | `P0_download_geo.sh` | GSE21257, GSE162454 (+ GSE152048, GSE16088 used for the inferCNV reference) | via GEO FTP |
| 3 | `P0_download_geo_expand.sh` | **GSE42352, GSE39055, GSE33382, GSE87624 + GPL10295 platform** | ⚠️ see note below |
| 4 | `P0_download_v2.sh` | Supplementary / retry pass for cohorts that failed on the first run | idempotent |
| 5 | `P0_build_probe_map.R` | *(no download)* builds the probe → symbol maps required by `P1_analysis.R` and `P6c_external_expand.R` | needs `GPL10295_family.soft.gz`; GPL14951 route uses `illuminaHumanv4.db` |
| 6 | `P0_build_clinical.js` | *(no download)* builds `data/processed/TARGET-OS_clinical.csv` from the GDC cases JSON | Node.js: `node P0_build_clinical.js` |
| 7 | `P7c_download_maf.py` | GDC masked somatic mutation MAF files (167 files / 142 cases, WXS, open tier) | used later by `06_figures_tables/P7c_tmb.R` |

> ⚠️ **Reproducibility note (honest disclosure).** The four cohorts fetched by
> `P0_download_geo_expand.sh` were originally obtained **manually** during the study, and no
> download script for them existed in the working directory. That script was written
> afterwards for this repository, and **every URL in it was verified to return HTTP 200** on
> 2026-10-02. The file layout it produces matches the paths the analysis scripts read. If
> you are re-running the pipeline from scratch, run it *before* the platform-map step.
>
> Two practical points: the official `GSE42352_RAW.tar` (2.7 GB of CEL files) is **not**
> needed — the analysis reads the series matrix and the platform annotation only. And the
> `GSE87624_series_matrix.txt.gz` file is only ~2.9 KB (a metadata shell with no expression
> values), so the expression data must come from `GSE87624_Human_masked.txt.gz`; both are
> downloaded by the script.

## Two data-handling decisions worth knowing

**1. Clinical variables come from the files, not from the papers.**
Covariates (EFS/OS time and event, metastatic status, necrosis rate, chemotherapy regimen)
were parsed from the GEO `characteristics_ch1` fields of the series matrices and from the
GDC clinical entity export — **not** transcribed from published tables. Published
descriptions and matrix metadata disagree in places; the manuscript's numbers follow the
files. If you re-analyse these cohorts, do the same.

**2. Probe annotation is platform-specific.**
GPL10295 (used by GSE21257 and GSE33382, which are the same platform) is mapped via
`GPL10295_probe_symbol.tsv`. GSE39055 is on GPL14951 and is mapped via
`illuminaHumanv4.db`, whose `PROBEID` keys are already in the `ILMN_` format used by the
series matrix. That route maps **all 14 signature genes** completely, which is why it is
preferred over downloading the 2.4 GB GPL14951 platform file.
