# -*- coding: utf-8 -*-
"""P10_build_repo.py — 把 code/ 整理成可公开的代码仓库包

产物: repo/OS_matrix_dualKO_code/   (git 仓库根，可直接 push 到 GitHub/Zenodo)
规则: 复制而非移动（原 code/ 保持不动，避免破坏正在运行的流水线）
"""
import os, re, shutil, subprocess, sys, json

ROOT = r"D:/projects/OS_matrix_dualKO"
SRC  = os.path.join(ROOT, "code")
REPO = os.path.join(ROOT, "repo", "OS_matrix_dualKO_code")

# 分阶段目录 ← 脚本清单（顺序即执行顺序）
LAYOUT = {
    "01_data_download": [
        "P0_download_gdc.sh", "P0_download_geo.sh", "P0_download_v2.sh",
        "P0_build_clinical.js", "P7c_download_maf.py",
    ],
    "02_bulk_transcriptome": [
        "P1_analysis.R", "P2_hub_ML.R", "P3_prerecheck.R", "P8_pool_enrichment_test.R",
    ],
    "03_single_cell": [
        "P4_scRNA.R",
        "P4d_cellchat.R", "P4d_replot.R",
        "P4e_cellchat_percellchat.R", "P4e_merge.R", "P4e_post.R", "P4e_replot_fix.R",
        "P4f_make_gene_order.R",
        "P4g_infercnv.R", "P4g_post_infercnv.R",
        "P4g_score.R", "P4g_score2.R", "P4g_score2b.R", "P4g_score2c.R",
    ],
    "04_dual_knockout": [
        "P5_dualKO.R", "P5_final_drg.R", "P5_post.R",
    ],
    "05_prognostic_signature": [
        "P6_signature.R", "P6b_signature_OS.R", "P6c_external_expand.R",
    ],
    "06_figures_tables": [
        "P7a_roc_curves.R", "P7a2_roc_individual.R",
        "P7b_estimate.R", "P7b_fig1_design.R", "P7b_make_composites.sh",
        "P7b_fix_order.py", "P7c_tmb.R", "P7c_redraw_km.R", "P7c_split_figures.py",
        "P8g_fig_improve.R",
    ],
    "07_manuscript": [
        "P7b_md_to_docx.py", "P7b_embed_figures_tables.py",
        "P8a_rewrite_front.py", "P8b_rewrite_methods.py", "P8c_rewrite_results.py",
        "P8d_rewrite_discussion.py", "P8e_revise_v3.py", "P8h_caption_sync.py",
        "P8i_abstract_conclusion.py", "P8i2_abstract_final.py",
        "P8i3_abstract_CN_EN.py", "P8i4_abstract_EN_trim.py",
        "P9_en_trim.py", "rewrite_title.py",
    ],
    # 环境安装与一次性诊断脚本：保留以保证可追溯，但不属于分析管线
    "misc": [
        "install_pkg.cmd", "install_r_pkg.sh", "install_cellchat.cmd",
        "install_cellchat.sh", "install_rjags.cmd", "install_rjags.sh",
        "_diag_loaddb.R", "_diag_txdb.R", "_final_stats.R", "_p5_probe3.R",
        "_install_jags.cmd", "_install_make.cmd",
        "_test_msys2.cmd", "_test_rtools.R", "_test_rtools.cmd", "_test_rtools2.cmd",
    ],
}

def build():
    if os.path.isdir(REPO):
        shutil.rmtree(REPO)
    os.makedirs(REPO)
    placed, missing = [], []
    for folder, files in LAYOUT.items():
        os.makedirs(os.path.join(REPO, folder), exist_ok=True)
        for f in files:
            s = os.path.join(SRC, f)
            if not os.path.exists(s):
                missing.append(f); continue
            shutil.copy2(s, os.path.join(REPO, folder, f))
            placed.append((folder, f))

    # 未分类的脚本（防漏）
    known = {f for fs in LAYOUT.values() for f in fs}
    unclassified = sorted(f for f in os.listdir(SRC)
                          if os.path.isfile(os.path.join(SRC, f)) and f not in known)
    if unclassified:
        os.makedirs(os.path.join(REPO, "misc"), exist_ok=True)
        for f in unclassified:
            shutil.copy2(os.path.join(SRC, f), os.path.join(REPO, "misc", f))
            placed.append(("misc", f))

    print(f"已归位 {len(placed)} 个文件到 {len(LAYOUT)} 个阶段目录")
    if missing: print("源文件缺失:", missing)
    if unclassified: print("自动归入 misc/ 的未分类文件:", unclassified)
    return placed

if __name__ == "__main__":
    build()
