#!/usr/bin/env bash
# P7b_make_composites.sh — 生成主图 F1-F10（PDF 矢量拼图）
set -u
PY="C:/Users/Administrator/.workbuddy/binaries/python/versions/3.13.12/python.exe"
MERGE="C:/Users/Administrator/.workbuddy/skills/figure-merge/scripts/figure_merge.py"
cd "D:/projects/OS_matrix_dualKO/results/figures" || exit 1

run () {  # run <out> <comma-list>
  echo "=== $1 ==="
  "$PY" "$MERGE" --files "$2" -o "$1" 2>&1 | tail -4
}

# F1 研究设计总览 —— 由 R 另行绘制（流程图），此处占位由 F1_design.pdf 提供
if [ -f F1_design.pdf ]; then cp F1_design.pdf F1.pdf; echo "=== F1 <- F1_design.pdf ==="; fi

# F3 单细胞图谱与双枢纽定位
run F3.pdf "P4_QC_violin_before.png,P4_QC_violin_after.png,P4_UMAP_clusters.png,P4_UMAP_celltypes.png,P4_UMAP_samples.png,P4_marker_dotplot.png,P4_hub_featureplot.png,P4_hub_violin_celltype.png"

# F4 inferCNV 恶性判定 + CellChat 通讯
run F4.pdf "P4g_infercnv_heatmap_obs.png,P4g_infercnv_malignancy_summary.png,P4g_hub_vs_infercnv_malignant.png,P4_cellchat_circle_weight.png,P4_cellchat_heatmap_weight.png,P4e_hub_outgoing_by_sample.png,P4e_stroma_vs_proliferation_ratio.png"

# F5 双 KO 扰动景观
run F5.pdf "P5_dualKO_Z_landscape.png,P5_DRG_venn.png,P5_axis_markers_perturbation.png"

# F6 双 KO 富集对比
run F6.pdf "P5_GO_RUNX2_specific.png,P5_GO_BUB1_specific.png,P5_GO_RUNX2_all.png,P5_GO_BUB1_all.png"

# F7 签名系数 + 训练 KM + 训练 timeROC
run F7.pdf "P6_signature_coefs.png,P6_KM_train.png,P7_timeROC_TARGET_train.png"

# F8 外部验证 GSE39055
run F8.pdf "P6c_KM_GSE39055.png,P7_timeROC_GSE39055.png,P6c_forest_Cindex.png"

# F9 转移维度三队列
run F9.pdf "P6c_met_boxplots.png,P6_KM_GSE21257.png"

# F10 ESTIMATE + TMB
run F10.pdf "P7_estimate_boxplots.png,P7_tmb_scatter.png"

echo
echo "=== 结果 ==="
ls -la F*.pdf 2>/dev/null
