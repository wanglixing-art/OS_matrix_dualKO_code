#!/usr/bin/env bash
# ============================================================
# P0_download_geo_expand.sh
# 四个"扩展"队列 + 两个平台注释文件的下载脚本。
#
# 为什么单独一个脚本：GSE42352 / GSE39055 / GSE33382 / GSE87624 在项目早期
# 是手工下载的，P0_download_geo.sh 并不覆盖它们。此脚本补上这个缺口，
# 所有 URL 均已实测返回 HTTP 200（2026-10-02）。
#
# 用法:  bash P0_download_geo_expand.sh
# 依赖:  curl / wget
# 体积:  合计约 330 MB（不含可选的 2.7 GB RAW.tar）
# ============================================================
set -u
ROOT="${1:-D:/projects/OS_matrix_dualKO}"
FTP="https://ftp.ncbi.nlm.nih.gov/geo"
RAWD="$ROOT/data/raw"
GSE42352="$RAWD/GSE42352"       # DEG 队列（84 tumour vs 15 MSC/osteoblast）
EXPAND="$RAWD/GEO_expand"       # 外部验证队列（P6c 使用）

mkdir -p "$GSE42352" "$EXPAND"

dl () {  # dl <url> <dest>
  local url="$1" dest="$2"
  if [ -s "$dest" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "[skip] $(basename "$dest")  已存在（FORCE=1 可强制重下）"; return 0
  fi
  echo "[get ] $(basename "$dest")"
  curl -fL --retry 3 --retry-delay 5 -C - -o "$dest" "$url" \
    && echo "  ok $(du -h "$dest" | cut -f1)" \
    || { echo "  FAIL $url"; return 1; }
}

echo "=== 1/2  GSE42352 (differential expression cohort) ==="
# 只需要 series matrix + 平台注释；官方还提供 2.7 GB 的 GSE42352_RAW.tar（CEL 原始文件），
# 本分析不读取 CEL 层，故不下载。
dl "$FTP/series/GSE42nnn/GSE42352/matrix/GSE42352_series_matrix.txt.gz" \
   "$GSE42352/GSE42352_series_matrix.txt.gz"
dl "$FTP/platforms/GPL10nnn/GPL10295/soft/GPL10295_family.soft.gz" \
   "$GSE42352/GPL10295_family.soft.gz"

echo
echo "=== 2/2  外部验证队列 (GSE39055 / GSE33382 / GSE87624) ==="
# GSE39055: n=37，GPL14951 (Illumina HumanHT-12 WG-DASL V4.0)，EFS 端点
dl "$FTP/series/GSE39nnn/GSE39055/matrix/GSE39055_series_matrix.txt.gz" \
   "$EXPAND/GSE39055_series_matrix.txt.gz"
# GSE33382: 与 GSE21257 完全同平台 (GPL10295)，转移 <5 年
dl "$FTP/series/GSE33nnn/GSE33382/matrix/GSE33382_series_matrix.txt.gz" \
   "$EXPAND/GSE33382_series_matrix.txt.gz"
# GSE87624: 原发 vs 转移 RNA-seq。
# 注意 series matrix 只有约 2.9 KB（元数据壳，无表达值）；
# 表达数据在 suppl 的 Human_masked 文件中，两个都要下。
dl "$FTP/series/GSE87nnn/GSE87624/matrix/GSE87624_series_matrix.txt.gz" \
   "$EXPAND/GSE87624_series_matrix.txt.gz"
dl "$FTP/series/GSE87nnn/GSE87624/suppl/GSE87624_Human_masked.txt.gz" \
   "$EXPAND/GSE87624_Human_masked.txt.gz"

echo
echo "=== 探针注释映射 ==="
echo "运行以下命令生成探针->symbol 映射表（P1 与 P6c 依赖）:"
echo "  Rscript P0_build_probe_map.R \"$ROOT\""
echo
echo "完成。下一步: Rscript ../01_data_download/P0_build_probe_map.R"
