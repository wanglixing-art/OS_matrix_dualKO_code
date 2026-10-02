#!/bin/bash
# P0: Download TARGET-OS open RNA-seq count files from GDC API
OUT=/d/projects/OS_matrix_dualKO/data/raw/TARGET-OS/files
mkdir -p "$OUT"
ok=0; fail=0
for id in $(C:/Users/Administrator/.workbuddy/binaries/node/versions/22.22.2-3/node.exe -e "const j=require('D:/projects/OS_matrix_dualKO/data/raw/TARGET-OS/gdc_files_manifest.json');j.data.hits.forEach(h=>console.log(h.file_id))"); do
  if [ -s "$OUT/$id.tsv.gz" ] || [ -s "$OUT/$id.tsv" ]; then ok=$((ok+1)); continue; fi
  curl -s --max-time 300 --retry 3 --retry-delay 5 "https://api.gdc.cancer.gov/data/$id" -o "$OUT/$id.tsv"
  if [ -s "$OUT/$id.tsv" ] && head -1 "$OUT/$id.tsv" | grep -q "gene_id"; then
    ok=$((ok+1)); echo "OK $id"
  else
    fail=$((fail+1)); echo "FAIL $id"; rm -f "$OUT/$id.tsv"
  fi
done
echo "DONE ok=$ok fail=$fail" > /d/projects/OS_matrix_dualKO/logs/gdc_download_status.txt
