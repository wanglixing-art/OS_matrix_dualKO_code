#!/bin/bash
# P0 v2: parallel range downloader for slow per-connection links
# Each file: split into N parts (10) downloaded concurrently, then concat.
# File queue is processed sequentially to cap total connections at ~10.

FTP=https://ftp.ncbi.nlm.nih.gov/geo/series
GDC=https://api.gdc.cancer.gov/data
BASE=/d/projects/OS_matrix_dualKO/data/raw
LOG=/d/projects/OS_matrix_dualKO/logs/p0v2.log
NPAR=10

pdl() {  # pdl <url> <outfile>
  local url=$1 out=$2
  local size clen parts i off end
  if [ -s "$out.done" ]; then echo "SKIP $(basename $out)" >> "$LOG"; return 0; fi
  clen=$(curl -sI --max-time 60 "$url" | grep -i content-length | tail -1 | tr -dc '0-9')
  size=${clen:-0}
  if [ "$size" -lt 100000 ]; then
    curl -s --max-time 600 --retry 3 "$url" -o "$out" && echo "SMALLOK $(basename $out) ($size B)" >> "$LOG" && touch "$out.done" && return 0
    echo "SMALLFAIL $(basename $out)" >> "$LOG"; return 1
  fi
  parts=$NPAR
  local pids=() chunk=$(( size / parts + 1 ))
  for i in $(seq 0 $((parts-1))); do
    off=$(( i * chunk )); end=$(( off + chunk - 1 ))
    [ $end -ge $size ] && end=$(( size - 1 ))
    curl -s --max-time 7200 --retry 5 --retry-delay 5 -r "$off-$end" "$url" -o "$out.part$i" &
    pids+=($!)
  done
  local ok=1
  for p in "${pids[@]}"; do wait "$p" || ok=0; done
  # verify part sizes
  local got=0
  for i in $(seq 0 $((parts-1))); do
    [ -f "$out.part$i" ] && got=$(( got + $(stat -c %s "$out.part$i") ))
  done
  if [ "$ok" = "1" ] && [ "$got" -eq "$size" ]; then
    cat $(for i in $(seq 0 $((parts-1))); do echo "$out.part$i"; done) > "$out"
    rm -f "$out".part*
    touch "$out.done"
    echo "OK $(basename $out) ($size B)" >> "$LOG"
  else
    echo "FAIL $(basename $out) got=$got want=$size" >> "$LOG"; return 1
  fi
}

echo "=== P0v2 start $(date) ===" >> "$LOG"

# --- GDC: 88 count files, file-level parallelism 8 at a time ---
mkdir -p "$BASE/TARGET-OS/files"
ids=$(C:/Users/Administrator/.workbuddy/binaries/node/versions/22.22.2-3/node.exe -e "
const j=require('D:/projects/OS_matrix_dualKO/data/raw/TARGET-OS/gdc_files_manifest.json');
j.data.hits.forEach(h=>console.log(h.file_id))")
echo "$ids" | xargs -P 8 -I {} bash -c '
  id="{}"; out="'"$BASE"'/TARGET-OS/files/{}.tsv"
  if [ -s "$out.done" ]; then exit 0; fi
  for a in 1 2 3; do
    curl -s --max-time 900 --retry 3 "'"$GDC"'/$id" -o "$out"
    if [ -s "$out" ] && head -1 "$out" | grep -q "gene_id"; then touch "$out.done"; echo "GDOK $id" >> '"$LOG"'; exit 0; fi
    sleep 3
  done
  echo "GDFAIL $id" >> '"$LOG"'
'
echo "=== GDC done $(date) ===" >> "$LOG"

# --- GEO series matrices (small) ---
pdl "$FTP/GSE16nnn/GSE16088/matrix/GSE16088_series_matrix.txt.gz"     "$BASE/GSE16088/series_matrix.txt.gz"
pdl "$FTP/GSE21nnn/GSE21257/matrix/GSE21257_series_matrix.txt.gz"     "$BASE/GSE21257/series_matrix.txt.gz"
pdl "$FTP/GSE162nnn/GSE162454/matrix/GSE162454_series_matrix.txt.gz"  "$BASE/GSE162454/series_matrix.txt.gz"
pdl "$FTP/GSE152nnn/GSE152048/matrix/GSE152048_series_matrix.txt.gz"  "$BASE/GSE152048/series_matrix.txt.gz"

# --- scRNA raw data (range-parallel) ---
pdl "$FTP/GSE162nnn/GSE162454/suppl/GSE162454_RAW.tar" "$BASE/GSE162454/GSE162454_RAW.tar"
for bc in BC2 BC22 BC10 BC17; do
  pdl "$FTP/GSE152nnn/GSE152048/suppl/GSE152048_${bc}.matrix.tar.gz" "$BASE/GSE152048/GSE152048_${bc}.matrix.tar.gz"
done

echo "=== P0v2 ALLDONE $(date) ===" >> "$LOG"
