#!/bin/bash
# P0: Download GEO datasets (bulk + scRNA) with retry loop
BASE=/d/projects/OS_matrix_dualKO/data/raw
FTP=https://ftp.ncbi.nlm.nih.gov/geo/series

dl() {  # dl <url> <outfile>
  local url=$1 out=$2 n=0
  while [ $n -lt 3 ]; do
    curl -s --max-time 3600 -C - "$url" -o "$out" && return 0
    n=$((n+1)); sleep 5
  done
  return 1
}

# bulk datasets
dl "$FTP/GSE16nnn/GSE16088/suppl/GSE16088_RAW.tar"          "$BASE/GSE16088/GSE16088_RAW.tar"   && echo OK16088
dl "$FTP/GSE21nnn/GSE21257/suppl/GSE21257_RAW.tar"          "$BASE/GSE21257/GSE21257_RAW.tar"   && echo OK21257
dl "$FTP/GSE162nnn/GSE162454/suppl/GSE162454_RAW.tar"       "$BASE/GSE162454/GSE162454_RAW.tar" && echo OK162454
dl "$FTP/GSE16nnn/GSE16088/matrix/GSE16088_series_matrix.txt.gz"     "$BASE/GSE16088/series_matrix.txt.gz"     && echo OK16088meta
dl "$FTP/GSE21nnn/GSE21257/matrix/GSE21257_series_matrix.txt.gz"     "$BASE/GSE21257/series_matrix.txt.gz"     && echo OK21257meta
dl "$FTP/GSE162nnn/GSE162454/matrix/GSE162454_series_matrix.txt.gz"  "$BASE/GSE162454/series_matrix.txt.gz"    && echo OK162454meta

# GSE152048: only the 4 samples used in the metastasis paper (BC2/BC22 primary, BC10/BC17 metastasis)
for bc in BC2 BC22 BC10 BC17; do
  dl "$FTP/GSE152nnn/GSE152048/suppl/GSE152048_${bc}.matrix.tar.gz" "$BASE/GSE152048/GSE152048_${bc}.matrix.tar.gz" && echo "OK152048_$bc"
done
dl "$FTP/GSE152nnn/GSE152048/matrix/GSE152048_series_matrix.txt.gz" "$BASE/GSE152048/series_matrix.txt.gz" && echo OK152048meta

echo ALLDONE
