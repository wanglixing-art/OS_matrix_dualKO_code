#!/bin/bash
# ============================================================
# Rtools45 (MSYS2) 环境下的 R 包安装器 —— 支持 CRAN / Bioconductor 多仓库
# 用法: install_r_pkg.sh <包名> [更多包名...]
#
# 必须在 Rtools45 的 MSYS2 bash 中运行，PATH 需含
#   /d/rtools45/x86_64-w64-mingw32.static.posix/bin  (g++/make 真实来源)
# ============================================================
set -u

RT="/d/rtools45"
RSCRIPT="/d/R/R-4.6.1/bin/x64/Rscript.exe"
RCMD="/d/R/R-4.6.1/bin/x64/Rcmd.exe"
LIB="D:/R/R-4.6.1/library"
PKGSRC="$RT/tmp/pkgsrc"
TUNA="https://mirrors.tuna.tsinghua.edu.cn"
CRAN="$TUNA/CRAN"
BIOCVER="3.23"

export PATH="$RT/x86_64-w64-mingw32.static.posix/bin:$RT/ucrt64/bin:$RT/usr/bin:/d/R/R-4.6.1/bin/x64:$PATH"
export TMPDIR="$RT/tmp"
export MSYSTEM=UCRT64
mkdir -p "$PKGSRC"

is_installed() { "$RSCRIPT" -e "quit(status=!requireNamespace('$1', quietly=TRUE))" >/dev/null 2>&1; }

# 找包源码 tarball 的下载 URL（跨仓库）
find_tarball_url() {
  local PKG="$1"
  local repos="$BIOCVER/bioc $BIOCVER/data/experiment $BIOCVER/data/annotation"
  for r in $repos; do
    local pk="$TUNA/bioconductor/packages/$r/src/contrib/PACKAGES"
    local line
    line=$(timeout 40 curl -s "$pk" 2>/dev/null | grep -A1 "^Package: ${PKG}$" | grep "^Version:" | head -1)
    if [ -n "$line" ]; then
      local ver=$(echo "$line" | sed 's/Version: *//')
      echo "$TUNA/bioconductor/packages/$r/src/contrib/${PKG}_${ver}.tar.gz"
      return 0
    fi
  done
  return 1
}

FAILED=""
for PKG in "$@"; do
  echo "========== [$PKG] =========="
  if is_installed "$PKG"; then echo "[$PKG] already installed, skip"; continue; fi

  ok=0
  # ---- 1) CRAN 二进制 ----
  "$RSCRIPT" -e "options(repos=c(CRAN='$CRAN')); install.packages('$PKG', lib='$LIB', type='win.binary', quiet=TRUE)" >/dev/null 2>&1
  is_installed "$PKG" && ok=1

  # ---- 2) CRAN 源码 ----
  if [ "$ok" -eq 0 ]; then
    echo "[$PKG] -> CRAN source"
    "$RSCRIPT" -e "options(repos=c(CRAN='$CRAN')); install.packages('$PKG', lib='$LIB', type='source', INSTALL_opts='--no-multiarch', quiet=TRUE)" >/dev/null 2>&1
    is_installed "$PKG" && ok=1
  fi

  # ---- 3) Bioconductor 源码（多仓库解析 tarball）----
  if [ "$ok" -eq 0 ]; then
    echo "[$PKG] -> Bioconductor source"
    cd "$PKGSRC" || exit 1
    URL=$(find_tarball_url "$PKG")
    if [ -n "$URL" ]; then
      echo "  url: $URL"
      curl -sL -o "${PKG}.tar.gz" "$URL"
      if [ -s "${PKG}.tar.gz" ] && tar tzf "${PKG}.tar.gz" >/dev/null 2>&1; then
        rm -rf "$PKG"; tar xzf "${PKG}.tar.gz"
        "$RCMD" INSTALL --no-multiarch -l "$LIB" "$PKG" 2>&1 | tail -18
        is_installed "$PKG" && ok=1
      fi
    else
      echo "  tarball URL not found in Bioc $BIOCVER"
    fi
  fi

  if [ "$ok" -eq 1 ]; then echo "[$PKG] OK"; else echo "[$PKG] FAILED"; FAILED="$FAILED $PKG"; fi
done

echo "=============================="
if [ -n "$FAILED" ]; then echo "FAILED PACKAGES:$FAILED"; exit 1; else echo "ALL PACKAGES OK"; fi
