#!/bin/bash
# 安装 CellChat（GitHub 源码，需编译）+ 其缺失依赖
set -u
RT="/d/rtools45"
RSCRIPT="/d/R/R-4.6.1/bin/x64/Rscript.exe"
RCMD="/d/R/R-4.6.1/bin/x64/Rcmd.exe"
LIB="D:/R/R-4.6.1/library"
PKGSRC="$RT/tmp/pkgsrc"
TUNA="https://mirrors.tuna.tsinghua.edu.cn"
CRAN="$TUNA/CRAN"
BIOCVER="3.23"

export PATH="$RT/x86_64-w64-mingw32.static.posix/bin:$RT/ucrt64/bin:$RT/usr/bin:/d/R/R-4.6.1/bin/x64:/c/Program Files/JAGS/JAGS-4.3.2/x64/bin:$PATH"
export TMPDIR="$RT/tmp"
export MSYSTEM=UCRT64
mkdir -p "$PKGSRC"

is_installed() { "$RSCRIPT" -e "quit(status=!requireNamespace('$1', quietly=TRUE))" >/dev/null 2>&1; }

install_one() {
  local PKG="$1"
  if is_installed "$PKG"; then echo "[$PKG] already installed"; return 0; fi
  # CRAN 二进制
  "$RSCRIPT" -e "options(repos=c(CRAN='$CRAN')); install.packages('$PKG', lib='$LIB', type='win.binary', quiet=TRUE)" >/dev/null 2>&1
  is_installed "$PKG" && return 0
  # CRAN 源码
  "$RSCRIPT" -e "options(repos=c(CRAN='$CRAN')); install.packages('$PKG', lib='$LIB', type='source', INSTALL_opts='--no-multiarch', quiet=TRUE)" >/dev/null 2>&1
  is_installed "$PKG" && return 0
  # Bioc 源码
  for r in "$BIOCVER/bioc" "$BIOCVER/data/experiment" "$BIOCVER/data/annotation"; do
    local pk="$TUNA/bioconductor/packages/$r/src/contrib/PACKAGES"
    local ver=$(timeout 40 curl -s "$pk" 2>/dev/null | grep -A1 "^Package: ${PKG}$" | grep "^Version:" | head -1 | sed 's/Version: *//')
    if [ -n "$ver" ]; then
      cd "$PKGSRC" || return 1
      curl -sL -o "${PKG}.tar.gz" "$TUNA/bioconductor/packages/$r/src/contrib/${PKG}_${ver}.tar.gz"
      if tar tzf "${PKG}.tar.gz" >/dev/null 2>&1; then
        rm -rf "$PKG"; tar xzf "${PKG}.tar.gz"
        "$RCMD" INSTALL --no-multiarch -l "$LIB" "$PKG" >/dev/null 2>&1
        is_installed "$PKG" && return 0
      fi
    fi
  done
  return 1
}

echo "########## CellChat dependencies ##########"
for D in NMF ggalluvial sna network pacote reshape2 BiocNeighbors; do
  echo "--- dep: $D ---"
  install_one "$D" || echo "[dep $D] FAILED (may still be optional)"
done

echo "########## CellChat itself ##########"
if is_installed CellChat; then
  echo "CellChat already installed"
else
  cd "$PKGSRC" || exit 1
  rm -rf CellChat CellChat-main CellChat-master
  echo "downloading CellChat from GitHub via gh-proxy ..."
  curl -sL --max-time 600 -o CellChat.tar.gz "https://gh-proxy.com/https://github.com/jinworks/CellChat/archive/refs/heads/main.tar.gz"
  if ! tar tzf CellChat.tar.gz >/dev/null 2>&1; then
    echo "main branch failed, trying master ..."
    curl -sL --max-time 600 -o CellChat.tar.gz "https://gh-proxy.com/https://github.com/jinworks/CellChat/archive/refs/heads/master.tar.gz"
  fi
  if tar tzf CellChat.tar.gz >/dev/null 2>&1; then
    tar xzf CellChat.tar.gz
    SRC=$(ls -d CellChat-* 2>/dev/null | head -1)
    echo "extracted: $SRC"
    "$RCMD" INSTALL --no-multiarch -l "$LIB" "$SRC" 2>&1 | tail -25
    if is_installed CellChat; then echo "[CellChat] OK"; else echo "[CellChat] FAILED"; fi
  else
    echo "[CellChat] download/extract FAILED"
  fi
fi
