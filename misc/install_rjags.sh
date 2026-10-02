#!/bin/bash
# 在 JAGS 已安装的前提下编译安装 rjags（需能找到 libjags-4.dll）
set -u
RT="/d/rtools45"
RSCRIPT="/d/R/R-4.6.1/bin/x64/Rscript.exe"
RCMD="/d/R/R-4.6.1/bin/x64/Rcmd.exe"
LIB="D:/R/R-4.6.1/library"
PKGSRC="$RT/tmp/pkgsrc"
CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"

export PATH="$RT/x86_64-w64-mingw32.static.posix/bin:$RT/ucrt64/bin:$RT/usr/bin:/d/R/R-4.6.1/bin/x64:/c/Program Files/JAGS/JAGS-4.3.2/x64/bin:$PATH"
export TMPDIR="$RT/tmp"
export MSYSTEM=UCRT64
export JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.2"
# rjags 配置脚本靠 JAGS_ROOT / JAGS_HOME 定位
export JAGS_ROOT="$JAGS_HOME"
mkdir -p "$PKGSRC"

echo "=== verify JAGS dll visible ==="
ls -la "/c/Program Files/JAGS/JAGS-4.3.2/x64/bin/libjags-4.dll"

for PKG in rjags "$@"; do
  echo "========== [$PKG] =========="
  if "$RSCRIPT" -e "quit(status=!requireNamespace('$PKG', quietly=TRUE))" >/dev/null 2>&1; then
    echo "[$PKG] already installed, skip"; continue
  fi
  cd "$PKGSRC" || exit 1
  # rjags 走源码编译（需要 JAGS）
  "$RSCRIPT" -e "options(repos=c(CRAN='$CRAN')); install.packages('$PKG', lib='$LIB', type='source', INSTALL_opts='--no-multiarch', quiet=TRUE)" 2>&1 | tail -12
  if "$RSCRIPT" -e "quit(status=!requireNamespace('$PKG', quietly=TRUE))" >/dev/null 2>&1; then
    echo "[$PKG] OK"
  else
    echo "[$PKG] FAILED"
  fi
done
