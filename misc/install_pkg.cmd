@echo off
REM Rtools45 MSYS2 R package installer launcher
REM Usage: install_pkg.cmd <pkg1> [pkg2 ...]
setlocal
set "MSYSTEM=UCRT64"
set "CHERE_INVOKING=1"
D:\rtools45\usr\bin\bash.exe -lc "/d/projects/OS_matrix_dualKO/code/install_r_pkg.sh %*"
endlocal
