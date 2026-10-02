@echo off
setlocal
set "MSYSTEM=UCRT64"
set "CHERE_INVOKING=1"
call "D:\rtools45\msys2_shell.cmd" -defterm -no-start -ucrt64 -here -c "cd /d/projects/OS_matrix_dualKO && R=/d/R/R-4.6.1/bin/x64/Rscript.exe && $R code/_test_rtools.R"
endlocal
