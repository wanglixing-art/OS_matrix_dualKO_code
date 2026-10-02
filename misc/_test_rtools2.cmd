@echo off
setlocal
set "PATH=D:\rtools45\usr\bin;D:\rtools45\x86_64-w64-mingw32.static.posix\bin;D:\R\R-4.6.1\bin\x64;%PATH%"
set "MSYS=noconv"
set "MSYS2_ARG_CONV_EXCL=*"
set "BINPREF=D:/rtools45/x86_64-w64-mingw32.static.posix/bin/"
set "MAKEFLAGS=--no-print-directory"
"D:\R\R-4.6.1\bin\x64\Rscript.exe" "D:\projects\OS_matrix_dualKO\code\_test_rtools.R" > "D:\projects\OS_matrix_dualKO\logs\_test_rtools2.log" 2>&1
echo EXIT=%ERRORLEVEL%
type "D:\projects\OS_matrix_dualKO\logs\_test_rtools2.log"
endlocal
