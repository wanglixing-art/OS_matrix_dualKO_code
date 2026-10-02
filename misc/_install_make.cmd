@echo off
set "PATH=D:\rtools45\usr\bin;%PATH%"
D:\rtools45\usr\bin\pacman.exe -S --noconfirm --needed mingw-w64-ucrt-x86_64-make
echo PACMAN_EXIT=%ERRORLEVEL%
