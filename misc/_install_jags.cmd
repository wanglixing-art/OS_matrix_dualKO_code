@echo off
REM Silent install JAGS 4.3.2
D:\rtools45\tmp\jags.exe /S /D=C:\Program Files\JAGS
echo JAGS_INSTALL_EXIT=%ERRORLEVEL%
timeout /t 20 /nobreak >nul
if exist "C:\Program Files\JAGS\JAGS-4.3.2\bin\jags.exe" (echo JAGS_FOUND) else (echo JAGS_MISSING)
if exist "C:\Program Files\JAGS\JAGS-4.3.2\x64\bin\jags.exe" (echo JAGS_X64_FOUND) else (echo JAGS_X64_MISSING)
dir /b "C:\Program Files\JAGS" 2>nul
