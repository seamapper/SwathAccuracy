@echo off
cd /d "%~dp0"
set "PYTHON_EXE=%LOCALAPPDATA%\miniforge3\python.exe"
if not exist "%PYTHON_EXE%" (
  echo Miniforge Python not found at %PYTHON_EXE%
  pause
  exit /b 1
)
echo Using Python: %PYTHON_EXE%
"%PYTHON_EXE%" build_exe.py
pause
