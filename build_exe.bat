@echo off
cd /d "%~dp0"
"%~dp0..\.venv\Scripts\python.exe" build_exe.py
pause
