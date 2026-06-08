@echo off
setlocal enabledelayedexpansion

REM Build a Windows .exe launcher for the full Celltaminate workflow.
REM Run from the repository root on a Windows machine with Python installed:
REM   windows\build_windows_exe.bat

cd /d "%~dp0\.."

python --version >nul 2>&1
if errorlevel 1 (
  echo ERROR: Python was not found. Install Python 3.9 or newer from https://www.python.org/downloads/windows/.
  exit /b 1
)

python -m pip install --upgrade pip
python -m pip install --upgrade pyinstaller

python -m PyInstaller ^
  --onefile ^
  --name CelltaminateRunner ^
  scripts\celltaminate_full_workflow.py

if errorlevel 1 (
  echo ERROR: PyInstaller build failed.
  exit /b 1
)

echo.
echo Built: dist\CelltaminateRunner.exe
echo Test with: dist\CelltaminateRunner.exe --help
endlocal
