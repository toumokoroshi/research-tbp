@echo off
setlocal enabledelayedexpansion


rem #############################################
rem # [2] oneAPI環境をセットアップ
rem #############################################
echo Setting up Intel oneAPI environment...
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

if errorlevel 1 (
    echo ERROR: Failed to set up oneAPI environment.
    pause
    goto :eof
)

rem #############################################
rem # [3] ビルドディレクトリをクリーンアップ
rem #############################################
echo.
echo Deleting old build directory...
if exist build rmdir /s /q build
if exist CMakeCache.txt del CMakeCache.txt

rem #############################################
rem # [4] プロジェクトを構成 ^& ビルド
rem #############################################
echo.
echo Configuring the project with CMake...
cmake --preset=vs-intel

if errorlevel 1 (
    echo ERROR: CMake configuration failed.
    pause
    goto :eof
)

echo.
echo Building the project...
cmake --build build

if errorlevel 1 (
    echo ERROR: Build failed.
    pause
    goto :eof
)

echo.
echo --- Build process finished successfully. ---
pause