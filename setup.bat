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

rem #############################################
rem # [4] プロジェクトを構成 ^& ビルド
rem #############################################
echo.
echo Configuring the project with CMake...
cmake --preset=default

@REM if errorlevel 1 (
@REM     echo ERROR: CMake configuration failed.
@REM     pause
@REM     goto :eof
@REM )

@REM echo.
@REM echo Building the project...
@REM cmake --build build

@REM if errorlevel 1 (
@REM     echo ERROR: Build failed.
@REM     pause
@REM     goto :eof
@REM )

@REM echo.
@REM echo --- Build process finished successfully. ---
@REM pause