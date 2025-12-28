@echo off
setlocal enabledelayedexpansion

rem ============================================
rem TBP Project - Main Setup Script
rem ============================================
rem Usage: scripts\setup.bat
rem 
rem This script will:
rem   1. Check dependencies (git, cmake, ninja, oneAPI, MinGW)
rem   2. Install vcpkg and Boost if needed
rem   3. Set up environment variables (BOOST_ROOT)
rem   4. Configure and build the project
rem ============================================

echo.
echo ============================================================
echo   TBP Project Setup
echo ============================================================
echo   This script will set up your development environment
echo   and build the project.
echo ============================================================
echo.

rem Get script directory and project root
set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."
cd /d "!PROJECT_ROOT!"

echo Project root: !PROJECT_ROOT!
echo.

rem ---------------------------------------------
rem [1] Check Dependencies
rem ---------------------------------------------
echo ============================================================
echo   Step 1: Checking Dependencies
echo ============================================================
echo.

call "!SCRIPT_DIR!check_dependencies.bat"
if errorlevel 1 (
    echo.
    echo [ERROR] Dependency check failed. Please install missing required dependencies.
    pause
    exit /b 1
)

rem Variables are exported from check_dependencies.bat:
rem   HAS_ONEAPI, HAS_NINJA, HAS_MINGW, HAS_BOOST, HAS_VCPKG

rem ---------------------------------------------
rem [2] Install vcpkg and Boost if needed
rem ---------------------------------------------
echo.
echo ============================================================
echo   Step 2: Setting up Boost Library
echo ============================================================
echo.

if "!HAS_BOOST!"=="0" (
    echo Boost is not installed. Installing via vcpkg...
    echo.
    call "!SCRIPT_DIR!install_vcpkg_boost.bat"
    if errorlevel 1 (
        echo.
        echo [ERROR] Failed to install Boost.
        pause
        exit /b 1
    )
) else (
    echo Boost is already available.
    echo BOOST_ROOT: !BOOST_ROOT!
)

rem Ensure BOOST_ROOT is set for this session
if not defined BOOST_ROOT (
    set "VCPKG_ROOT=!PROJECT_ROOT!\external\vcpkg"
    set "BOOST_ROOT=!VCPKG_ROOT!\installed\x64-windows"
)

rem ---------------------------------------------
rem [3] Select Build Configuration
rem ---------------------------------------------
echo.
echo ============================================================
echo   Step 3: Select Build Configuration
echo ============================================================
echo.

set "BUILD_PRESET="
set "USE_ONEAPI=0"

rem Determine available options
if "!HAS_ONEAPI!"=="1" if "!HAS_NINJA!"=="1" (
    set "ONEAPI_AVAILABLE=1"
) else (
    set "ONEAPI_AVAILABLE=0"
)

if "!ONEAPI_AVAILABLE!"=="1" if "!HAS_MINGW!"=="1" (
    rem Both options available - ask user
    echo Available build configurations:
    echo   1. Intel oneAPI + Ninja [Recommended - best performance]
    echo   2. MinGW-gcc + Ninja [Alternative - good compatibility]
    echo.
    set /p "BUILD_CHOICE=Select configuration [1 or 2]: "
)

rem Process choice using GOTO to avoid nested if parsing issues
if "!ONEAPI_AVAILABLE!"=="0" if "!HAS_MINGW!"=="0" goto :no_compiler
if "!ONEAPI_AVAILABLE!"=="0" goto :use_mingw
if "!HAS_MINGW!"=="0" goto :use_oneapi

rem Both available - check user choice
if "!BUILD_CHOICE!"=="2" goto :use_mingw
goto :use_oneapi

:use_oneapi
echo Using Intel oneAPI + Ninja configuration.
set "BUILD_PRESET=ninja-intel"
set "USE_ONEAPI=1"
goto :compiler_selected

:use_mingw
echo Using MinGW-gcc configuration.
set "BUILD_PRESET=mingw-gcc"
set "USE_ONEAPI=0"
goto :compiler_selected

:no_compiler
echo [ERROR] No suitable compiler found!
echo         Please install Intel oneAPI or MinGW-gcc.
pause
exit /b 1

:compiler_selected

echo.
echo Selected preset: !BUILD_PRESET!
echo.

rem ---------------------------------------------
rem [4] Set up Compiler Environment
rem ---------------------------------------------
echo ============================================================
echo   Step 4: Setting up Compiler Environment
echo ============================================================
echo.

if "!USE_ONEAPI!"=="1" (
    echo Setting up Visual Studio environment for Intel oneAPI...
    
    rem Intel oneAPI requires MSVC linker - call vcvars64.bat first
    set "VCVARS_PATH="
    if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
    ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"
    ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
    ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
    ) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
    ) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"
    ) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
    ) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
        set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
    )
    
    if defined VCVARS_PATH (
        echo   Calling: !VCVARS_PATH!
        call "!VCVARS_PATH!"
    ) else (
        echo [WARNING] vcvars64.bat not found. Intel oneAPI may not work correctly.
    )
    
    echo Setting up Intel oneAPI environment...
    call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
    if errorlevel 1 (
        echo [WARNING] Failed to set up oneAPI environment.
        echo           Falling back to MinGW-gcc if available.
        if "!HAS_MINGW!"=="1" (
            set "BUILD_PRESET=mingw-gcc"
            set "USE_ONEAPI=0"
        ) else (
            echo [ERROR] No fallback compiler available.
            pause
            exit /b 1
        )
    )
) else (
    echo Using MinGW-gcc compiler.
    echo g++ path:
    where g++
)

rem ---------------------------------------------
rem [5] Clean and Configure
rem ---------------------------------------------
echo.
echo ============================================================
echo   Step 5: Configuring Project with CMake
echo ============================================================
echo.

rem Ask about cleaning build directory
if exist "!PROJECT_ROOT!\build" (
    echo Build directory exists.
    set /p "CLEAN_BUILD=Do you want to clean and reconfigure? (Y/N): "
    if /i "!CLEAN_BUILD!"=="Y" (
        echo Deleting old build directory...
        rmdir /s /q "!PROJECT_ROOT!\build"
        if exist "!PROJECT_ROOT!\CMakeCache.txt" del "!PROJECT_ROOT!\CMakeCache.txt"
    )
)

echo.
echo Configuring project with preset: !BUILD_PRESET!
echo BOOST_ROOT: !BOOST_ROOT!
echo.

cmake --preset=!BUILD_PRESET!
if errorlevel 1 (
    echo.
    echo [ERROR] CMake configuration failed.
    echo.
    echo Troubleshooting tips:
    echo   1. Ensure BOOST_ROOT is set correctly: !BOOST_ROOT!
    echo   2. Check that Boost headers exist at: !BOOST_ROOT!\include\boost
    echo   3. Try running this script again after fixing issues.
    pause
    exit /b 1
)

rem ---------------------------------------------
rem [6] Build Project
rem ---------------------------------------------
echo.
echo ============================================================
echo   Step 6: Building Project
echo ============================================================
echo.

cmake --build build
if errorlevel 1 (
    echo.
    echo [ERROR] Build failed.
    pause
    exit /b 1
)

rem ---------------------------------------------
rem [7] Summary
rem ---------------------------------------------
echo.
echo ============================================================
echo   Setup Complete!
echo ============================================================
echo.
echo   Build Configuration: !BUILD_PRESET!
echo   BOOST_ROOT: !BOOST_ROOT!
echo   Build Directory: !PROJECT_ROOT!\build
echo   Executables: !PROJECT_ROOT!\build\bin
echo.
echo   To rebuild later:
echo     cmake --build build
echo.
echo   To reconfigure:
echo     cmake --preset=!BUILD_PRESET!
echo.
echo ============================================================
echo.

pause
endlocal
exit /b 0
