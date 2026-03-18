@echo off
setlocal enabledelayedexpansion

rem ============================================
rem TBP Project - Main Setup Script
rem ============================================
rem Usage: tools\setup.bat [options]
rem 
rem Options:
rem   --unattended, -y    Non-interactive mode (auto-install all)
rem   --with-oneapi       Install Intel oneAPI (default behavior)
rem   --without-oneapi    Skip Intel oneAPI installation
rem   --skip-python       Skip Python environment setup
rem   --python-venv       Install Python dependencies into .venv
rem   --skip-test         Skip running tests after build
rem   --help, -h          Show this help message
rem 
rem This script will:
rem   1. Check dependencies (git, cmake, ninja, oneAPI, MinGW)
rem   2. Install vcpkg and Boost if needed
rem   3. Set up environment variables (BOOST_ROOT)
rem   4. Configure and build the project
rem   5. (Optional) Set up Python environment
rem   6. (Optional) Run tests to verify build
rem ============================================

rem ---------------------------------------------
rem Parse command line arguments
rem ---------------------------------------------
set "UNATTENDED=0"
set "WITH_ONEAPI=1"
set "SKIP_PYTHON=0"
set "PYTHON_VENV=0"
set "SKIP_TEST=0"

:parse_args
if "%~1"=="" goto :args_done
if /i "%~1"=="--unattended" set "UNATTENDED=1" & shift & goto :parse_args
if /i "%~1"=="-y" set "UNATTENDED=1" & shift & goto :parse_args
if /i "%~1"=="--with-oneapi" set "WITH_ONEAPI=1" & shift & goto :parse_args
if /i "%~1"=="--without-oneapi" set "WITH_ONEAPI=0" & shift & goto :parse_args
if /i "%~1"=="--skip-python" set "SKIP_PYTHON=1" & shift & goto :parse_args
if /i "%~1"=="--python-venv" set "PYTHON_VENV=1" & shift & goto :parse_args
if /i "%~1"=="--skip-test" set "SKIP_TEST=1" & shift & goto :parse_args
if /i "%~1"=="--help" goto :show_help
if /i "%~1"=="-h" goto :show_help
echo [WARNING] Unknown option: %~1
shift
goto :parse_args

:show_help
echo.
echo Usage: tools\setup.bat [options]
echo.
echo Options:
echo   --unattended, -y    Non-interactive mode (auto-install all)
echo   --with-oneapi       Install Intel oneAPI (default, best performance)
echo   --without-oneapi    Skip Intel oneAPI installation
echo   --skip-python       Skip Python environment setup
echo   --python-venv       Install Python dependencies into .venv
echo   --skip-test         Skip running tests after build
echo   --help, -h          Show this help message
echo.
echo Examples:
echo   tools\setup.bat                    Interactive setup
echo   tools\setup.bat --unattended       Fully automated setup (oneAPI preferred)
echo   tools\setup.bat -y --without-oneapi Auto setup without Intel oneAPI
echo   tools\setup.bat -y --python-venv   Auto setup with isolated Python env
echo   tools\setup.bat -y --skip-test     Auto setup without tests
echo.
exit /b 0

:args_done

echo.
echo ============================================================
echo   TBP Project Setup
echo ============================================================
echo   This script will set up your development environment
echo   and build the project.
if "!UNATTENDED!"=="1" (
echo   Mode: UNATTENDED (non-interactive)
) else (
echo   Mode: Interactive
)
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

rem Pass UNATTENDED and WITH_ONEAPI flags to check_dependencies.bat
call "!SCRIPT_DIR!check_dependencies.bat" !UNATTENDED! !WITH_ONEAPI!
if errorlevel 1 (
    echo.
    echo [ERROR] Dependency check failed. Please install missing required dependencies.
    if "!UNATTENDED!"=="0" pause
    exit /b 1
)

rem Variables are exported from check_dependencies.bat:
rem   HAS_ONEAPI, HAS_NINJA, HAS_MINGW, HAS_BOOST, HAS_VCPKG

rem ---------------------------------------------
rem [1.5] Add winget portable packages to PATH (e.g., Ninja)
rem ---------------------------------------------
rem winget portable packages are not added to PATH automatically
set "WINGET_PACKAGES=%LOCALAPPDATA%\Microsoft\WinGet\Packages"
if exist "!WINGET_PACKAGES!" (
    rem Search for Ninja in winget packages
    for /d %%D in ("!WINGET_PACKAGES!\Ninja-build.Ninja_*") do (
        if exist "%%D\ninja.exe" (
            echo Adding Ninja to PATH from winget: %%D
            set "PATH=%%D;!PATH!"
            set "HAS_NINJA=1"
        )
    )
)

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
    call "!SCRIPT_DIR!install_vcpkg_boost.bat" !UNATTENDED!
    if errorlevel 1 (
        echo.
        echo [ERROR] Failed to install Boost.
        if "!UNATTENDED!"=="0" pause
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

rem In unattended mode, auto-select best available compiler
if "!UNATTENDED!"=="1" (
    if "!ONEAPI_AVAILABLE!"=="1" (
        echo [AUTO] Selecting Intel oneAPI + Ninja
        goto :use_oneapi
    ) else if "!HAS_MINGW!"=="1" if "!HAS_NINJA!"=="1" (
        echo [AUTO] Selecting MinGW-gcc + Ninja (Intel oneAPI not available)
        goto :use_mingw_ninja
    ) else if "!HAS_MINGW!"=="1" (
        echo [AUTO] Selecting MinGW-gcc + MinGW Makefiles (Ninja not available)
        goto :use_mingw_make
    ) else (
        goto :no_compiler
    )
)

rem Interactive mode - ask user if both are available
if "!ONEAPI_AVAILABLE!"=="1" if "!HAS_MINGW!"=="1" (
    rem Both options available - ask user
    echo Available build configurations:
    echo   1. Intel oneAPI + Ninja [Recommended]
    if "!HAS_NINJA!"=="1" (
        echo   2. MinGW-gcc + Ninja [Alternative]
    ) else (
        echo   2. MinGW-gcc + MinGW Makefiles [No Ninja required]
    )
    echo.
    set /p "BUILD_CHOICE=Select configuration [1 or 2]: "
)

rem Process choice using GOTO to avoid nested if parsing issues
if "!ONEAPI_AVAILABLE!"=="0" if "!HAS_MINGW!"=="0" goto :no_compiler
if "!ONEAPI_AVAILABLE!"=="0" if "!HAS_NINJA!"=="1" goto :use_mingw_ninja
if "!ONEAPI_AVAILABLE!"=="0" goto :use_mingw_make
if "!HAS_MINGW!"=="0" goto :use_oneapi

rem Both available - check user choice
if "!BUILD_CHOICE!"=="2" if "!HAS_NINJA!"=="1" goto :use_mingw_ninja
if "!BUILD_CHOICE!"=="2" goto :use_mingw_make
goto :use_oneapi

:use_oneapi
echo Using Intel oneAPI + Ninja configuration.
set "BUILD_PRESET=ninja-intel"
set "USE_ONEAPI=1"
goto :compiler_selected

:use_mingw_ninja
echo Using MinGW-gcc + Ninja configuration.
set "BUILD_PRESET=mingw-gcc"
set "USE_ONEAPI=0"
goto :compiler_selected

:use_mingw_make
echo Using MinGW-gcc + MinGW Makefiles configuration.
set "BUILD_PRESET=mingw-gcc-make"
set "USE_ONEAPI=0"
goto :compiler_selected

:no_compiler
echo [ERROR] No suitable compiler found!
echo         Please install Intel oneAPI or MinGW-gcc.
if "!UNATTENDED!"=="0" pause
exit /b 1
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
    
    rem Search for vcvars64.bat in various VS installations
    rem VS 2022+ Insiders/Preview (Program Files)
    for /d %%V in ("C:\Program Files\Microsoft Visual Studio\*") do (
        for /d %%E in ("%%V\*") do (
            if exist "%%E\VC\Auxiliary\Build\vcvars64.bat" (
                if not defined VCVARS_PATH (
                    set "VCVARS_PATH=%%E\VC\Auxiliary\Build\vcvars64.bat"
                )
            )
        )
    )
    rem VS 2022 Standard editions (Program Files)
    if not defined VCVARS_PATH (
        if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
        ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"
        ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
        ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
        )
    )
    rem VS 2022 (Program Files x86)
    if not defined VCVARS_PATH (
        if exist "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
        )
    )
    rem VS 2019 (Program Files x86)
    if not defined VCVARS_PATH (
        if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
        ) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"
        ) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
        ) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
            set "VCVARS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
        )
    )
    
    if defined VCVARS_PATH (
        echo   Calling: !VCVARS_PATH!
        call "!VCVARS_PATH!"
        if errorlevel 1 (
            echo [WARNING] vcvars64.bat returned an error.
        )
    ) else (
        echo [WARNING] vcvars64.bat not found. Intel oneAPI may not work correctly.
        echo          Please install Visual Studio Build Tools with C++ workload.
    )
    
    echo Setting up Intel oneAPI environment...
    call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
    if errorlevel 1 (
        echo [WARNING] Failed to set up oneAPI environment.
        echo           Falling back to MinGW-gcc if available.
        if "!HAS_MINGW!"=="1" (
            if "!HAS_NINJA!"=="1" (
                set "BUILD_PRESET=mingw-gcc"
            ) else (
                set "BUILD_PRESET=mingw-gcc-make"
            )
            set "USE_ONEAPI=0"
        ) else (
            echo [ERROR] No fallback compiler available.
            if "!UNATTENDED!"=="0" pause
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

rem Ask about cleaning build directory (skip in unattended mode - keep existing)
if exist "!PROJECT_ROOT!\build" (
    if "!UNATTENDED!"=="1" (
        echo Build directory exists. Keeping existing configuration.
    ) else (
        echo Build directory exists.
        set /p "CLEAN_BUILD=Do you want to clean and reconfigure? (Y/N): "
        if /i "!CLEAN_BUILD!"=="Y" (
            echo Deleting old build directory...
            rmdir /s /q "!PROJECT_ROOT!\build"
            if exist "!PROJECT_ROOT!\CMakeCache.txt" del "!PROJECT_ROOT!\CMakeCache.txt"
        )
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
    if "!UNATTENDED!"=="0" pause
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
    if "!UNATTENDED!"=="0" pause
    exit /b 1
)

rem ---------------------------------------------
rem [7] Python Environment Setup (Optional)
rem ---------------------------------------------
if "!SKIP_PYTHON!"=="0" (
    echo.
    echo ============================================================
    echo   Step 7: Setting up Python Environment (Optional)
    echo ============================================================
    echo.
    
    where python >nul 2>&1
    if errorlevel 1 (
        echo [SKIP] Python is not installed. Skipping Python setup.
    ) else (
        if exist "!PROJECT_ROOT!\scripts\requirements.txt" (
            if "!PYTHON_VENV!"=="1" (
                set "VENV_DIR=!PROJECT_ROOT!\.venv"
                if not exist "!VENV_DIR!\Scripts\python.exe" (
                    echo Creating virtual environment: !VENV_DIR!
                    python -m venv "!VENV_DIR!"
                )
                if exist "!VENV_DIR!\Scripts\python.exe" (
                    echo Installing Python dependencies into .venv...
                    "!VENV_DIR!\Scripts\python.exe" -m pip install -r "!PROJECT_ROOT!\scripts\requirements.txt" --quiet
                    if errorlevel 1 (
                        echo [WARNING] Failed to install some Python packages in .venv.
                        echo          You can install them manually later:
                        echo          .venv\Scripts\python -m pip install -r scripts\requirements.txt
                    ) else (
                        echo [OK] Python dependencies installed in .venv.
                    )
                ) else (
                    echo [WARNING] Failed to create .venv. Falling back to global pip.
                    python -m pip install -r "!PROJECT_ROOT!\scripts\requirements.txt" --quiet
                    if errorlevel 1 (
                        echo [WARNING] Failed to install some Python packages.
                        echo          You can install them manually later:
                        echo          pip install -r scripts\requirements.txt
                    ) else (
                        echo [OK] Python dependencies installed.
                    )
                )
            ) else (
                echo Installing Python dependencies...
                python -m pip install -r "!PROJECT_ROOT!\scripts\requirements.txt" --quiet
                if errorlevel 1 (
                    echo [WARNING] Failed to install some Python packages.
                    echo          You can install them manually later:
                    echo          pip install -r scripts\requirements.txt
                ) else (
                    echo [OK] Python dependencies installed.
                )
            )
        ) else (
            echo [SKIP] No requirements.txt found.
        )
    )
) else (
    echo.
    echo [SKIP] Python setup skipped (--skip-python).
)

rem ---------------------------------------------
rem [8] Run Tests (Optional)
rem ---------------------------------------------
if "!SKIP_TEST!"=="0" (
    echo.
    echo ============================================================
    echo   Step 8: Running Tests
    echo ============================================================
    echo.
    
    where ctest >nul 2>&1
    if errorlevel 1 (
        echo [SKIP] CTest not found. Skipping tests.
    ) else (
        echo Running tests to verify build...
        ctest --test-dir build --output-on-failure -C Debug
        if errorlevel 1 (
            echo.
            echo [WARNING] Some tests failed. Check the output above.
        ) else (
            echo [OK] All tests passed.
        )
    )
) else (
    echo.
    echo [SKIP] Tests skipped (--skip-test).
)

rem ---------------------------------------------
rem [9] Summary
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
echo   To run tests:
echo     ctest --test-dir build -C Debug
echo.
echo ============================================================
echo.

if "!UNATTENDED!"=="0" pause
endlocal
exit /b 0
