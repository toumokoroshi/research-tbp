@echo off
setlocal enabledelayedexpansion

rem ============================================
rem TBP Project - Dependency Checker
rem ============================================
rem This script checks for required dependencies
rem and reports their status.
rem ============================================

echo ============================================
echo   TBP Project - Dependency Check
echo ============================================
echo.

set "MISSING_REQUIRED=0"
set "HAS_ONEAPI=0"
set "HAS_NINJA=0"
set "HAS_MINGW=0"
set "HAS_BOOST=0"
set "HAS_VCPKG=0"

rem ---------------------------------------------
rem [1] Git Check (Required)
rem ---------------------------------------------
echo [1/7] Checking Git...
where git >nul 2>&1
if errorlevel 1 (
    echo     [MISSING] Git is not installed or not in PATH.
    echo              Please install Git from: https://git-scm.com/
    set "MISSING_REQUIRED=1"
) else (
    for /f "tokens=3" %%v in ('git --version') do set "GIT_VERSION=%%v"
    echo     [OK] Git !GIT_VERSION!
)

rem ---------------------------------------------
rem [2] CMake Check (Required)
rem ---------------------------------------------
echo [2/7] Checking CMake...
where cmake >nul 2>&1
if errorlevel 1 (
    echo     [MISSING] CMake is not installed or not in PATH.
    echo              Please install CMake from: https://cmake.org/
    set "MISSING_REQUIRED=1"
) else (
    for /f "tokens=3" %%v in ('cmake --version ^| findstr /i "cmake version"') do set "CMAKE_VERSION=%%v"
    echo     [OK] CMake !CMAKE_VERSION!
)

rem ---------------------------------------------
rem [3] Ninja Check (Optional - for Intel oneAPI)
rem ---------------------------------------------
echo [3/7] Checking Ninja...
where ninja >nul 2>&1
if errorlevel 1 (
    echo     [NOT FOUND] Ninja is not installed or not in PATH.
    echo                 Ninja is optional but recommended for faster builds.
    echo                 Install via: winget install Ninja-build.Ninja
    set "HAS_NINJA=0"
) else (
    for /f "tokens=1" %%v in ('ninja --version') do set "NINJA_VERSION=%%v"
    echo     [OK] Ninja !NINJA_VERSION!
    set "HAS_NINJA=1"
)

rem ---------------------------------------------
rem [4] Intel oneAPI Check (Optional)
rem ---------------------------------------------
echo [4/7] Checking Intel oneAPI...
if exist "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" (
    echo     [OK] Intel oneAPI found at: C:\Program Files ^(x86^)\Intel\oneAPI
    set "HAS_ONEAPI=1"
    
    rem Check if icx-cl is available after setvars
    where icx-cl >nul 2>&1
    if errorlevel 1 (
        echo          Note: Run setvars.bat to add icx-cl to PATH
    ) else (
        echo          icx-cl is available in PATH
    )
) else (
    echo     [NOT FOUND] Intel oneAPI is not installed.
    echo                 oneAPI is optional but provides better performance.
    echo                 Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html
    set "HAS_ONEAPI=0"
)

rem ---------------------------------------------
rem [5] MinGW-gcc Check (Fallback compiler)
rem ---------------------------------------------
echo [5/7] Checking MinGW-gcc...
where g++ >nul 2>&1
if errorlevel 1 (
    echo     [NOT FOUND] MinGW g++ is not installed or not in PATH.
    echo                 MinGW is used as fallback when oneAPI is not available.
    echo                 Install via: winget install -e --id mingw.mingw-w64-ucrt-x86_64
    set "HAS_MINGW=0"
) else (
    for /f "tokens=*" %%v in ('g++ --version ^| findstr /i "g++"') do set "GCC_VERSION=%%v"
    echo     [OK] !GCC_VERSION!
    set "HAS_MINGW=1"
)

rem ---------------------------------------------
rem [6] vcpkg Check
rem ---------------------------------------------
echo [6/7] Checking vcpkg...
set "VCPKG_ROOT=%~dp0..\external\vcpkg"
if exist "!VCPKG_ROOT!\vcpkg.exe" (
    echo     [OK] vcpkg found at: !VCPKG_ROOT!
    set "HAS_VCPKG=1"
) else (
    echo     [NOT FOUND] vcpkg is not installed in external\vcpkg
    echo                 Run setup.bat to install vcpkg automatically.
    set "HAS_VCPKG=0"
)

rem ---------------------------------------------
rem [7] Boost Check
rem ---------------------------------------------
echo [7/7] Checking Boost...
set "HAS_BOOST=0"

rem First check BOOST_ROOT environment variable
if defined BOOST_ROOT (
    if not "!BOOST_ROOT!"=="" (
        if exist "!BOOST_ROOT!\boost" (
            echo     [OK] Boost found at: !BOOST_ROOT! ^(from environment variable^)
            set "HAS_BOOST=1"
        ) else if exist "!BOOST_ROOT!\include\boost" (
            echo     [OK] Boost found at: !BOOST_ROOT! ^(from environment variable^)
            set "HAS_BOOST=1"
        ) else (
            echo     [WARNING] BOOST_ROOT is set but Boost headers not found at: !BOOST_ROOT!
        )
    )
)

rem If not found via BOOST_ROOT, check vcpkg
if "!HAS_BOOST!"=="0" (
    set "VCPKG_BOOST=!VCPKG_ROOT!\installed\x64-windows\include\boost"
    if exist "!VCPKG_BOOST!" (
        echo     [OK] Boost found via vcpkg
        set "HAS_BOOST=1"
    ) else (
        echo     [NOT FOUND] Boost is not installed.
        echo                 Run setup.bat to install Boost via vcpkg.
    )
)

rem ---------------------------------------------
rem Summary
rem ---------------------------------------------
echo.
echo ============================================
echo   Dependency Check Summary
echo ============================================

if "!MISSING_REQUIRED!"=="1" (
    echo [ERROR] Required dependencies are missing!
    echo         Please install Git and CMake before continuing.
    exit /b 1
)

echo Required:
echo   - Git:   OK
echo   - CMake: OK
echo.
echo Optional:
if "!HAS_NINJA!"=="1" (
    echo   - Ninja:      OK
) else (
    echo   - Ninja:      NOT FOUND
)
if "!HAS_ONEAPI!"=="1" (
    echo   - Intel oneAPI: OK
) else (
    echo   - Intel oneAPI: NOT FOUND
)
if "!HAS_MINGW!"=="1" (
    echo   - MinGW-gcc:  OK
) else (
    echo   - MinGW-gcc:  NOT FOUND
)
echo.
echo Libraries:
if "!HAS_VCPKG!"=="1" (
    echo   - vcpkg:      OK
) else (
    echo   - vcpkg:      NOT FOUND ^(will be installed by setup.bat^)
)
if "!HAS_BOOST!"=="1" (
    echo   - Boost:      OK
) else (
    echo   - Boost:      NOT FOUND ^(will be installed by setup.bat^)
)

echo.

rem Determine available build configuration
if "!HAS_ONEAPI!"=="1" if "!HAS_NINJA!"=="1" (
    echo [RECOMMENDED] Intel oneAPI + Ninja build is available.
)
if "!HAS_MINGW!"=="1" (
    echo [AVAILABLE] MinGW-gcc build is available.
)
if "!HAS_ONEAPI!"=="0" if "!HAS_MINGW!"=="0" (
    echo [WARNING] No compiler found! Please install Intel oneAPI or MinGW-gcc.
    exit /b 1
)

echo ============================================
echo.

rem Export variables for parent script
endlocal & (
    set "HAS_ONEAPI=%HAS_ONEAPI%"
    set "HAS_NINJA=%HAS_NINJA%"
    set "HAS_MINGW=%HAS_MINGW%"
    set "HAS_BOOST=%HAS_BOOST%"
    set "HAS_VCPKG=%HAS_VCPKG%"
)

exit /b 0
