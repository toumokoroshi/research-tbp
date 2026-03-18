@echo off
setlocal enabledelayedexpansion

rem ============================================
rem TBP Project - Dependency Checker
rem ============================================
rem This script checks for required dependencies
rem and reports their status.
rem 
rem Usage: check_dependencies.bat [UNATTENDED] [WITH_ONEAPI]
rem   UNATTENDED=1 : Auto-install without prompts
rem   WITH_ONEAPI=1 : Install Intel oneAPI (default)
rem   WITH_ONEAPI=0 : Skip Intel oneAPI installation
rem ============================================

rem Parse arguments
set "UNATTENDED=0"
set "WITH_ONEAPI=1"
if "%~1"=="1" set "UNATTENDED=1"
if "%~2"=="1" set "WITH_ONEAPI=1"
if "%~2"=="0" set "WITH_ONEAPI=0"

echo ============================================
echo   TBP Project - Dependency Check
if "!UNATTENDED!"=="1" echo   (Unattended Mode)
if "!WITH_ONEAPI!"=="1" echo   (With Intel oneAPI(ちょっと時間かかるかもよどんまい。（10～30分）))
echo ============================================
echo.

set "MISSING_REQUIRED=0"
set "HAS_ONEAPI=0"
set "HAS_NINJA=0"
set "HAS_MINGW=0"
set "HAS_WINSDK=0"
set "HAS_BOOST=0"
set "HAS_VCPKG=0"
set "HAS_WINGET=0"

rem ---------------------------------------------
rem [0/9] winget Check (Required for auto-install)
rem ---------------------------------------------
echo [0/9] Checking winget...
call :ensure_winget
if "!HAS_WINGET!"=="1" (
    echo     [OK] winget is available.
) else (
    echo     [WARNING] winget is not available.
    echo               Auto-install is disabled; install App Installer manually if needed.
)

rem ---------------------------------------------
rem [1/9] Git Check (Required)
rem ---------------------------------------------
echo [1/9] Checking Git...
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
rem [2/9] CMake Check (Required)
rem ---------------------------------------------
echo [2/9] Checking CMake...
where cmake >nul 2>&1
if errorlevel 1 (
    echo     [MISSING] CMake is not installed or not in PATH.
    if "!UNATTENDED!"=="1" (
        set "INSTALL_CMAKE=Y"
    ) else (
        set /p "INSTALL_CMAKE=    Install CMake now? [Y/n]: "
        if /i "!INSTALL_CMAKE!"=="" set "INSTALL_CMAKE=Y"
    )
    if /i "!INSTALL_CMAKE!"=="Y" (
        if "!HAS_WINGET!"=="0" (
            echo     [ERROR] Cannot auto-install CMake because winget is unavailable.
            set "MISSING_REQUIRED=1"
        ) else (
            echo     Installing CMake...
            winget install Kitware.CMake --silent --accept-source-agreements --accept-package-agreements
            if errorlevel 1 (
                echo     [ERROR] CMake installation failed.
                set "MISSING_REQUIRED=1"
            ) else (
                echo     [OK] CMake installed.
                rem Add CMake to session PATH immediately
                set "CMAKE_PATH="
                if exist "C:\Program Files\CMake\bin\cmake.exe" (
                    set "CMAKE_PATH=C:\Program Files\CMake\bin"
                ) else if exist "%LOCALAPPDATA%\Programs\CMake\bin\cmake.exe" (
                    set "CMAKE_PATH=%LOCALAPPDATA%\Programs\CMake\bin"
                )
                if defined CMAKE_PATH (
                    echo          Adding to session PATH: !CMAKE_PATH!
                    set "PATH=!CMAKE_PATH!;!PATH!"
                    rem Add to permanent user PATH
                    for /f "tokens=2*" %%A in ('reg query "HKCU\Environment" /v Path 2^>nul') do set "USER_PATH=%%B"
                    if not defined USER_PATH set "USER_PATH="
                    echo !USER_PATH! | findstr /i /c:"CMake" >nul 2>&1
                    if errorlevel 1 (
                        setx PATH "!CMAKE_PATH!;!USER_PATH!"
                        echo          Added to permanent user PATH.
                    )
                ) else (
                    echo          [WARNING] CMake path not found. Restart terminal to use.
                    set "MISSING_REQUIRED=1"
                )
            )
        )
    ) else (
        set "MISSING_REQUIRED=1"
    )
) else (
    for /f "tokens=3" %%v in ('cmake --version ^| findstr /i "cmake version"') do set "CMAKE_VERSION=%%v"
    echo     [OK] CMake !CMAKE_VERSION!
)

rem ---------------------------------------------
rem [3/9] Ninja Check (Required for vs-intel)
rem ---------------------------------------------
echo [3/9] Checking Ninja...
set "NINJA_FOUND=0"

rem First check if ninja is already in PATH
where ninja >nul 2>&1
if not errorlevel 1 (
    set "NINJA_FOUND=1"
)

rem If not found, check winget portable packages directory
if "!NINJA_FOUND!"=="0" (
    set "WINGET_PACKAGES=%LOCALAPPDATA%\Microsoft\WinGet\Packages"
    for /d %%D in ("!WINGET_PACKAGES!\Ninja-build.Ninja_*") do (
        if exist "%%D\ninja.exe" (
            echo     [INFO] Found Ninja in winget packages: %%D
            set "PATH=%%D;!PATH!"
            set "NINJA_FOUND=1"
        )
    )
)

rem If not found, check VS bundled Ninja
if "!NINJA_FOUND!"=="0" (
    for /d %%V in ("C:\Program Files\Microsoft Visual Studio\*") do (
        for /d %%E in ("%%V\*") do (
            set "VS_NINJA=%%E\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja"
            if exist "!VS_NINJA!\ninja.exe" (
                echo     [INFO] Found Ninja bundled with VS: !VS_NINJA!
                set "PATH=!VS_NINJA!;!PATH!"
                set "NINJA_FOUND=1"
            )
        )
    )
)

if "!NINJA_FOUND!"=="0" (
    echo     [NOT FOUND] Ninja is not installed or not in PATH.
    echo                 Ninja is required for Intel oneAPI builds.
    if "!UNATTENDED!"=="1" (
        set "INSTALL_NINJA=Y"
    ) else (
        set /p "INSTALL_NINJA=    Install Ninja now? [Y/n]: "
        if /i "!INSTALL_NINJA!"=="" set "INSTALL_NINJA=Y"
    )
    if /i "!INSTALL_NINJA!"=="Y" (
        if "!HAS_WINGET!"=="0" (
            echo     [ERROR] Cannot auto-install Ninja because winget is unavailable.
            set "HAS_NINJA=0"
        ) else (
            echo     Installing Ninja...
            winget install Ninja-build.Ninja --silent --accept-source-agreements --accept-package-agreements
            if errorlevel 1 (
                echo     [ERROR] Ninja installation failed.
                set "HAS_NINJA=0"
            ) else (
                echo     [OK] Ninja installed.
                rem Add winget package path immediately
                set "NINJA_PATH="
                for /d %%D in ("%LOCALAPPDATA%\Microsoft\WinGet\Packages\Ninja-build.Ninja_*") do (
                    if exist "%%D\ninja.exe" (
                        set "NINJA_PATH=%%D"
                        set "PATH=%%D;!PATH!"
                    )
                )
                if defined NINJA_PATH (
                    echo          Adding to session PATH: !NINJA_PATH!
                    rem Add to permanent user PATH
                    for /f "tokens=2*" %%A in ('reg query "HKCU\Environment" /v Path 2^>nul') do set "USER_PATH=%%B"
                    if not defined USER_PATH set "USER_PATH="
                    echo !USER_PATH! | findstr /i /c:"Ninja" >nul 2>&1
                    if errorlevel 1 (
                        setx PATH "!NINJA_PATH!;!USER_PATH!"
                        echo          Added to permanent user PATH.
                    )
                )
                set "HAS_NINJA=1"
            )
        )
    ) else (
        set "HAS_NINJA=0"
    )
) else (
    for /f "tokens=1" %%v in ('ninja --version') do set "NINJA_VERSION=%%v"
    echo     [OK] Ninja !NINJA_VERSION!
    set "HAS_NINJA=1"
)

rem ---------------------------------------------
rem [4/9] Intel oneAPI Check (Recommended)
rem ---------------------------------------------
echo [4/9] Checking Intel oneAPI...
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
    echo                 Intel oneAPI is recommended for best performance.
    if "!UNATTENDED!"=="1" if "!WITH_ONEAPI!"=="0" (
        rem Explicitly disabled by --without-oneapi
        echo     [SKIP] Skipping Intel oneAPI installation in unattended mode.
        echo            Use default behavior or --with-oneapi to install.
        set "HAS_ONEAPI=0"
    ) else (
        rem Interactive mode OR unattended with oneAPI enabled
        if "!WITH_ONEAPI!"=="1" (
            set "INSTALL_ONEAPI=Y"
            echo     [AUTO] Installing Intel oneAPI (enabled by default)
        ) else (
            set /p "INSTALL_ONEAPI=    Install Intel oneAPI Base Toolkit now? [Y/n]: "
            if /i "!INSTALL_ONEAPI!"=="" set "INSTALL_ONEAPI=Y"
        )
        if /i "!INSTALL_ONEAPI!"=="Y" (
            if "!HAS_WINGET!"=="0" (
                echo     [ERROR] Cannot auto-install Intel oneAPI because winget is unavailable.
                set "HAS_ONEAPI=0"
            ) else (
                echo.
                echo     Installing Intel oneAPI Base Toolkit...
                echo     This may take 10-30 minutes...
                winget install Intel.oneAPI.BaseToolkit --silent --accept-source-agreements --accept-package-agreements
                if errorlevel 1 (
                    echo     [ERROR] Intel oneAPI installation failed.
                    echo              Manual download: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html
                    set "HAS_ONEAPI=0"
                ) else (
                    echo     [OK] Intel oneAPI installed.
                    set "HAS_ONEAPI=1"
                )
            )
        ) else (
            set "HAS_ONEAPI=0"
        )
    )
)

rem ---------------------------------------------
rem [5/9] MinGW-gcc Check (Fallback compiler)
rem ---------------------------------------------
echo [5/9] Checking MinGW-gcc...
where g++ >nul 2>&1
if errorlevel 1 (
    echo     [NOT FOUND] MinGW g++ is not installed or not in PATH.
    echo                 MinGW is used as fallback when oneAPI is not available.
    if "!UNATTENDED!"=="1" (
        set "INSTALL_MINGW=Y"
    ) else (
        set /p "INSTALL_MINGW=    Install MinGW-gcc now? [Y/n]: "
        if /i "!INSTALL_MINGW!"=="" set "INSTALL_MINGW=Y"
    )
    if /i "!INSTALL_MINGW!"=="Y" (
        if "!HAS_WINGET!"=="0" (
            echo     [ERROR] Cannot auto-install MinGW because winget is unavailable.
            set "HAS_MINGW=0"
        ) else (
            echo     Installing MinGW-gcc via winget...
            winget install mingw.mingw-w64-ucrt-x86_64 --silent --accept-source-agreements --accept-package-agreements
            if errorlevel 1 (
                echo     [ERROR] MinGW installation failed.
                set "HAS_MINGW=0"
            ) else (
                echo     [OK] MinGW-gcc installed.
                rem Try to find and add to PATH immediately
                set "MINGW_PATH="
                if exist "C:\mingw64\bin\g++.exe" (
                    set "MINGW_PATH=C:\mingw64\bin"
                ) else if exist "C:\Program Files\mingw64\bin\g++.exe" (
                    set "MINGW_PATH=C:\Program Files\mingw64\bin"
                ) else (
                    rem Check winget package location
                    for /d %%D in ("%LOCALAPPDATA%\Microsoft\WinGet\Packages\mingw*") do (
                        if exist "%%D\mingw64\bin\g++.exe" (
                            set "MINGW_PATH=%%D\mingw64\bin"
                        )
                    )
                )
                if defined MINGW_PATH (
                    echo          Adding to session PATH: !MINGW_PATH!
                    set "PATH=!MINGW_PATH!;!PATH!"
                    rem Add to permanent user PATH
                    for /f "tokens=2*" %%A in ('reg query "HKCU\Environment" /v Path 2^>nul') do set "USER_PATH=%%B"
                    if not defined USER_PATH set "USER_PATH="
                    echo !USER_PATH! | findstr /i /c:"mingw" >nul 2>&1
                    if errorlevel 1 (
                        setx PATH "!MINGW_PATH!;!USER_PATH!"
                        echo          Added to permanent user PATH.
                    )
                    set "HAS_MINGW=1"
                ) else (
                    echo          [WARNING] MinGW path not found. Restart terminal to use.
                    set "HAS_MINGW=0"
                )
            )
        )
    ) else (
        set "HAS_MINGW=0"
    )
) else (
    for /f "tokens=*" %%v in ('g++ --version ^| findstr /i "g++"') do set "GCC_VERSION=%%v"
    echo     [OK] !GCC_VERSION!
    set "HAS_MINGW=1"
)

rem ---------------------------------------------
rem [6/9] Windows SDK Check (Needed only when vcpkg Boost install is required)
rem ---------------------------------------------
echo [6/9] Checking Windows SDK (rc.exe)...
set "WINSDK_RC_EXE="
set "EXISTING_VS_FOUND=0"

rem First, check if any existing VS installation has vcvars64.bat
rem This covers VS 2019, 2022, 2026, etc.
for /d %%V in ("C:\Program Files\Microsoft Visual Studio\*") do (
    for /d %%E in ("%%V\*") do (
        if exist "%%E\VC\Auxiliary\Build\vcvars64.bat" (
            echo     [INFO] Found existing VS installation: %%V\%%~nxE
            set "EXISTING_VS_FOUND=1"
        )
    )
)
rem Also check Program Files (x86) for older installations
for /d %%V in ("C:\Program Files (x86)\Microsoft Visual Studio\*") do (
    for /d %%E in ("%%V\*") do (
        if exist "%%E\VC\Auxiliary\Build\vcvars64.bat" (
            echo     [INFO] Found existing VS installation: %%V\%%~nxE
            set "EXISTING_VS_FOUND=1"
        )
    )
)

rem Search for rc.exe in Windows Kits
for /d %%D in ("C:\Program Files (x86)\Windows Kits\10\bin\10.*") do (
    if exist "%%D\x64\rc.exe" (
        set "WINSDK_RC_EXE=%%D\x64\rc.exe"
    )
)

if defined WINSDK_RC_EXE (
    echo     [OK] Windows SDK found: !WINSDK_RC_EXE!
    set "HAS_WINSDK=1"
) else (
    rem Check if we have existing VS but no Windows SDK
    if "!EXISTING_VS_FOUND!"=="1" (
        echo     [WARNING] VS Build Tools found but Windows SDK is missing.
        echo              Please install Windows SDK via Visual Studio Installer.
        set "HAS_WINSDK=0"
    ) else (
        echo     [MISSING] Windows SDK is not installed or rc.exe not found.
        echo              vcpkg requires Windows SDK to build packages.
        echo.
        if "!UNATTENDED!"=="1" (
            set "INSTALL_VSBT=Y"
        ) else (
            set /p "INSTALL_VSBT=    Install Visual Studio Build Tools now? [Y/n]: "
            if /i "!INSTALL_VSBT!"=="" set "INSTALL_VSBT=Y"
        )
        if /i "!INSTALL_VSBT!"=="Y" (
            if "!HAS_WINGET!"=="0" (
                echo     [WARNING] Cannot auto-install Visual Studio Build Tools because winget is unavailable.
                set "HAS_WINSDK=0"
            ) else (
                echo.
                echo     Installing Visual Studio Build Tools with required components...
                echo     This may take several minutes...
                winget install Microsoft.VisualStudio.2022.BuildTools --accept-source-agreements --accept-package-agreements --override "--quiet --wait --add Microsoft.VisualStudio.Workload.VCTools --add Microsoft.VisualStudio.Component.VC.Tools.x86.x64 --add Microsoft.VisualStudio.Component.Windows11SDK.22621"
                if errorlevel 1 (
                    echo     [ERROR] Installation failed. Please install manually.
                    set "HAS_WINSDK=0"
                ) else (
                    echo     [OK] Visual Studio Build Tools installed successfully.
                    rem Try to find and set up vcvars64.bat path for this session
                    set "VCVARS_FOUND=0"
                    for /d %%V in ("C:\Program Files\Microsoft Visual Studio\*") do (
                        for /d %%E in ("%%V\*") do (
                            if exist "%%E\VC\Auxiliary\Build\vcvars64.bat" (
                                set "VCVARS64_PATH=%%E\VC\Auxiliary\Build\vcvars64.bat"
                                set "VCVARS_FOUND=1"
                            )
                        )
                    )
                    if "!VCVARS_FOUND!"=="1" (
                        echo          Found vcvars64.bat: !VCVARS64_PATH!
                        rem Re-check for Windows SDK after installation
                        for /d %%D in ("C:\Program Files (x86)\Windows Kits\10\bin\10.*") do (
                            if exist "%%D\x64\rc.exe" (
                                set "WINSDK_RC_EXE=%%D\x64\rc.exe"
                                set "HAS_WINSDK=1"
                            )
                        )
                        if "!HAS_WINSDK!"=="1" (
                            echo          Windows SDK detected: !WINSDK_RC_EXE!
                        ) else (
                            echo          [WARNING] Windows SDK not detected yet. Restart may be required.
                        )
                    ) else (
                        echo          [WARNING] vcvars64.bat not found. Restart may be required.
                        set "HAS_WINSDK=0"
                    )
                    if "!HAS_WINSDK!"=="0" (
                        for /d %%D in ("C:\Program Files (x86)\Windows Kits\10\bin\10.*") do (
                            if exist "%%D\x64\rc.exe" (
                                set "WINSDK_RC_EXE=%%D\x64\rc.exe"
                                set "HAS_WINSDK=1"
                            )
                        )
                    )
                )
            )
        ) else (
            echo.
            echo              Manual install via Visual Studio Installer ^(Individual components^):
            echo                - Windows 11 SDK ^(10.0.22621^)
            echo                - MSVC v143 - VS 2022 C++ x64/x86 build tools
            set "HAS_WINSDK=0"
        )
    )
)

rem ---------------------------------------------
rem [7/9] vcpkg Check
rem ---------------------------------------------
echo [7/9] Checking vcpkg...
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
rem [8/9] Boost Check
rem ---------------------------------------------
echo [8/9] Checking Boost...
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
rem [9/9] Final requirement checks
rem ---------------------------------------------
if "!HAS_BOOST!"=="0" if "!HAS_WINSDK!"=="0" (
    echo.
    echo [ERROR] Boost is missing and Windows SDK is not available.
    echo         Windows SDK is required for vcpkg to build Boost on Windows.
    set "MISSING_REQUIRED=1"
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
    echo         Please install missing dependencies before continuing.
    echo.
    if "!HAS_BOOST!"=="0" if "!HAS_WINSDK!"=="0" (
        echo         [^^!] Windows SDK is required for vcpkg to build packages.
        echo             Install via Visual Studio Installer.
    )
    exit /b 1
)

echo Required:
echo   - Git:       OK
echo   - CMake:     OK
if "!HAS_BOOST!"=="0" (
if "!HAS_WINSDK!"=="1" (
    echo   - Windows SDK: OK
) else (
    echo   - Windows SDK: MISSING
)
) else (
    echo   - Windows SDK: SKIP ^(Boost already available^)
)
echo.
echo Optional:
if "!HAS_WINGET!"=="1" (
    echo   - winget:     OK
) else (
    echo   - winget:     NOT AVAILABLE ^(auto-install disabled^)
)
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
    set "HAS_WINGET=%HAS_WINGET%"
    set "HAS_ONEAPI=%HAS_ONEAPI%"
    set "HAS_NINJA=%HAS_NINJA%"
    set "HAS_MINGW=%HAS_MINGW%"
    set "HAS_WINSDK=%HAS_WINSDK%"
    set "HAS_BOOST=%HAS_BOOST%"
    set "HAS_VCPKG=%HAS_VCPKG%"
)

exit /b 0

:ensure_winget
set "HAS_WINGET=0"
where winget >nul 2>&1
if errorlevel 1 goto :winget_missing
winget --info >nul 2>&1
if not errorlevel 1 (
    set "HAS_WINGET=1"
    exit /b 0
)

:winget_missing
echo     [INFO] Attempting to repair/install App Installer (winget)...
powershell -NoProfile -ExecutionPolicy Bypass -Command "try { Add-AppxPackage -RegisterByFamilyName -MainPackage 'Microsoft.DesktopAppInstaller_8wekyb3d8bbwe' -ErrorAction Stop; exit 0 } catch { exit 1 }" >nul 2>&1
if errorlevel 1 (
    powershell -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference='Stop'; $out=Join-Path $env:TEMP 'Microsoft.DesktopAppInstaller.msixbundle'; Invoke-WebRequest 'https://aka.ms/getwinget' -OutFile $out; Add-AppxPackage -Path $out" >nul 2>&1
)

where winget >nul 2>&1
if errorlevel 1 exit /b 1
winget --info >nul 2>&1
if errorlevel 1 exit /b 1
set "HAS_WINGET=1"
exit /b 0
