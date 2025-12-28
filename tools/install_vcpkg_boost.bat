@echo off
setlocal enabledelayedexpansion

rem ============================================
rem TBP Project - vcpkg & Boost Installer
rem ============================================
rem This script installs vcpkg and Boost to
rem the external/ directory and sets up
rem environment variables.
rem ============================================

echo ============================================
echo   TBP Project - vcpkg ^& Boost Installer
echo ============================================
echo.

rem Get the project root directory
set "PROJECT_ROOT=%~dp0.."
set "EXTERNAL_DIR=!PROJECT_ROOT!\external"
set "VCPKG_ROOT=!EXTERNAL_DIR!\vcpkg"

rem ---------------------------------------------
rem [1] Create external directory if needed
rem ---------------------------------------------
if not exist "!EXTERNAL_DIR!" (
    echo Creating external directory...
    mkdir "!EXTERNAL_DIR!"
)

rem ---------------------------------------------
rem [2] Install vcpkg if not present
rem ---------------------------------------------
echo.
echo [1/4] Checking vcpkg installation...

if exist "!VCPKG_ROOT!\vcpkg.exe" (
    echo       vcpkg is already installed.
) else (
    echo       Installing vcpkg to: !VCPKG_ROOT!
    echo       This may take a few minutes...
    echo.
    
    rem Clone vcpkg repository
    git clone https://github.com/microsoft/vcpkg.git "!VCPKG_ROOT!"
    if errorlevel 1 (
        echo [ERROR] Failed to clone vcpkg repository.
        exit /b 1
    )
    
    rem Bootstrap vcpkg
    echo.
    echo       Running vcpkg bootstrap...
    call "!VCPKG_ROOT!\bootstrap-vcpkg.bat" -disableMetrics
    if errorlevel 1 (
        echo [ERROR] Failed to bootstrap vcpkg.
        exit /b 1
    )
    
    echo       [OK] vcpkg installed successfully.
)

rem ---------------------------------------------
rem [3] Install Boost (odeint) via vcpkg
rem ---------------------------------------------
echo.
echo [2/4] Checking Boost installation...

set "BOOST_ODEINT_INCLUDE=!VCPKG_ROOT!\installed\x64-windows\include\boost\numeric\odeint.hpp"
if exist "!BOOST_ODEINT_INCLUDE!" (
    echo       Boost.odeint is already installed via vcpkg.
) else (
    echo       Installing Boost.odeint via vcpkg...
    echo       ^(Installing only required components for faster setup^)
    echo       This may take 2-5 minutes depending on your system.
    echo.
    
    "!VCPKG_ROOT!\vcpkg.exe" install boost-odeint:x64-windows
    if errorlevel 1 (
        echo [ERROR] Failed to install Boost.odeint.
        exit /b 1
    )
    
    echo       [OK] Boost.odeint installed successfully.
)

rem ---------------------------------------------
rem [4] Set BOOST_ROOT environment variable
rem ---------------------------------------------
echo.
echo [3/4] Setting up BOOST_ROOT environment variable...

set "NEW_BOOST_ROOT=!VCPKG_ROOT!\installed\x64-windows"

rem Check if BOOST_ROOT is already set
if defined BOOST_ROOT (
    echo       Current BOOST_ROOT: !BOOST_ROOT!
    echo       New BOOST_ROOT:     !NEW_BOOST_ROOT!
    echo.
    
    rem Compare paths (case-insensitive)
    if /i "!BOOST_ROOT!"=="!NEW_BOOST_ROOT!" (
        echo       BOOST_ROOT is already correctly configured.
        goto :skip_setx
    )
    
    echo       BOOST_ROOT is set to a different path.
    set /p "OVERWRITE=       Do you want to update BOOST_ROOT? (Y/N): "
    if /i "!OVERWRITE!" neq "Y" (
        echo       Keeping existing BOOST_ROOT.
        goto :skip_setx
    )
)

echo       Registering BOOST_ROOT as user environment variable...
echo       Path: !NEW_BOOST_ROOT!

rem Set for current session
set "BOOST_ROOT=!NEW_BOOST_ROOT!"

rem Set permanently using setx (user environment variable)
setx BOOST_ROOT "!NEW_BOOST_ROOT!"
if errorlevel 1 (
    echo [WARNING] Failed to set BOOST_ROOT permanently.
    echo          Please set BOOST_ROOT manually to: !NEW_BOOST_ROOT!
) else (
    echo       [OK] BOOST_ROOT registered as user environment variable.
    echo.
    echo       NOTE: The environment variable has been set permanently.
    echo             New terminal windows will have BOOST_ROOT available.
)

:skip_setx

rem ---------------------------------------------
rem [5] Summary
rem ---------------------------------------------
echo.
echo [4/4] Installation Summary
echo ============================================
echo   vcpkg:     !VCPKG_ROOT!
echo   Boost:     !NEW_BOOST_ROOT!\include\boost
echo   BOOST_ROOT: !NEW_BOOST_ROOT!
echo ============================================
echo.
echo Installation complete!
echo.

endlocal & (
    set "BOOST_ROOT=%NEW_BOOST_ROOT%"
    set "VCPKG_ROOT=%VCPKG_ROOT%"
)

exit /b 0
