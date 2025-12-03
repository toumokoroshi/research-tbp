@echo off
set PATH=%PATH%;C:\Program Files (x86)\Intel\oneAPI\compiler\2025.2\bin
echo Deleting old exe... > debug_log.txt
if exist build\bin\Debug\trajectory_SALI_v2.exe del build\bin\Debug\trajectory_SALI_v2.exe
echo Starting build... >> debug_log.txt
cmake --build build --config Debug --target trajectory_SALI >> debug_log.txt 2>&1
if %errorlevel% neq 0 (
    echo Build failed >> debug_log.txt
    exit /b 1
)
echo Build success. Running app... >> debug_log.txt
build\bin\Debug\trajectory_SALI_v2.exe >> debug_log.txt 2>&1
echo App finished. >> debug_log.txt
