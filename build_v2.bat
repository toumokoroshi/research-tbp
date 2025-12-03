@echo off
set PATH=%PATH%;C:\Program Files (x86)\Intel\oneAPI\compiler\2025.2\bin
taskkill /F /IM trajectory_SALI.exe 2>nul
taskkill /F /IM gnuplot.exe 2>nul
if exist build\bin\Debug\trajectory_SALI.exe del build\bin\Debug\trajectory_SALI.exe
cmake --build build --config Debug --target trajectory_SALI
if %errorlevel% neq 0 (
    echo Build failed
    exit /b 1
)
echo Build success. Running app with input '4'...
build\bin\Debug\trajectory_SALI.exe
echo App finished.
