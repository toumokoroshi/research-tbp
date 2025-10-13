rem #############################################
echo Setting up Intel oneAPI environment...
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
echo.
echo oneAPI environment has been set for this terminal session.


echo.
echo Checking compiler version...
icx --version