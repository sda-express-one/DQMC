@echo off
setlocal enabledelayedexpansion

set COMPILER=g++
set LINKER=g++
set output=DQMC.exe
set CFLAGS=-c -fopenmp -Wall -Wextra -Werror -O3
set LFLAGS=

REM Create obj directory if it doesn't exist
if not exist build mkdir build

cd src

for %%f in (*.cpp) do (
    echo   Compiling %%f...
    %COMPILER% %CFLAGS% "%%f" -o "%%~nf.o"
    if !errorlevel! neq 0 (
        echo Failed to compile %%f
        pause
        exit /b 1
    )
)

cd utils 
for %%f in (*.cpp) do (
    echo   Compiling %%f...
    %COMPILER% %CFLAGS% "%%f" -o "..\%%~nf.o"
    if !errorlevel! neq 0 (
        echo Failed to compile %%f
        pause
        exit /b 1
    )
)

cd ..

echo.
echo Step 2: Linking object files...
%LINKER% *.o %LFLAGS% -o %OUTPUT%

if %errorlevel% equ 0 (
    echo.
    echo ============================
    echo Build successful!
    echo Output: %OUTPUT%
    echo ============================
) else (
    echo.
    echo ============================
    echo Linking failed!
    echo ============================
)

move %OUTPUT% ..\build\%OUTPUT%

del -f "*.o"

pause