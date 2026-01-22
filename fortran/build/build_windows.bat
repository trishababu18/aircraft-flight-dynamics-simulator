@echo off
REM Windows build script for F-16 Flight Simulator

echo.
echo ================================================================================
echo Building F-16 Flight Simulator
echo ================================================================================
echo.

REM Create ws2_32.def if it doesn't exist
@REM if not exist ws2_32.def (
@REM     echo Creating ws2_32.def...
@REM     (
@REM         echo LIBRARY ws2_32.dll
@REM         echo EXPORTS
@REM         echo WSAStartup
@REM         echo WSACleanup
@REM         echo socket
@REM         echo bind
@REM         echo sendto
@REM         echo recvfrom
@REM         echo closesocket
@REM         echo htons
@REM         echo inet_addr
@REM         echo ioctlsocket
@REM     ) > ws2_32.def
@REM     echo Created ws2_32.def
@REM     echo.
@REM )

@REM REM Check if we need to create ws2_32 import library
@REM if not exist libws2_32.a (
@REM     echo Creating ws2_32 import library from ws2_32.def...
@REM     dlltool -d ws2_32.def -l libws2_32.a -D ws2_32.dll
@REM     if errorlevel 1 (
@REM         echo ERROR: Failed to create import library
@REM         echo Make sure dlltool is in your PATH
@REM         goto :error
@REM     )
@REM     echo Successfully created libws2_32.a
@REM     echo.
@REM )

echo Cleaning up old build files...
del /Q main.exe *.o *.mod 2>NUL

echo.
echo Compiling Fortran modules...
echo   - babu.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c babu.f90
if errorlevel 1 goto :error

echo   - linalg_mod.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c linalg_mod.f90
if errorlevel 1 goto :error

echo   - json.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c json.f90
if errorlevel 1 goto :error

echo   - jsonx.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c jsonx.f90
if errorlevel 1 goto :error

echo   - database_m.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c database_m.f90
if errorlevel 1 goto :error

echo   - udp_windows_m.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c udp_windows_m.f90
if errorlevel 1 goto :error

echo   - connection_m.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c connection_m.f90
if errorlevel 1 goto :error

echo   - sim.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c sim.f90
if errorlevel 1 goto :error

echo   - main.f90
gfortran -fdefault-real-8 -ffree-line-length-none -c main.f90
if errorlevel 1 goto :error

echo.
echo Linking executable...
@REM gfortran -fdefault-real-8 -ffree-line-length-none babu.o linalg_mod.o json.o jsonx.o database_m.o udp_windows_m.o connection_m.o sim.o main.o -o main.exe -L. -lws2_32

gfortran -fdefault-real-8 -ffree-line-length-none babu.o linalg_mod.o json.o jsonx.o database_m.o udp_windows_m.o connection_m.o sim.o main.o -o main.exe -lws2_32
if errorlevel 1 goto :error

echo.
echo ================================================================================
echo BUILD SUCCESSFUL!k
echo ================================================================================
echo.
echo Run the simulator with: main.exe input.json
echo.
goto :end

:error
echo.
echo ================================================================================
echo BUILD FAILED!
echo ================================================================================
echo.
pause
exit /b 1
:end