@echo off
cd /d "%~dp0"
echo =======================================================================
echo          Building GenomiX C++ Clinical Platform (Standalone Binary)
echo =======================================================================

set GXX="C:\Program Files\CodeBlocks\MinGW\bin\g++.exe"
if not exist %GXX% (
    set GXX=g++
)

if not exist bin mkdir bin

echo Compiling C++ modules (statically linked for maximum portability)...
%GXX% -std=c++17 -O3 -Iinclude ^
    cpp_src/z_algorithm.cpp ^
    cpp_src/fasta_parser.cpp ^
    cpp_src/tandem_analyzer.cpp ^
    cpp_src/statistical_model.cpp ^
    cpp_src/disease_registry.cpp ^
    cpp_src/classifier.cpp ^
    cpp_src/auth_system.cpp ^
    cpp_src/pathology_report.cpp ^
    cpp_src/web_server.cpp ^
    cpp_src/main.cpp ^
    -static -static-libgcc -static-libstdc++ ^
    -lws2_32 ^
    -o bin/genomix.exe

if %ERRORLEVEL% equ 0 (
    echo [SUCCESS] Statically linked binary created: bin\genomix.exe
) else (
    echo [FAILED] Compilation errors encountered.
)
