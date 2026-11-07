# Works with mingw32-make / GNU make on Windows (MinGW-w64) and on Linux/macOS (CLI only).
CXX      = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Iinclude

SRCS = cpp_src/z_algorithm.cpp \
       cpp_src/fasta_parser.cpp \
       cpp_src/tandem_analyzer.cpp \
       cpp_src/statistical_model.cpp \
       cpp_src/disease_registry.cpp \
       cpp_src/classifier.cpp \
       cpp_src/auth_system.cpp \
       cpp_src/pathology_report.cpp \
       cpp_src/web_server.cpp \
       cpp_src/main.cpp

ifeq ($(OS),Windows_NT)
  TARGET  = bin/genomix.exe
  LDFLAGS = -static -static-libgcc -static-libstdc++ -lws2_32 -pthread
  MKBIN   = @if not exist bin mkdir bin
  RM_CMD  = @if exist bin\genomix.exe del /Q bin\genomix.exe
else
  # web_server.cpp uses WinSock, so the Linux/macOS build swaps in a stub (Web UI is Windows-only).
  TARGET  = bin/genomix
  SRCS   := $(filter-out cpp_src/web_server.cpp,$(SRCS)) cpp_src/web_server_stub.cpp
  LDFLAGS = -pthread
  MKBIN   = @mkdir -p bin
  RM_CMD  = @rm -f bin/genomix
endif

all: $(TARGET)

$(TARGET): $(SRCS)
	$(MKBIN)
	$(CXX) $(CXXFLAGS) $(SRCS) $(LDFLAGS) -o $(TARGET)

clean:
	$(RM_CMD)
