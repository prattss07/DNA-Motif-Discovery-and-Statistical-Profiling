CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Iinclude
LDFLAGS = -static -static-libgcc -static-libstdc++ -lws2_32
TARGET = bin/genomix.exe

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

all: $(TARGET)

$(TARGET): $(SRCS)
	@if not exist bin mkdir bin
	$(CXX) $(CXXFLAGS) $(SRCS) $(LDFLAGS) -o $(TARGET)

clean:
	@if exist bin\genomix.exe del /Q bin\genomix.exe
