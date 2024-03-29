all: kmc kmc_dump kmc_tools

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
UNAME_P := $(shell uname -p)

KMC_MAIN_DIR = kmc_core
KMC_CLI_DIR = kmc_CLI
KMC_API_DIR = kmc_api
KMC_DUMP_DIR = kmc_dump
KMC_TOOLS_DIR = kmc_tools

OUT_BIN_DIR = bin
OUT_INCLUDE_DIR = include
CPU_FLAGS =
STATIC_CFLAGS = 
STATIC_LFLAGS = 

CC 	= gcc
CXX = g++

ifeq ($(UNAME_M),aarch64)
	ARCH_LFLAGS = -pthread	
else
	ARCH_LFLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
endif

CXXFLAGS	= -fPIC -Wall -O3 -fsigned-char -m64 -std=c++14
LDFLAGS		= -lz -lm $(ARCH_LFLAGS) -std=c++14

KMC_CLI_OBJS = \
$(KMC_CLI_DIR)/kmc.o

KFF_OBJS = \
$(KMC_MAIN_DIR)/kff_writer.o

KMC_CORE_OBJS = \
$(KMC_MAIN_DIR)/mem_disk_file.o \
$(KMC_MAIN_DIR)/rev_byte.o \
$(KMC_MAIN_DIR)/bkb_writer.o \
$(KMC_MAIN_DIR)/cpu_info.o \
$(KMC_MAIN_DIR)/bkb_reader.o \
$(KMC_MAIN_DIR)/fastq_reader.o \
$(KMC_MAIN_DIR)/timer.o \
$(KMC_MAIN_DIR)/develop.o \
$(KMC_MAIN_DIR)/kb_completer.o \
$(KMC_MAIN_DIR)/kb_storer.o \
$(KMC_MAIN_DIR)/kmer.o \
$(KMC_MAIN_DIR)/splitter.o \
$(KMC_MAIN_DIR)/kb_collector.o \
$(KMC_MAIN_DIR)/kmc_runner.o

ifeq ($(UNAME_M),aarch64)
	RADULS_OBJS = \
	$(KMC_MAIN_DIR)/raduls_neon.o
else
	RADULS_OBJS = \
	$(KMC_MAIN_DIR)/raduls_sse41.o \
	$(KMC_MAIN_DIR)/raduls_avx2.o \
	$(KMC_MAIN_DIR)/raduls_avx.o
endif

LIB_KMC_CORE = $(OUT_BIN_DIR)/libkmc_core.a

KMC_DUMP_OBJS = \
$(KMC_DUMP_DIR)/nc_utils.o \
$(KMC_DUMP_DIR)/kmc_dump.o

KMC_API_OBJS = \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

KMC_API_SRC_FILES = $(wildcard $(KMC_API_DIR)/*.cpp)

KMC_TOOLS_OBJS = \
$(KMC_TOOLS_DIR)/kmer_file_header.o \
$(KMC_TOOLS_DIR)/kmc_tools.o \
$(KMC_TOOLS_DIR)/nc_utils.o \
$(KMC_TOOLS_DIR)/parameters_parser.o \
$(KMC_TOOLS_DIR)/parser.o \
$(KMC_TOOLS_DIR)/tokenizer.o \
$(KMC_TOOLS_DIR)/fastq_filter.o \
$(KMC_TOOLS_DIR)/fastq_reader.o \
$(KMC_TOOLS_DIR)/fastq_writer.o \
$(KMC_TOOLS_DIR)/percent_progress.o \
$(KMC_TOOLS_DIR)/kff_info_reader.o

$(KMC_MAIN_DIR)/raduls_sse41.o: $(KMC_MAIN_DIR)/raduls_sse41.cpp
	$(CXX) $(CXXFLAGS) -msse4.1 -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx.o: $(KMC_MAIN_DIR)/raduls_avx.cpp
	$(CXX) $(CXXFLAGS) -mavx -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx2.o: $(KMC_MAIN_DIR)/raduls_avx2.cpp
	$(CXX) $(CXXFLAGS) -mavx2 -c $< -o $@

$(KMC_MAIN_DIR)/raduls_neon.o: $(KMC_MAIN_DIR)/raduls_neon.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


libkmc_core.a: $(KMC_CORE_OBJS) $(RADULS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS)
	mkdir -p $(OUT_INCLUDE_DIR)
	cp $(KMC_MAIN_DIR)/kmc_runner.h $(OUT_INCLUDE_DIR)/kmc_runner.h
	mkdir -p $(OUT_BIN_DIR)
	ar rcs $(OUT_BIN_DIR)/$@ $^

kmc: $(KMC_CLI_OBJS) $(KMC_CORE_OBJS) $(RADULS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $(OUT_BIN_DIR)/$@ $^ $(LDFLAGS)

kmc_dump: $(KMC_DUMP_OBJS) $(KMC_API_OBJS)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $(OUT_BIN_DIR)/$@ $^ $(LDFLAGS)

kmc_tools: $(KMC_TOOLS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $(OUT_BIN_DIR)/$@ $^ $(LDFLAGS)

clean:
	-rm -f $(KMC_MAIN_DIR)/*.o
	-rm -f $(KMC_CLI_DIR)/*.o
	-rm -f $(KMC_API_DIR)/*.o
	-rm -f $(KMC_DUMP_DIR)/*.o
	-rm -f $(KMC_TOOLS_DIR)/*.o
	-rm -rf $(OUT_BIN_DIR)
	-rm -rf $(OUT_INCLUDE_DIR)
