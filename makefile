all: kmc kmc_dump kmc_tools py_kmc_api

KMC_BIN_DIR = bin
KMC_MAIN_DIR = kmer_counter
KMC_CLI_DIR = kmc-CLI
KMC_API_DIR = kmc_api
KMC_DUMP_DIR = kmc_dump
KMC_TOOLS_DIR = kmc_tools
PY_KMC_API_DIR = py_kmc_api

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14
CLINK	= -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14

KMC_TOOLS_CFLAGS	= -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14
KMC_TOOLS_CLINK	= -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14


LIB_KMC_CORE = $(KMC_BIN_DIR)/libkmc_core.a

KMC_CLI_OBJS = \
$(KMC_CLI_DIR)/kmc.o

KMC_CORE_OBJS = \
$(KMC_MAIN_DIR)/mmer.o \
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
$(KMC_MAIN_DIR)/kff_writer.o \
$(KMC_MAIN_DIR)/kmc_runner.o

RADULS_OBJS = \
$(KMC_MAIN_DIR)/raduls_sse2.o \
$(KMC_MAIN_DIR)/raduls_sse41.o \
$(KMC_MAIN_DIR)/raduls_avx2.o \
$(KMC_MAIN_DIR)/raduls_avx.o

KMC_LIBS = \
$(KMC_MAIN_DIR)/libs/libz.a \
$(KMC_MAIN_DIR)/libs/libbz2.a

KMC_DUMP_OBJS = \
$(KMC_DUMP_DIR)/nc_utils.o \
$(KMC_DUMP_DIR)/kmc_dump.o

KMC_API_OBJS = \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

KMC_API_SRC_FILES = $(wildcard $(KMC_API_DIR)/*.cpp)
PY_KMC_API_OBJS = $(patsubst $(KMC_API_DIR)/%.cpp,$(PY_KMC_API_DIR)/%.o,$(KMC_API_SRC_FILES))

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
$(KMC_TOOLS_DIR)/kff_info_reader.o \
$(KMC_MAIN_DIR)/kff_writer.o

KMC_TOOLS_LIBS = \
$(KMC_TOOLS_DIR)/libs/libz.a \
$(KMC_TOOLS_DIR)/libs/libbz2.a

$(KMC_CLI_OBJS) $(KMC_CORE_OBJS) $(KMC_DUMP_OBJS) $(KMC_API_OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(KMC_TOOLS_OBJS): %.o: %.cpp
	$(CC) $(KMC_TOOLS_CFLAGS) -c $< -o $@

$(KMC_MAIN_DIR)/raduls_sse2.o: $(KMC_MAIN_DIR)/raduls_sse2.cpp
	$(CC) $(CFLAGS) -msse2 -c $< -o $@
$(KMC_MAIN_DIR)/raduls_sse41.o: $(KMC_MAIN_DIR)/raduls_sse41.cpp
	$(CC) $(CFLAGS) -msse4.1 -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx.o: $(KMC_MAIN_DIR)/raduls_avx.cpp
	$(CC) $(CFLAGS) -mavx -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx2.o: $(KMC_MAIN_DIR)/raduls_avx2.cpp
	$(CC) $(CFLAGS) -mavx2 -c $< -o $@

$(LIB_KMC_CORE): $(KMC_CORE_OBJS) $(RADULS_OBJS)
	-mkdir -p $(KMC_BIN_DIR)
	ar rcs $@ $^

kmc: $(KMC_CLI_OBJS) $(LIB_KMC_CORE)
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(CLINK) -o $(KMC_BIN_DIR)/$@ $^ $(KMC_LIBS)

kmc_dump: $(KMC_DUMP_OBJS) $(KMC_API_OBJS)
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(CLINK) -o $(KMC_BIN_DIR)/$@ $^

kmc_tools: $(KMC_TOOLS_OBJS) $(KMC_API_OBJS)
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(KMC_TOOLS_CLINK) -o $(KMC_BIN_DIR)/$@ $^ $(KMC_TOOLS_LIBS)

$(PY_KMC_API_DIR)/%.o: $(KMC_API_DIR)/%.cpp
	$(CC) -c -fPIC -Wall -O3 -m64 -std=c++11 $^ -o $@

py_kmc_api: $(PY_KMC_API_OBJS) $(PY_KMC_API_OBJS)
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) -fPIC -Wall -shared -std=c++11 -O3 $(PY_KMC_API_DIR)/py_kmc_api.cpp $(PY_KMC_API_OBJS) \
	-I $(KMC_API_DIR) \
	-I $(PY_KMC_API_DIR)/libs/pybind11/include \
	-I `python3 -c "import sysconfig;print(sysconfig.get_paths()['include'])"` \
	-o $(KMC_BIN_DIR)/$@`python3-config --extension-suffix`

clean:
	-rm -f $(KMC_MAIN_DIR)/*.o
	-rm -f $(KMC_API_DIR)/*.o
	-rm -f $(KMC_DUMP_DIR)/*.o
	-rm -f $(KMC_TOOLS_DIR)/*.o
	-rm -f $(PY_KMC_API_DIR)/*.o
	-rm -f $(PY_KMC_API_DIR)/*.so
	-rm -rf bin
