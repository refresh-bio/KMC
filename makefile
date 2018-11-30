all: kmc
	
KMC_BIN_DIR = bin
KMC_MAIN_DIR = kmer_counter
KMC_API_DIR = kmc_api
KMC_DUMP_DIR = kmc_dump
KMC_TOOLS_DIR = kmc_tools

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11
CLINK	= -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11

KMC_TOOLS_CFLAGS	= -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14
KMC_TOOLS_CLINK	= -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14

KMC_OBJS = \
$(KMC_MAIN_DIR)/kmer_counter.o \
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
$(KMC_MAIN_DIR)/kb_collector.o 
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

KMC_TOOLS_OBJS = \
$(KMC_TOOLS_DIR)/kmc_header.o \
$(KMC_TOOLS_DIR)/kmc_tools.o \
$(KMC_TOOLS_DIR)/nc_utils.o \
$(KMC_TOOLS_DIR)/parameters_parser.o \
$(KMC_TOOLS_DIR)/parser.o \
$(KMC_TOOLS_DIR)/tokenizer.o \
$(KMC_TOOLS_DIR)/fastq_filter.o \
$(KMC_TOOLS_DIR)/fastq_reader.o \
$(KMC_TOOLS_DIR)/fastq_writer.o \
$(KMC_TOOLS_DIR)/percent_progress.o

KMC_TOOLS_LIBS = \
$(KMC_TOOLS_DIR)/libs/libz.a \
$(KMC_TOOLS_DIR)/libs/libbz2.a 

$(KMC_OBJS) $(KMC_DUMP_OBJS) $(KMC_API_OBJS): %.o: %.cpp
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
	
kmc: $(KMC_OBJS) $(RADULS_OBJS)
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(CLINK) -o $(KMC_BIN_DIR)/$@ $^ $(KMC_LIBS)
kmc_dump: $(KMC_DUMP_OBJS) $(KMC_API_OBJS)
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(CLINK) -o $(KMC_BIN_DIR)/$@ $^
	
kmc_tools: $(KMC_TOOLS_OBJS) $(KMC_API_OBJS)
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(KMC_TOOLS_CLINK) -o $(KMC_BIN_DIR)/$@ $^ $(KMC_TOOLS_LIBS)
	
clean:
	-rm $(KMC_MAIN_DIR)/*.o
	-rm $(KMC_API_DIR)/*.o
	-rm $(KMC_DUMP_DIR)/*.o
	-rm $(KMC_TOOLS_DIR)/*.o
	-rm -rf bin

all: kmc kmc_dump kmc_tools
