if HAVE_AVX2
RADULS=src/kmc_core/raduls_avx2.cpp
endif
if HAVE_AVX
RADULS=src/kmc_core/raduls_avx.cpp
endif
if HAVE_SSE41
RADULS=src/kmc_core/raduls_sse41.cpp
endif
# only one of these will every be defined as true, per the configure - allows non-recursive if
bin_PROGRAMS = kmc kmc_dump kmc_tools

kmc_CXXFLAGS = $(AM_CXXFLAGS) $(KMC_CXXFLAGS)
kmc_dump_CXXFLAGS = $(AM_CXXFLAGS) $(KMC_CXXFLAGS)
kmc_tools_CXXFLAGS = $(AM_CXXFLAGS) $(KMC_CXXFLAGS)
kmc_LDFLAGS = $(AM_LDFLAGS) $(KMC_LDFLAGS)
kmc_dump_LDFLAGS = $(AM_LDFLAGS) $(KMC_LDFLAGS)
kmc_tools_LDFLAGS = $(AM_LDFLAGS) $(KMC_LDFLAGS)

KMC_CORE_SRC=src/kmc_core/mem_disk_file.cpp src/kmc_core/rev_byte.cpp src/kmc_core/bkb_writer.cpp src/kmc_core/bkb_reader.cpp src/kmc_core/fastq_reader.cpp src/kmc_core/timer.cpp src/kmc_core/develop.cpp src/kmc_core/kb_completer.cpp src/kmc_core/kb_storer.cpp src/kmc_core/kmer.cpp src/kmc_core/splitter.cpp src/kmc_core/kb_collector.cpp src/kmc_core/kmc_runner.cpp
KMC_TOOLS_SRC=src/kmc_tools/kmer_file_header.cpp src/kmc_tools/nc_utils.cpp src/kmc_tools/parameters_parser.cpp src/kmc_tools/parser.cpp src/kmc_tools/tokenizer.cpp src/kmc_tools/fastq_filter.cpp src/kmc_tools/fastq_reader.cpp src/kmc_tools/fastq_writer.cpp src/kmc_tools/percent_progress.cpp src/kmc_tools/kff_info_reader.cpp
KMC_API_SRC=src/kmc_api/mmer.cpp src/kmc_api/kmc_file.cpp src/kmc_api/kmer_api.cpp
KFF_SRC=src/kmc_core/kff_writer.cpp

kmc_SOURCES = ${RADULS} ${KMC_CORE_SRC} ${KMC_API_SRC} ${KFF_SRC} src/kmc_cli/kmc.cpp

kmc_dump_SOURCES = ${KMC_API_SRC} src/kmc_dump/nc_utils.cpp src/kmc_dump/kmc_dump.cpp

kmc_tools_SOURCES = ${KMC_TOOLS_SRC} ${KMC_API_SRC} ${KFF_SRC} src/kmc_tools/kmc_tools.cpp

lib_LIBRARIES = libkmc_core.a
libkmc_core_a_SOURCES = ${RADULS} ${KMC_CORE_SRC} ${KMC_API_SRC} ${KFF_SRC}
