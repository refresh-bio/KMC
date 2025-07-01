ACLOCAL_AMFLAGS=-I m4
RADULS=src/kmc_core/raduls_impl.h \
	   src/kmc_core/intr_copy.h \
	   src/kmc_core/raduls.h \
	   src/kmc_core/radix.h \
	   src/kmc_core/small_sort.h \
	   src/kmc_core/small_k_buf.h

if HAVE_AVX2
RADULS+=src/kmc_core/raduls_avx2.cpp
endif
if HAVE_AVX
RADULS+=src/kmc_core/raduls_avx.cpp
endif
if HAVE_SSE41
RADULS+=src/kmc_core/raduls_sse41.cpp
endif
# only one of these will every be defined as true, per the configure - allows non-recursive if
bin_PROGRAMS = kmc kmc_dump kmc_tools

kmc_CXXFLAGS = $(AM_CXXFLAGS) $(KMC_CXXFLAGS)
kmc_dump_CXXFLAGS = $(AM_CXXFLAGS) $(KMC_CXXFLAGS)
kmc_tools_CXXFLAGS = $(AM_CXXFLAGS) $(KMC_CXXFLAGS)
kmc_LDFLAGS = $(AM_LDFLAGS) $(KMC_LDFLAGS)
kmc_dump_LDFLAGS = $(AM_LDFLAGS) $(KMC_LDFLAGS)
kmc_tools_LDFLAGS = $(AM_LDFLAGS) $(KMC_LDFLAGS)

KMC_CORE_SRC=src/kmc_core/mem_disk_file.cpp \
			 src/kmc_core/mem_disk_file.h \
			 src/kmc_core/rev_byte.cpp \
			 src/kmc_core/rev_byte.h \
			 src/kmc_core/bkb_writer.cpp \
			 src/kmc_core/bkb_writer.h \
			 src/kmc_core/bkb_reader.cpp \
			 src/kmc_core/bkb_reader.h \
			 src/kmc_core/fastq_reader.cpp \
			 src/kmc_core/fastq_reader.h \
			 src/kmc_core/timer.cpp \
			 src/kmc_core/timer.h \
			 src/kmc_core/develop.cpp \
			 src/kmc_core/develop.h \
			 src/kmc_core/kb_completer.cpp \
			 src/kmc_core/kb_completer.h \
			 src/kmc_core/kb_storer.cpp \
			 src/kmc_core/kb_storer.h \
			 src/kmc_core/kmer.cpp \
			 src/kmc_core/kmer.h \
			 src/kmc_core/splitter.cpp \
			 src/kmc_core/splitter.h \
			 src/kmc_core/kb_collector.cpp \
			 src/kmc_core/kb_collector.h \
			 src/kmc_core/kmc_runner.cpp \
			 src/kmc_core/kmc_runner.h \
			 src/kmc_core/first_dispatch.h \
	   		 src/kmc_core/meta_oper.h \
	   		 src/kmc_core/queues.h \
	   		 src/kmc_core/critical_error_handler.h \
	   		 src/kmc_core/thread_cancellation_exception.h \
			 src/kmc_core/params.h \
			 src/kmc_core/binary_reader.h \
			 src/kmc_core/bkb_uncompactor.h \
			 src/kmc_core/bkb_merger.h \
			 src/kmc_core/bam_utils.h \
			 src/kmc_core/bkb_sorter.h \
			 src/kmc_core/bkb_subbin.h \
			 src/kmc_core/defs.h \
			 src/kmc_core/exception_aware_thread.h \
			 src/kmc_core/kb_sorter.h \
			 src/kmc_core/kb_reader.h \
			 src/kmc_core/s_mapper.h \
			 src/kmc_core/libs/ntHash/ntHashWrapper.h \
			 src/kmc_core/libs/ntHash/ntHash.hpp \
			 src/kmc_core/tmp_files_owner.h \
			 src/kmc_core/kmc.h \
			 src/kmc_core/percent_progress.h \
			 src/kmc_core/kxmer_set.h

KMC_TOOLS_SRC=src/kmc_tools/kmer_file_header.cpp \
			  src/kmc_tools/kmer_file_header.h \
			  src/kmc_tools/nc_utils.cpp \
			  src/kmc_tools/nc_utils.h \
			  src/kmc_tools/parameters_parser.cpp \
			  src/kmc_tools/parameters_parser.h \
			  src/kmc_tools/parser.cpp \
			  src/kmc_tools/parser.h \
			  src/kmc_tools/tokenizer.cpp \
			  src/kmc_tools/tokenizer.h \
			  src/kmc_tools/fastq_filter.cpp \
			  src/kmc_tools/fastq_filter.h \
			  src/kmc_tools/fastq_reader.cpp \
			  src/kmc_tools/fastq_reader.h \
			  src/kmc_tools/fastq_writer.cpp \
			  src/kmc_tools/fastq_writer.h \
			  src/kmc_tools/percent_progress.cpp \
			  src/kmc_tools/percent_progress.h \
			  src/kmc_tools/kff_info_reader.cpp \
			  src/kmc_tools/kff_info_reader.h \
			  src/kmc_tools/check_kmer.h \
			  src/kmc_tools/kmer.h \
			  src/kmc_tools/defs.h \
			  src/kmc_tools/meta_oper.h \
			  src/kmc_tools/expression_node.h \
			  src/kmc_tools/operations.h \
			  src/kmc_tools/bundle.h \
			  src/kmc_tools/db_reader_factory.h \
			  src/kmc_tools/db_writer.h \
			  src/kmc_tools/db_writer_factory.h \
			  src/kmc_tools/dump_writer.h \
			  src/kmc_tools/kmc1_db_reader.h \
			  src/kmc_tools/kmc1_db_writer.h \
			  src/kmc_tools/kmc2_db_reader.h \
			  src/kmc_tools/kff_db_reader.h \
			  src/kmc_tools/kff_db_writer.h \
			  src/kmc_tools/kff_random_access.h \
			  src/kmc_tools/kff_kmc2_reader_utils.h \
			  src/kmc_tools/queues.h \
			  src/kmc_tools/output_parser.h \
			  src/kmc_tools/timer.h \
			  src/kmc_tools/histogram_writer.h \
			  src/kmc_tools/config.h
KMC_API_SRC=src/kmc_api/mmer.cpp \
			src/kmc_api/mmer.h \
			src/kmc_api/kmc_file.cpp \
			src/kmc_api/kmc_file.h \
			src/kmc_api/kmer_api.cpp \
			src/kmc_api/kmer_api.h \
			src/kmc_api/kmer_defs.h
KFF_SRC=src/kmc_core/kff_writer.cpp \
	    src/kmc_core/kff_writer.h

kmc_SOURCES = ${RADULS} ${KMC_CORE_SRC} ${KMC_API_SRC} ${KFF_SRC} src/kmc_cli/kmc.cpp

kmc_dump_SOURCES = ${KMC_API_SRC} src/kmc_dump/nc_utils.cpp src/kmc_dump/nc_utils.h src/kmc_dump/kmc_dump.cpp

kmc_tools_SOURCES = ${KMC_TOOLS_SRC} ${KMC_API_SRC} ${KFF_SRC} src/kmc_tools/kmc_tools.cpp

#lib_LIBRARIES = libkmc_core.a
include_HEADERS = src/kmc_core/kmc_runner.h
#libkmc_core_a_SOURCES = ${RADULS} ${KMC_CORE_SRC} ${KMC_API_SRC} ${KFF_SRC}
lib_LTLIBRARIES = libkmc_core.la
libkmc_core_la_SOURCES = ${RADULS} ${KMC_CORE_SRC} ${KMC_API_SRC} ${KFF_CRC}
