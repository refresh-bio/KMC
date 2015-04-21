cmake_minimum_required (VERSION 3.0 FATAL_ERROR)

project (KMC)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set (KMC_API_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kmc_api)
set (KMC_MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter)
set (KMC_DUMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kmc_dump)

find_package (OpenMP QUIET)
if (NOT OPENMP_FOUND)
	message (FATAL_ERROR "KMC requires OpenMP. Please install an OpenMP compatible compiler.")
endif (NOT OPENMP_FOUND)


### KMC
option (KMC_ADD_ASM_LIB "Use ASM library" ON)

include_directories (${KMC_MAIN_DIR})
set (KMC_INCLUDES
	${KMC_MAIN_DIR}/asmlib_wrapper.h
	${KMC_MAIN_DIR}/bkb_merger.h
	${KMC_MAIN_DIR}/bkb_reader.h
	${KMC_MAIN_DIR}/bkb_sorter.h
	${KMC_MAIN_DIR}/bkb_subbin.h
	${KMC_MAIN_DIR}/bkb_uncompactor.h
	${KMC_MAIN_DIR}/bkb_writer.h
	${KMC_MAIN_DIR}/fastq_reader.h
	${KMC_MAIN_DIR}/kb_collector.h
	${KMC_MAIN_DIR}/kb_completer.h
	${KMC_MAIN_DIR}/kb_reader.h
	${KMC_MAIN_DIR}/kb_sorter.h
	${KMC_MAIN_DIR}/kb_storer.h
	${KMC_MAIN_DIR}/kmc.h
	${KMC_MAIN_DIR}/kmer.h
	${KMC_MAIN_DIR}/kxmer_set.h
	${KMC_MAIN_DIR}/mem_disk_file.h
	${KMC_MAIN_DIR}/meta_oper.h
	${KMC_MAIN_DIR}/mmer.h
	${KMC_MAIN_DIR}/params.h
	${KMC_MAIN_DIR}/queues.h
	${KMC_MAIN_DIR}/radix.h
	${KMC_MAIN_DIR}/rev_byte.h
	${KMC_MAIN_DIR}/s_mapper.h
	${KMC_MAIN_DIR}/splitter.h
	${KMC_MAIN_DIR}/stdafx.h
	${KMC_MAIN_DIR}/targetver.h
	${KMC_MAIN_DIR}/timer.h
	)
set (KMC_SOURCES
	${KMC_MAIN_DIR}/kmer_counter.cpp
	${KMC_MAIN_DIR}/mmer.cpp
	${KMC_MAIN_DIR}/mem_disk_file.cpp
	${KMC_MAIN_DIR}/rev_byte.cpp
	${KMC_MAIN_DIR}/bkb_writer.cpp
	${KMC_MAIN_DIR}/bkb_reader.cpp
	${KMC_MAIN_DIR}/fastq_reader.cpp
	${KMC_MAIN_DIR}/timer.cpp
	${KMC_MAIN_DIR}/radix.cpp
	${KMC_MAIN_DIR}/kb_completer.cpp
	${KMC_MAIN_DIR}/kb_storer.cpp
	${KMC_MAIN_DIR}/kmer.cpp
	)

if (WIN32 OR APPLE)
	set (KMC_LIBRARIES ${KMC_MAIN_DIR}/libs/libz.1.2.5.dylib ${KMC_MAIN_DIR}/libs/libbz2.1.0.5.dylib)
	if (${KMC_ADD_ASM_LIB})
		set (KMC_LIBRARIES ${KMC_LIBRARIES} ${KMC_MAIN_DIR}/libs/libamac64.a)
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_ASMLIB")
	endif()
else (WIN32 OR APPLE)
	set (KMC_LIBRARIES ${KMC_MAIN_DIR}/libs/libz.a ${KMC_MAIN_DIR}/libs/libbz2.a)
	if (${KMC_ADD_ASM_LIB})
		set (KMC_LIBRARIES ${KMC_LIBRARIES} ${KMC_MAIN_DIR}/libs/alibelf64.a)
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_ASMLIB")
	endif ()
endif (WIN32 OR APPLE)

add_executable (kmc ${KMC_SOURCES} ${KMC_INCLUDES})
target_link_libraries (kmc ${KMC_LIBRARIES})

### KMC DUMP
set (KMC_DUMP_SOURCES
	${KMC_DUMP_DIR}/nc_utils.cpp
	${KMC_DUMP_DIR}/kmc_dump.cpp
	${KMC_API_DIR}/mmer.cpp
	${KMC_API_DIR}/kmc_file.cpp
	${KMC_API_DIR}/kmer_api.cpp
	)

add_executable (kmc_dump ${KMC_DUMP_SOURCES})
