cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

project (KMC C CXX)
find_package (OpenMP REQUIRED)
find_package (Threads REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set (KMC_API_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kmc_api)
set (KMC_COUNTER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter)
set (KMC_DUMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kmc_dump)
set (KMC_DUMP_SAMPLE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kmc_sample_dump)

### KMC
option (KMC_ADD_ASM_LIB "Use ASM library" ON)
option (KMC_ADD_DEVELOP "Add KMC developer tools" OFF)

set (KMC_COUNTER_INCLUDES
	${KMC_COUNTER_DIR}/asmlib_wrapper.h
	${KMC_COUNTER_DIR}/bkb_merger.h
	${KMC_COUNTER_DIR}/bkb_reader.h
	${KMC_COUNTER_DIR}/bkb_sorter.h
	${KMC_COUNTER_DIR}/bkb_subbin.h
	${KMC_COUNTER_DIR}/bkb_uncompactor.h
	${KMC_COUNTER_DIR}/bkb_writer.h
	${KMC_COUNTER_DIR}/defs.h
	${KMC_COUNTER_DIR}/fastq_reader.h
	${KMC_COUNTER_DIR}/kb_collector.h
	${KMC_COUNTER_DIR}/kb_completer.h
	${KMC_COUNTER_DIR}/kb_reader.h
	${KMC_COUNTER_DIR}/kb_sorter.h
	${KMC_COUNTER_DIR}/kb_storer.h
	${KMC_COUNTER_DIR}/kmc.h
	${KMC_COUNTER_DIR}/kmer.h
	${KMC_COUNTER_DIR}/kxmer_set.h
	${KMC_COUNTER_DIR}/mem_disk_file.h
	${KMC_COUNTER_DIR}/meta_oper.h
	${KMC_COUNTER_DIR}/mmer.h
	${KMC_COUNTER_DIR}/params.h
	${KMC_COUNTER_DIR}/queues.h
	${KMC_COUNTER_DIR}/radix.h
	${KMC_COUNTER_DIR}/rev_byte.h
	${KMC_COUNTER_DIR}/s_mapper.h
	${KMC_COUNTER_DIR}/splitter.h
	${KMC_COUNTER_DIR}/stdafx.h
	${KMC_COUNTER_DIR}/targetver.h
	${KMC_COUNTER_DIR}/timer.h)
set (KMC_COUNTER_SOURCES
	${KMC_COUNTER_DIR}/bkb_writer.cpp
	${KMC_COUNTER_DIR}/bkb_reader.cpp
	${KMC_COUNTER_DIR}/fastq_reader.cpp
	${KMC_COUNTER_DIR}/kb_completer.cpp
	${KMC_COUNTER_DIR}/kb_storer.cpp
	${KMC_COUNTER_DIR}/kmer.cpp
	${KMC_COUNTER_DIR}/kmer_counter.cpp
	${KMC_COUNTER_DIR}/mmer.cpp
	${KMC_COUNTER_DIR}/mem_disk_file.cpp
	${KMC_COUNTER_DIR}/radix.cpp
	${KMC_COUNTER_DIR}/rev_byte.cpp
	${KMC_COUNTER_DIR}/stdafx.cpp
	${KMC_COUNTER_DIR}/timer.cpp)

if (KMC_ADD_DEVELOP)
	set (KMC_COUNTER_INCLUDES ${KMC_COUNTER_INCLUDES} ${KMC_COUNTER_DIR}/develop.h)
	set (KMC_COUNTER_SOURCES ${KMC_COUNTER_SOURCES} ${KMC_COUNTER_DIR}/develop.cpp)
endif (KMC_ADD_DEVELOP)

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

add_executable (kmc ${KMC_COUNTER_SOURCES} ${KMC_COUNTER_INCLUDES})
target_link_libraries (kmc ${KMC_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

### KMC API
set (KMC_API_INCLUDES
	${KMC_API_DIR}/kmc_file.h
	${KMC_API_DIR}/kmer_api.h
	${KMC_API_DIR}/kmer_defs.h
	${KMC_API_DIR}/mmer.h
	${KMC_API_DIR}/stdafx.h)
set (KMC_API_SOURCES
	${KMC_API_DIR}/kmc_file.cpp
	${KMC_API_DIR}/kmer_api.cpp
	${KMC_API_DIR}/mmer.cpp)
add_library (kmc_api ${KMC_API_INCLUDES} ${KMC_API_SOURCES})

### KMC DUMP
set (KMC_DUMP_INCLUDES
	${KMC_DUMP_DIR}/nc_utils.h
	${KMC_DUMP_DIR}/stdafx.h
	${KMC_DUMP_DIR}/targetver.h)
set (KMC_DUMP_SOURCES
	${KMC_DUMP_DIR}/kmc_dump.cpp
	${KMC_DUMP_DIR}/nc_utils.cpp
	${KMC_DUMP_DIR}/stdafx.cpp)

add_executable (kmc_dump ${KMC_DUMP_INCLUDES} ${KMC_DUMP_SOURCES})

### KMC DUMP SAMPLE
set (KMC_DUMP_SAMPLE_INCLUDES
	${KMC_DUMP_SAMPLE_DIR}/stdafx.h
	${KMC_DUMP_SAMPLE_DIR}/targetver.h)
set (KMC_DUMP_SAMPLE_SOURCES
	${KMC_DUMP_SAMPLE_DIR}/kmc_dump_sample.cpp
	${KMC_DUMP_SAMPLE_DIR}/stdafx.cpp)

add_executable (kmc_dump ${KMC_DUMP_SAMPLE_INCLUDES} ${KMC_DUMP_SAMPLE_SOURCES})
