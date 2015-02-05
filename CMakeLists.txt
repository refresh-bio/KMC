cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

project (KMC)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

option (ASMLIB "Use ASM library" ON)

set (KMC_BIN_DIR bin)
set (KMC_MAIN_DIR kmer_counter)
set (KMC_API_DIR kmc_api)
set (KMC_DUMP_DIR kmc_dump)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${KMC_BIN_DIR})

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
set (KMC_DUMP_OBJECTS
	${KMC_DUMP_DIR}/nc_utils.cpp
	${KMC_API_DIR}/mmer.cpp
	${KMC_DUMP_DIR}/kmc_dump.cpp
	${KMC_API_DIR}/kmc_file.cpp
	${KMC_API_DIR}/kmer_api.cpp
	)
if (APPLE OR WIN32)
	set (KMC_LIBRARIES
		${CMAKE_SOURCE_DIR}/${KMC_MAIN_DIR}/libs/libz.1.2.5.dylib
		${CMAKE_SOURCE_DIR}/${KMC_MAIN_DIR}/libs/libbz2.1.0.5.dylib
		)
	if (${ASMLIB})
		set (KMC_LIBARIES ${KMC_LIBRARIES} ${CMAKE_SOURCE_DIR}/${KMC_MAIN_DIR}/libs/libamac64.a)
	else ()
		set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -DDISABLE_ASMLIB)
	endif ()
else ()
	set (KMC_LIBRARIES
		${CMAKE_SOURCE_DIR}/${KMC_MAIN_DIR}/libs/libz.a
		${CMAKE_SOURCE_DIR}/${KMC_MAIN_DIR}/libs/libbz2.a
		)
	if (${ASMLIB})
		set (KMC_LIBRARIES ${KMC_LIBRARIES} ${CMAKE_SOURCE_DIR}/${KMC_MAIN_DIR}/libs/alibelf64.a)
	else ()
		set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -DDISABLE_ASMLIB)
	endif ()
endif ()

add_executable (KMC ${KMC_SOURCES})
target_link_libraries (KMC ${KMC_LIBRARIES})
