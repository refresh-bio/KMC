set (KMC_CORE_INCLUDES
	asmlib_wrapper.h
	bkb_merger.h
	bkb_reader.h
	bkb_sorter.h
	bkb_subbin.h
	bkb_uncompactor.h
	bkb_writer.h
	defs.h
	fastq_reader.h
	kb_collector.h
	kb_completer.h
	kb_reader.h
	kb_sorter.h
	kb_storer.h
	kmc.h
	kmer.h
	kxmer_set.h
	mem_disk_file.h
	meta_oper.h
	mmer.h
	params.h
	queues.h
	radix.h
	rev_byte.h
	s_mapper.h
	splitter.h
	stdafx.h
	targetver.h
	timer.h)
set (KMC_CORE_SOURCES
	bkb_reader.cpp
	bkb_writer.cpp
	fastq_reader.cpp
	kb_completer.cpp
	kb_storer.cpp
	kmer.cpp
	mem_disk_file.cpp
	mmer.cpp
	radix.cpp
	rev_byte.cpp
	stdafx.cpp
	timer.cpp)

if (KMC_ADD_DEVELOP)
	set (KMC_CORE_INCLUDES ${KMC_CORE_COUNTER_INCLUDES} develop.h)
	set (KMC_CORE_COUNTER_SOURCES ${KMC_CORE_COUNTER_SOURCES} develop.cpp)
endif (KMC_ADD_DEVELOP)

if (WIN32 OR APPLE)
	set (KMC_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/libs/libz.1.2.5.dylib ${CMAKE_CURRENT_SOURCE_DIR}/libs/libbz2.1.0.5.dylib)
	if (${KMC_ADD_ASM_LIB})
		set (KMC_LIBRARIES ${KMC_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/libs/libamac64.a)
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_ASMLIB")
	endif()
else (WIN32 OR APPLE)
	set (KMC_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/libs/libz.a ${CMAKE_CURRENT_SOURCE_DIR}/libs/libbz2.a)
	if (${KMC_ADD_ASM_LIB})
		set (KMC_LIBRARIES ${KMC_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/libs/alibelf64.a)
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_ASMLIB")
	endif ()
endif (WIN32 OR APPLE)

add_library (kmc_core
	${KMC_CORE_INCLUDES}
	${KMC_CORE_SOURCES}
	${CMAKE_THREAD_LIBS_INIT})
add_executable (kmer_counter kmer_counter.cpp)
if (${KMC_ADD_COUNTER})
	target_link_libraries (kmer_counter kmc_core ${KMC_LIBRARIES})
endif (${KMC_ADD_COUNTER})
