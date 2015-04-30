set (KMC_CORE_INCLUDES
	asmlib_wrapper.h
	bkb_merger.h
	bkb_reader.h
	bkb_sorter.h
	bkb_subbin.h
	bkb_uncompactor.h
	bkb_writer.h
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
	timer.h)

set (KMC_CORE_SOURCE_DIR source)
set (KMC_CORE_SOURCES
	${KMC_CORE_SOURCE_DIR}/bkb_reader.cpp
	${KMC_CORE_SOURCE_DIR}/bkb_writer.cpp
	${KMC_CORE_SOURCE_DIR}/fastq_reader.cpp
	${KMC_CORE_SOURCE_DIR}/kb_completer.cpp
	${KMC_CORE_SOURCE_DIR}/kb_storer.cpp
	${KMC_CORE_SOURCE_DIR}/kmer.cpp
	${KMC_CORE_SOURCE_DIR}/mem_disk_file.cpp
	${KMC_CORE_SOURCE_DIR}/mmer.cpp
	${KMC_CORE_SOURCE_DIR}/radix.cpp
	${KMC_CORE_SOURCE_DIR}/rev_byte.cpp
	${KMC_CORE_SOURCE_DIR}/timer.cpp)

if (KMC_ENABLE_DEVELOP)
	set (KMC_CORE_INCLUDES ${KMC_CORE_COUNTER_INCLUDES} develop.h)
	set (KMC_CORE_COUNTER_SOURCES ${KMC_CORE_COUNTER_SOURCES} develop.cpp)
endif (KMC_ENABLE_DEVELOP)

set (KMER_COUNTER_EXTERNAL_INCLUDES
	${CMAKE_CURRENT_SOURCE_DIR}/external/asmlib.h
	${CMAKE_CURRENT_SOURCE_DIR}/external/bzlib.h
	${CMAKE_CURRENT_SOURCE_DIR}/external/bzlib_private.h
	${CMAKE_CURRENT_SOURCE_DIR}/external/zconf.h
	${CMAKE_CURRENT_SOURCE_DIR}/external/zlib.h)

if (WIN32 OR APPLE)
	set (KMER_COUNTER_EXTERNAL_LIBRARIES
		${CMAKE_CURRENT_SOURCE_DIR}/external/libz.1.2.5.dylib
		${CMAKE_CURRENT_SOURCE_DIR}/external/libbz2.1.0.5.dylib)
	if (${KMC_BUILD_WITH_ASM})
		list (APPEND KMER_COUNTER_EXTERNAL_LIBRARIES
			${CMAKE_CURRENT_SOURCE_DIR}/external/libamac64.a)
	endif()
else (WIN32 OR APPLE)
	set (KMER_COUNTER_EXTERNAL_LIBRARIES
		${CMAKE_CURRENT_SOURCE_DIR}/external/libz.a
		${CMAKE_CURRENT_SOURCE_DIR}/external/libbz2.a)
	if (${KMC_BUILD_WITH_ASM})
		list (APPEND KMER_COUNTER_EXTERNAL_LIBRARIES
			${CMAKE_CURRENT_SOURCE_DIR}/external/alibelf64.a)
	endif ()
endif (WIN32 OR APPLE)

add_library (kmc_core STATIC ${KMC_CORE_INCLUDES} ${KMC_CORE_SOURCES} ${CMAKE_THREAD_LIBS_INIT})
set_target_properties (kmc_core PROPERTIES VERSION ${KMC_VERSION})

if (${KMC_INSTALL_ARCHIVE})
	install (TARGETS kmc_core KMC
		INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}
		ARCHIVE DESTINATION ${INSTALL_ARCHIVE_DIR})
endif (${KMC_INSTALL_ARCHIVE})

if (${KMC_INSTALL_RUNTIME})
	add_executable (kmer_counter ${KMER_COUNTER_EXTERNAL_INCLUDES} ${KMC_CORE_SOURCE_DIR}/kmer_counter.cpp)
	target_link_libraries (kmer_counter kmc_core ${KMER_COUNTER_EXTERNAL_LIBRARIES})
	install (TARGETS kmer_counter RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR})
endif (${KMC_INSTALL_RUNTIME})
