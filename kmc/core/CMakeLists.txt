set (HEADERS
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
	kmc_typedefs.h
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

set (SOURCES
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
	timer.cpp)

if (KMC_ENABLE_DEVELOP)
	set (HEADERS ${HEADERS} develop.h)
	set (SOURCES ${SOURCES} develop.cpp)
endif (KMC_ENABLE_DEVELOP)

add_library (kmc_core STATIC ${HEADERS} ${SOURCES})
target_link_libraries (kmc_core asm bzip zlib)
set_target_properties (kmc_core PROPERTIES VERSION ${KMC_VERSION})
set_target_properties (kmc_core PROPERTIES SOVERSION ${KMC_SOVERSION})

install (TARGETS kmc_core
	EXPORT KMCTargets
	RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR}
	ARCHIVE DESTINATION ${INSTALL_ARCHIVE_DIR}
	LIBRARY DESTINATION ${INSTALL_LIBRARY_DIR}
	INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}/kmc/core)
