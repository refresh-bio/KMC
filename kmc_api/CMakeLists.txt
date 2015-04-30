set (KMC_API_INCLUDES
	kmc_file.h
	kmer_api.h
	mmer.h)

set (KMC_API_SOURCES
	kmc_file.cpp
	kmer_api.cpp
	mmer.cpp)

add_library (kmc_api STATIC ${KMC_API_INCLUDES} ${KMC_API_SOURCES})
set_target_properties (kmc_api PROPERTIES VERSION "${KMC_VERSION}")

if (${KMC_INSTALL_ARCHIVE})
	install (TARGETS kmc_api EXPORT KMC
		INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}
		ARCHIVE DESTINATION ${INSTALL_ARCHIVE_DIR})
endif (${KMC_INSTALL_ARCHIVE})
