set (HEADERS
	kmc_file.h
	kmer_api.h
	mmer.h)

set (SOURCE
	kmc_file.cpp
	kmer_api.cpp
	mmer.cpp)

add_library (kmc_api STATIC ${HEADERS} ${SOURCE})
set_target_properties (kmc_api PROPERTIES VERSION ${KMC_VERSION})
set_target_properties (kmc_api PROPERTIES SOVERSION ${KMC_SOVERSION})

install (TARGETS kmc_api
	EXPORT KMCTargets
	RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR}
	ARCHIVE DESTINATION ${INSTALL_ARCHIVE_DIR}
	LIBRARY DESTINATION ${INSTALL_LIBRARY_DIR}
	INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}/kmc/api)
