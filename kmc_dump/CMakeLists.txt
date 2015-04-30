set (KMC_DUMP_INCLUDES nc_utils.h)
set (KMC_DUMP_SOURCES
	kmc_dump.cpp
	nc_utils.cpp)

if (${KMC_BUILD_RUNTIME})
	add_executable (kmc_dump ${KMC_DUMP_INCLUDES} ${KMC_DUMP_SOURCES})
	target_link_libraries (kmc_dump kmc_api)
	install (TARGETS kmc_dump RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR}/kmc)
endif (${KMC_BUILD_RUNTIME})
