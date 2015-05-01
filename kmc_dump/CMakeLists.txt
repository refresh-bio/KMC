set (KMC_DUMP_INCLUDES nc_utils.h)
set (KMC_DUMP_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/kmc_dump.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/nc_utils.cpp)

if (${KMC_BUILD_RUNTIME})
	add_executable (kmc_dump ${KMC_DUMP_INCLUDES} ${KMC_DUMP_SOURCES})
	target_link_libraries (kmc_dump kmc_api)
	install (TARGETS kmc_dump RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR}/kmc)
endif (${KMC_BUILD_RUNTIME})
