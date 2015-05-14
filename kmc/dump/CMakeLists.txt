set (HEADERS nc_utils.h)
set (SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/kmc_dump.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/nc_utils.cpp)

add_executable (kmc_dump ${HEADERS} ${SOURCES})
target_link_libraries (kmc_dump kmc_api)
set_target_properties (kmc_dump PROPERTIES VERSION ${KMC_VERSION})
set_target_properties (kmc_dump PROPERTIES SOVERSION ${KMC_SOVERSION})

install (TARGETS kmc_dump
	EXPORT KMCTargets
	RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR}
	ARCHIVE DESTINATION ${INSTALL_ARCHIVE_DIR}/kmc/dump
	LIBRARY DESTINATION ${INSTALL_LIBRARY_DIR}/kmc/dump
	INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}/kmc/dump)
