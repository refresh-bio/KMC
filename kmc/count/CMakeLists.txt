add_executable (kmc_count kmc_count.cpp)
target_link_libraries (kmc_count kmc_core)
set_target_properties (kmc_count PROPERTIES VERSION ${KMC_VERSION})
set_target_properties (kmc_count PROPERTIES SOVERSION ${KMC_SOVERSION})

install (TARGETS kmc_count
	EXPORT KMCTargets
	RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR}
	ARCHIVE DESTINATION ${INSTALL_ARCHIVE_DIR}/kmc/count
	LIBRARY DESTINATION ${INSTALL_LIBRARY_DIR}/kmc/count
	INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}/kmc/count)
