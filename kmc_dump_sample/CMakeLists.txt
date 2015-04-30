set (KMC_DUMP_SAMPLE_SOURCES kmc_dump_sample.cpp)

if (${KMC_BUILD_RUNTIME})
	add_executable (kmc_dump_sample ${KMC_DUMP_SAMPLE_SOURCES})
	target_link_libraries (kmc_dump_sample kmc_api)
	install (TARGETS kmc_dump_sample RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR}/kmc)
endif (${KMC_BUILD_RUNTIME})
