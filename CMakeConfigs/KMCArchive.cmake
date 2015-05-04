if (${KMC_INSTALL_ARCHIVE})
	install (FILES ${KMC_SOURCE_DIR}/definitions.h DESTINATION ${INSTALL_INCLUDE_DIR}/kmc)
	export (KMC ${INSTALL_INCLUDE_DIR}/definitions.h)
endif (${KMC_INSTALL_ARCHIVE})
