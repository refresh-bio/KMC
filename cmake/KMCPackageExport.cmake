include (CMakePackageConfigHelpers)
write_basic_package_version_file (
	${CMAKE_CURRENT_BINARY_DIR}/cmake/KMCConfigVersion.cmake
	VERSION ${KMC_VERSION}
	COMPATIBILITY AnyNewerVersion)

export (EXPORT KMCTargets
	FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/KMCTargets.cmake)

configure_file (
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/KMCConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/KMCConfig.cmake
	COPYONLY)

install (FILES
	${CMAKE_CURRENT_BINARY_DIR}/cmake/KMCConfigVersion.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/KMCTargets.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/KMCConfig.cmake
	DESTINATION ${INSTALL_PACKAGE_DIR}
	COMPONENT devel)
