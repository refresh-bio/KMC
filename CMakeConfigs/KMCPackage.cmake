set (CPACK_PACKAGE_NAME ${PROJECT_NAME})
set (CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set (CPACK_PACKAGE_VENDOR ${PROJECT_VENDOR})
set (CPACK_PACKAGE_CONTACT ${PROJECT_AUTHOR_EMAIL})
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
set (CPACK_PACKAGE_DESCRIPTION_FILE ${PROJECT_SOURCE_DIR}/README.md)

set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set (CPACK_DEBIAN_PACKAGE_DEPENDS "cmake (>=2.8), gcc (>=4:4.8), g++ (>=4:4.8)")

if (${KMC_INSTALL_PACKAGE})
	install (TARGETS KMC
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		EXPORT KMCTargets)
	install (EXPORT KMCTargets
		FILE KMCTargets.cmake
		NAMESPACE KMC::
		DESTINATION ${INSTALL_PACKAGE_DIR})
	install (FILES "${CMAKE_CURRENT_BINARY_DIR}/CMakeConfigs/KMCConfig.cmake"
		DESTINATION ${INSTALL_PACKAGE_DIR})
endif (${KMC_INSTALL_PACKAGE})

include (CPack)
include (CMakePackageConfigHelpers)