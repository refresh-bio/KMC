cmake_minimum_required (VERSION 2.8 FATAL_ERROR)
include (CPack)
include (CMakePackageConfigHelpers)

project (kmc LANGUAGES C CXX)
set (KMC_VERSION_MAJOR 2)
set (KMC_VERSION_MINOR 2)
set (KMC_VERSION_PATCH 0)
set (KMC_VERSION_TWEAK 20150415)
set (KMC_VERSION ${KMC_VERSION_MAJOR}.${KMC_VERSION_MINOR}.${KMC_VERSION_PATCH})
message (STATUS "KMC Version: ${KMC_VERSION}")

add_definitions (${OpenMP_C_FLAGS})
add_definitions (${OpenMP_CXX_FLAGS})
add_definitions ("-std=c++11")

find_package (OpenMP REQUIRED)
find_package (Threads REQUIRED)

option (KMC_BUILD_WITH_ASM "Build KMC with ASM library included" ON)
option (KMC_INSTALL_RUNTIME "Build and install KMC executables" ON)
option (KMC_INSTALL_ARCHIVE "Build and install KMC headers and as a static library" OFF)
option (KMC_INSTALL_PACKAGE "Build and install KMC package configuration for FindPacakge" OFF)
option (KMC_ENABLE_DEBUG "Add KMC debug mode" OFF)
option (KMC_ENABLE_DEVELOP "Add KMC developer mode" OFF)
mark_as_advanced (KMC_BUILD_ARCHIVE KMC_ENABLE_DEBUG KMC_ENABLE_DEVELOP)

set (LIBRARY_SUFFIX "" CACHE STRING "Suffix to append the installation directory")
set (INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set (INSTALL_RUNTIME_DIR bin${LIBRARY_SUFFIX} CACHE PATH "Installation directory for executables")
set (INSTALL_ARCHIVE_DIR lib${LIBRARY_SUFFIX} CACHE PATH "Installation directory for static libraries")
set (INSTALL_PACKAGE_DIR lib${LIBRARY_SUFFIX}/cmake CACHE PATH "Installation directory for cmake package configuration")
mark_as_advanced (INSTALL_INCLUDE_DIR INSTALL_RUNTIME_DIR INSTALL_ARCHIVE_DIR)

include (CMakeConfigs/KMCPreprocessor.cmake)
configure_file (${CMAKE_SOURCE_DIR}/CMakeConfigs/definitions.h.in ${CMAKE_SOURCE_DIR}/kmc/definitions.h)
if (KMC_INSTALL_ARCHIVE)
	install (FILES ${KMC_SOURCE_DIR}/kmc/definitions.h DESTINATION ${INSTALL_INCLUDE_DIR}/kmc)
	export (KMC ${INSTALL_INCLUDE_DIR}/kmc/definitions.h)
endif (KMC_INSTALL_ARCHIVE)

add_subdirectory (kmc_api)
add_subdirectory (kmc_core)
add_subdirectory (kmc_counter)
add_subdirectory (kmc_dump)
add_subdirectory (kmc_dump_sample)

if (KMC_INSTALL_PACKAGE)
	configure_package_config_file (
		CMakeConfigs/KMCConfig.cmake.in CMakeConfigs/KMCConfig.cmake
		INSTALL_DESTINATION ${INSTALL_PACKAGE_DIR}
		PATH_VARS INCLUDE_INSTALL_DIR)
	install (TARGETS KMC
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		EXPORT KMCTargets)
	install (EXPORT KMCTargets
		FILE KMCTargets.cmake
		NAMESPACE KMC::
		DESTINATION ${INSTALL_PACKAGE_DIR})
	install (FILES "${CMAKE_CURRENT_BINARY_DIR}/CMakeConfigs/KMCConfig.cmake"
		DESTINATION ${INSTALL_PACKAGE_DIR})
endif (KMC_INSTALL_PACKAGE)
