cmake_minimum_required (VERSION 2.8.11 FATAL_ERROR)
project (KMC C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set (KMC_VERSION_MAJOR 2)
set (KMC_VERSION_MINOR 2)
set (KMC_VERSION_PATCH 0)
set (KMC_VERSION_TWEAK 20150415)
set (KMC_VERSION ${KMC_VERSION_MAJOR}.${KMC_VERSION_MINOR}.${KMC_VERSION_PATCH})
set (KMC_SOVERSION 1)

set (KMC_VENDOR "REFRESH Bioinformatics Group")
set (KMC_DESCRIPTION "a utility designed for counting k-mers in a set of reads from genome sequencing projects")
set (KMC_WEBSITE "http://sun.aei.polsl.pl/REFRESH/index.php?page=projects-kmc")
set (KMC_SOURCE "https://github.com/marekkokot/KMC")
set (KMC_AUTHOR_NAME "Marek Kokot")
set (KMC_AUTHOR_EMAIL "Marek.Kokot@polsl.pl")

find_package (OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")

option (KMC_BUILD_RUNTIME "BUild KMC executables - make take a while." ON)
option (KMC_BUILD_WITH_ASM "Build KMC with ASM library included" ON)
option (KMC_BUILD_WITH_DEBUG "Build KMC with debug mode on" OFF)
option (KMC_BUILD_WITH_DEVELOP "Build KMC with developer mode on" OFF)
mark_as_advanced (KMC_BUILD_WITH_DEBUG KMC_BUILD_WITH_DEVELOP)

set (INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set (INSTALL_LIBRARY_DIR lib CACHE PATH "Installation directory for shared libraries")
set (INSTALL_ARCHIVE_DIR lib CACHE PATH "Installation directory for static libraries")
set (INSTALL_RUNTIME_DIR bin CACHE PATH "Installation directory for executables")
set (INSTALL_PACKAGE_DIR lib/cmake CACHE PATH "Installation directory for cmake package configuration")
mark_as_advanced (INSTALL_INCLUDE_DIR INSTALL_RUNTIME_DIR INSTALL_ARCHIVE_DIR INSTALL_PACKAGE_DIR)

include (KMCBuildOptions)
include (KMCPackageOptions)
include (KMCExternalImports)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory (kmc)

include (KMCPackageExport)
