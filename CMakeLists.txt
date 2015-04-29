cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

project (KMC LANGUAGES C CXX)
set(PROJECT_VERSION_MAJOR 2)
set(PROJECT_VERSION_MINOR 2)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_VERSION_TWEAK 20150415)
set(PROJECT_VERSION
	${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}.${PROJECT_VERSION_TWEAK})

find_package (OpenMP REQUIRED)
find_package (Threads REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

option (KMC_COUNTER "Build the kmer counter" ON)
option (KMC_DISABLE_ASM_LIB "Use ASM library" OFF)
option (KMC_ENABLE_DEVELOP "Add KMC developer mode" OFF)
option (KMC_ENABLE_DEBUG "Add KMC debub mode" OFF)

include (cmake/config.cmake)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory (kmc_api)
add_subdirectory (kmc_dump)
add_subdirectory (kmc_dump_sample)
add_subdirectory (kmer_counter)
