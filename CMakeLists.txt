cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

project (KMC C CXX)
find_package (OpenMP REQUIRED)
find_package (Threads REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

option (KMC_ADD_COUNTER "Build the kmer counter" ON)
option (KMC_ADD_ASM_LIB "Use ASM library" ON)
option (KMC_ADD_DEVELOP "Add KMC developer tools" OFF)

add_subdirectory (kmc_api)
add_subdirectory (kmc_dump)
add_subdirectory (kmc_dump_sample)
add_subdirectory (kmer_counter)
