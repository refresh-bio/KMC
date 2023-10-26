# -*- cmake -*-

# Hard prerequisites
cmake_minimum_required(VERSION 3.16)

project(KMC
        VERSION 3.2.2
        LANGUAGES C CXX)

# FIXME:
# set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Sanity check our source directory to make sure that we are not trying to
# generate an in-tree build
if ((CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR) OR
    (CMAKE_SOURCE_DIR STREQUAL "${CMAKE_BINARY_DIR}/src"))
  message(FATAL_ERROR "In-source builds are not allowed.
Please create a directory and run cmake from there, passing the path
to this source directory as the last argument.
This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
Please delete them.")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FIXME: Add option to build pybindings

# External dependencies
add_subdirectory(3rd_party)

# Default build configuration
set(KMC_DEFAULT_BUILD_TYPE "RelWithDebInfo" CACHE STRING "KMC default build type")
if (NOT CMAKE_BUILD_TYPE)
  message("Setting default build configuration: ${KMC_DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE "${KMC_DEFAULT_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithAsserts RelWithDebInfo."
    FORCE)
endif()
  
# Default compiler flags
if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  message("Making Debug Configuration...")

  add_compile_options(-g3)
  add_definitions(-D_GLIBCXX_DEBUG)
else()
  message("Making Release Configuration...")

  if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_compile_options(-g1)
  else()
    add_compile_options(-g0)
  endif()

  add_compile_options(-O3)
  if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithAsserts" OR
      ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_definitions(-UNDEBUG)
  else()
    add_definitions(-DNDEBUG)
  endif()
endif()
add_compile_options(-Wall -fsigned-char)

add_subdirectory(kmc_api)
add_subdirectory(kmc_core)
add_subdirectory(kmc_CLI)
add_subdirectory(kmc_dump)
add_subdirectory(kmc_tools)
