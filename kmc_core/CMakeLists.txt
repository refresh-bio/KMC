# -*- cmake -*-

project(KMC_CORE
        LANGUAGES C CXX)

# determine & normalize architecture
set(_system_processor ${CMAKE_SYSTEM_PROCESSOR})
if (APPLE AND NOT "${CMAKE_OSX_ARCHITECTURES}" STREQUAL "")
  set(_system_processor ${CMAKE_OSX_ARCHITECTURES})
endif()
if (${_system_processor} STREQUAL "arm64")
  set(ARCH "arm64")
elseif (${_system_processor} STREQUAL "aarch64")
  set(ARCH "arm64")
elseif (${_system_processor} STREQUAL "x86_64")
  set(ARCH "x86_64")
elseif (${_system_processor} STREQUAL "AMD64")
  set(ARCH "x86_64")
else()
  message(FATAL_ERROR "unknown architecture " ${_system_processor})
endif()

if (${ARCH} STREQUAL "arm64")
  set(raduls_sources raduls_neon.cpp)
elseif (${ARCH} STREQUAL "x86_64")
  set(raduls_sources
    raduls_sse2.cpp
    raduls_sse41.cpp
    raduls_avx2.cpp
    raduls_avx.cpp)
  set_source_files_properties(raduls_sse2.cpp PROPERTIES COMPILE_OPTIONS "-msse2")
  set_source_files_properties(raduls_sse41.cpp PROPERTIES COMPILE_OPTIONS "-msse4.1")
  set_source_files_properties(raduls_avx.cpp PROPERTIES COMPILE_OPTIONS "-mavx")
  set_source_files_properties(raduls_avx2.cpp PROPERTIES COMPILE_OPTIONS "-mavx2")
endif()

add_library(kff_core OBJECT
  kff_writer.cpp)

add_library(kmc_core STATIC
  mem_disk_file.cpp rev_byte.cpp bkb_writer.cpp 
  cpu_info.cpp bkb_reader.cpp fastq_reader.cpp 
  timer.cpp develop.cpp kb_completer.cpp 
  kb_storer.cpp kmer.cpp splitter.cpp kb_collector.cpp kmc_runner.cpp
  ${raduls_sources}
)
target_link_libraries(kmc_core kmc_api kff_core zlib)
