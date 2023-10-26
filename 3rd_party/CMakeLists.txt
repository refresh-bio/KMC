# -*- cmake -*-
include(ExternalProject)

include(FetchContent)

FetchContent_Declare(cloudflare_zlib
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cloudflare")
FetchContent_GetProperties(cloudflare_zlib)
if (NOT cloudflare_zlib_POPULATED)
  FetchContent_Populate(cloudflare_zlib)
  add_subdirectory(${cloudflare_zlib_SOURCE_DIR} ${cloudflare_zlib_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
