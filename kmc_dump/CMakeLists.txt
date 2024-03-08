# -*- cmake -*-

project(KMC_DUMP
        LANGUAGES CXX)

add_executable(kmc_dump kmc_dump.cpp nc_utils.cpp)
target_link_libraries(kmc_dump kmc_api)

install(TARGETS kmc_dump
        DESTINATION bin
        COMPONENT kmc)
