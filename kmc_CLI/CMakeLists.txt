# -*- cmake -*-

project(KMC_CLI
        LANGUAGES CXX)

add_executable(kmc kmc.cpp)
target_link_libraries(kmc kmc_core)

install(TARGETS kmc
        DESTINATION bin
        COMPONENT kmc)
