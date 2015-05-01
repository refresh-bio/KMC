set (KMC_DUMP_SAMPLE_SOURCES kmc_dump_sample.cpp)

add_executable (kmc_dump_sample ${KMC_DUMP_SAMPLE_SOURCES})
target_link_libraries (kmc_dump_sample kmc_api)
