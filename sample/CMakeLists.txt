set (KMC_DUMP_SAMPLE_SOURCES kmc_dump_sample.cpp)

add_executable (kmc_dump_sample ${KMC_DUMP_SAMPLE_SOURCES})
set_target_properties (kmc_dump_sample PROPERTIES VERSION ${KMC_VERSION})
set_target_properties (kmc_dump_sample PROPERTIES SOVERSION ${KMC_SOVERSION})
target_link_libraries (kmc_dump_sample kmc_api)
