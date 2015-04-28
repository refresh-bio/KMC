set (KMC_DUMP_SAMPLE_INCLUDES
	stdafx.h
	targetver.h)
set (KMC_DUMP_SAMPLE_SOURCES
	kmc_dump_sample.cpp
	stdafx.cpp)

add_executable (kmc_dump_sample ${KMC_DUMP_SAMPLE_INCLUDES} ${KMC_DUMP_SAMPLE_SOURCES})
target_link_libraries (kmc_dump_sample kmc_api)
