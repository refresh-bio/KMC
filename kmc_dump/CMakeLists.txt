set (KMC_DUMP_INCLUDES
	nc_utils.h
	stdafx.h
	targetver.h)
set (KMC_DUMP_SOURCES
	kmc_dump.cpp
	nc_utils.cpp
	stdafx.cpp)

add_executable (kmc_dump ${KMC_DUMP_INCLUDES} ${KMC_DUMP_SOURCES})
target_link_libraries (kmc_dump kmc_api)

