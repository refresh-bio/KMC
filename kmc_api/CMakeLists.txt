set (KMC_API_INCLUDES
	kmc_file.h
	kmer_api.h
	kmer_defs.h
	mmer.h
	stdafx.h)
set (KMC_API_SOURCES
	kmc_file.cpp
	kmer_api.cpp
	mmer.cpp)
add_library (kmc_api ${KMC_API_INCLUDES} ${KMC_API_SOURCES})
