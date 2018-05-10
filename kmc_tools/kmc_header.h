/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _KMC_HEADER_H
#define _KMC_HEADER_H
#include "defs.h"
#include <string>
#include <iostream>

//************************************************************************************************************
// CKMC_header - represents header of KMC database.
//************************************************************************************************************
struct CKMC_header
{
public:
	uint32 kmer_len = 0;
	uint32 mode = 0;
	uint32 counter_size = 0;
	uint32 lut_prefix_len = 0;
	uint32 signature_len = 0; //only for kmc2
	uint32 min_count = 0;
	uint32 max_count = 0;
	uint64 total_kmers = 0;
	bool both_strands = true;
	uint32 db_version = 0;
	uint32 header_offset = 0;
	uint64 file_size = 0;

	uint32 no_of_bins = 0; //only for kmc2
	bool IsKMC2()const
	{
		return db_version == 0x200;
	}
	CKMC_header(std::string file_name);

private:
	template<typename T> void load_uint(FILE* file, T& res)
	{
		res = 0;
		for (uint32 i = 0; i < sizeof(T); ++i)
			res += (T)getc(file) << (i << 3);
	}
};

#endif


// ***** EOF