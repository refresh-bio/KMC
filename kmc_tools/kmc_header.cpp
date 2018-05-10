/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include "kmc_header.h"
#include <cstring>

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

CKMC_header::CKMC_header(std::string file_name)
{
	file_name += ".kmc_pre";
	FILE* file = my_fopen(file_name.c_str(), "rb");
	if (!file)
	{
		std::cerr << "Error: Cannot open file " << file_name << "\n";
		exit(1);
	}
	char marker[4];
	if (fread(marker, 1, 4, file) != 4)
	{
		std::cerr << "Error while reading start marker in " << file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCP", 4) != 0)
	{
		std::cerr << "Error: wrong start marker in " << file_name << "\n";
		exit(1);
	}

	my_fseek(file, -4, SEEK_END);
	if (fread(marker, 1, 4, file) != 4)
	{
		std::cerr << "Error while reading end marker in " << file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCP", 4) != 0)
	{
		std::cerr << "Error: wrong end marker in " << file_name << "\n";
		exit(1);
	}

	my_fseek(file, 0, SEEK_END);
	file_size = my_ftell(file);

	my_fseek(file, -8, SEEK_END);
	load_uint(file, header_offset);

	my_fseek(file, -12, SEEK_END);
	load_uint(file, db_version);

	my_fseek(file, 0LL - (header_offset + 8), SEEK_END);
	load_uint(file, kmer_len);
	load_uint(file, mode);
	load_uint(file, counter_size);
	load_uint(file, lut_prefix_len);
	if (IsKMC2())
		load_uint(file, signature_len);
	load_uint(file, min_count);
	load_uint(file, max_count);
	load_uint(file, total_kmers);
	uchar both_s_tmp;
	load_uint(file, both_s_tmp);
	both_strands = both_s_tmp == 1;
	both_strands = !both_strands;

	fclose(file);

	if (IsKMC2())
	{
		uint32 single_lut_size = (1ull << (2 * lut_prefix_len)) * sizeof(uint64);
		uint32 map_size = ((1 << 2 * signature_len) + 1) * sizeof(uint32);
		no_of_bins = (uint32)((file_size - sizeof(uint64) - 12 - header_offset - map_size) / single_lut_size);
	}
}

// ***** EOF