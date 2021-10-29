/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#include "kmer_file_header.h"
#include "kff_info_reader.h"
#include <cstring>
#include <set>

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

CKmerFileHeader::CKmerFileHeader(std::string file_name)
{
	if (is_kff_file(file_name))
	{
		kmer_file_type = KmerFileType::KFF1;
		read_from_kff_file(file_name);
		
		return;
	}

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
	uint64 file_size = my_ftell(file);

	my_fseek(file, -8, SEEK_END);
	load_uint(file, header_offset);

	my_fseek(file, -12, SEEK_END);
	load_uint(file, db_version);

	kmer_file_type = db_version == 0x200 ? KmerFileType::KMC2 : KmerFileType::KMC1;

	my_fseek(file, 0LL - (header_offset + 8), SEEK_END);
	load_uint(file, kmer_len);
	load_uint(file, mode);
	load_uint(file, counter_size);
	load_uint(file, lut_prefix_len);
	if (kmer_file_type == KmerFileType::KMC2)
		load_uint(file, signature_len);
	load_uint(file, min_count);
	uint32_t max_count_lo;
	load_uint(file, max_count_lo);
	load_uint(file, total_kmers);
	uchar both_s_tmp;
	load_uint(file, both_s_tmp);
	both_strands = both_s_tmp == 1;
	both_strands = !both_strands;

	fseek(file, 3, SEEK_CUR);
	uint32_t max_count_hi;
	load_uint(file, max_count_hi);
	max_count = (((uint64_t)max_count_hi) << 32) + max_count_lo;
	fclose(file);

	if (kmer_file_type == KmerFileType::KMC2)
	{
		uint32 single_lut_size = (1ull << (2 * lut_prefix_len)) * sizeof(uint64);
		uint32 map_size = ((1 << 2 * signature_len) + 1) * sizeof(uint32);
		no_of_bins = (uint32)((file_size - sizeof(uint64) - 12 - header_offset - map_size) / single_lut_size);
	}
}

bool CKmerFileHeader::is_kff_file(std::string& fname)
{
	auto file = my_fopen(fname.c_str(), "rb");
	if (!file)
	{
		file = my_fopen((fname + ".kff").c_str(), "rb");
		if (file)
			fname += ".kff";
		else
			return false;
	}
	
	char marker_start[3];
	fread(marker_start, 1, 3, file);
	char marker_end[3];
	my_fseek(file, -3, SEEK_END);
	fread(marker_end, 1, 3, file);

	if (strncmp("KFF", marker_start, 3) || strncmp("KFF", marker_end, 3))	
		return false;

	return true;
}

void CKmerFileHeader::read_from_kff_file(const std::string& fname)
{
	CKFFInfoReader kff_reader(fname);
	kff_file_struct = kff_reader.GetKffFileStruct();

	std::set<uint64_t> k_values;
	this->counter_size = 0;
	if (kff_file_struct.scopes.empty())
	{
		std::cerr << "Error: no not-empty scope was found in KFF file. Make sure KFF file is fully indexed\n";
		exit(1);
	}
	for (auto& s : kff_file_struct.scopes)
	{
		k_values.insert(s.kmer_size);
		if (!s.ordered)
		{
			std::cerr << "Error: kmc_tools requires all KFF sections to be ordered\n";
			exit(1);
		}

		for (auto& section: s.data_sections)
			if (!(section.type == KFFDataSectionType::RAW))
			{
				std::cerr << "Error: currently kmc_tools supports only raw KFF sections\n";
				exit(1);
			}

		if (s.data_size > this->counter_size)
			this->counter_size = s.data_size;
	}

	this->min_count = GetFromFooterOrDefault("min_count", 1);
	this->max_count = GetFromFooterOrDefault("max_count", std::numeric_limits<uint32>::max());
	this->both_strands = kff_file_struct.both_strands;

	if (k_values.size() == 1)
		kmer_len = *k_values.begin();
	else
	{
		std::cerr << "Error: only KFF files with single k value are supported. This file contains following k values:";		
		for (auto k : k_values)
			std::cerr << " " << k;
		std::cerr << "\n";
		exit(1);		
	}
}
// ***** EOF