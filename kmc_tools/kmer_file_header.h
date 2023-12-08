/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.2.3
  Date   : 2023-12-08
*/

#ifndef _KMER_FILE_HEADER_H
#define _KMER_FILE_HEADER_H
#include "defs.h"
#include "kff_info_reader.h"
#include <string>
#include <iostream>

//************************************************************************************************************
// CKmerFileHeader - represents header of k-mer database.
//************************************************************************************************************

enum class KmerFileType { KMC1, KMC2, KFF1 };

struct CKmerFileHeader
{
	bool is_kff_file(std::string& fname);
	void read_from_kff_file(const std::string& fname);
	CKFFFileStruct kff_file_struct;

	uint64_t GetFromFooterOrDefault(const std::string& name, uint64_t default_value)
	{
		const auto& m = kff_file_struct.footer;
		auto r = m.find(name);
		if (r != m.end())
			return r->second;
		return default_value;
	}
public:
	uint32 kmer_len = 0;
	uint32 mode = 0;
	uint32 counter_size = 0;
	uint32 lut_prefix_len = 0;
	uint32 signature_len = 0; //only for kmc2
	uint32 min_count = 0;
	uint64 max_count = 0;
	uint64 total_kmers = 0;
	bool both_strands = true;
	uint32 db_version = 0;
	uint32 header_offset = 0;
	
	uint32 no_of_bins = 0; //only for kmc2
	KmerFileType kmer_file_type;
	//bool IsKMC2() const
	//{
	//	return db_version == 0x200;
	//}

	KmerFileType GetType() const
	{
		return kmer_file_type;
	}

	uint8_t GetEncoding() const
	{
		switch (kmer_file_type)
		{
			case KmerFileType::KMC1:
			case KmerFileType::KMC2:
				return 0b00011011;
			case KmerFileType::KFF1:
				return kff_file_struct.encoding;
			default:
			{
				std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << "\t" << __LINE__ << "\n";
				exit(1);
			}
		}
	}
	CKmerFileHeader(std::string file_name);

	

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