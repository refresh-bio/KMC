/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/

#ifndef _CHECK_KMER_H
#define _CHECK_KMER_H

#include "config.h"
#include "kmer.h"
#include "../kmc_api/mmer.h"
#include <string>

template<unsigned SIZE>
class CKmerCheck
{
	CConfig& config;
	const CKmerFileHeader& header;
	const CInputDesc& input_desc;
	CKmer<SIZE> mask;
	uint64 max_prefix;
	uint32 record_size;
	FILE* prefix_file = nullptr, *suffix_file = nullptr;
	
	std::unique_ptr<uchar[]> rec;

	uint32 check_in_suffix_file(uint64 lower, uint64 upper, CKmer<SIZE>& kmer_suffix)
	{
		
		if (!upper)
		{			
			return 0;
		}
		upper--;

		if (upper < lower)	
			return 0;
		
		

		

		uint32 cutoff_range = input_desc.cutoff_max - input_desc.cutoff_min;

		auto read_at = [](FILE* file, uint64 pos, uchar* &tmp, uint32 record_size, uint32 counter_size) -> CKmer<SIZE>
		{
			my_fseek(file, pos, SEEK_SET);
			fread(tmp, 1, record_size, file);
			CKmer<SIZE> res;
			res.load(tmp, record_size - counter_size);
			return res;
		};

		while (upper >= lower)
		{
			uint64 middle = (upper + lower) / 2;
			uchar* tmp = rec.get();
			auto middle_suffix = read_at(suffix_file, 4 + record_size * middle, tmp, record_size, header.counter_size);
			if (middle_suffix < kmer_suffix)
			{
				lower = middle + 1;
			}
			else if (kmer_suffix < middle_suffix)
			{
				upper = middle - 1;
			}
			else
			{
				uint32 counter = 0;
				for (uint32 i = 0; i < header.counter_size; ++i)
				{
					counter += (((uint32)*tmp++) << (i << 3));
				}
				if (counter - input_desc.cutoff_min < cutoff_range)
				{					
					return counter;
				}
				break;
			}
		}
		
		return 0;
	}

	void get_lower_upper(uint64 prefix, CKmer<SIZE>& kmer, uint64& lower, uint64& upper)
	{
		if (header.kmer_file_type == KmerFileType::KMC1)
		{
			uint64 pos = 4 + sizeof(uint64)*prefix;
			my_fseek(prefix_file, pos, SEEK_SET);

			fread(&lower, sizeof(uint64), 1, prefix_file);
			if (prefix == max_prefix)
				upper = header.total_kmers;
			else
				fread(&upper, sizeof(uint64), 1, prefix_file);
		}
		else if (header.kmer_file_type == KmerFileType::KMC2)
		{
			uint32 sig_len = header.signature_len;
			CMmer cur_mmr(sig_len);


			uint32 pos = header.kmer_len * 2 - 2;
			for (uint32 i = 0; i < sig_len; ++i)
			{
				cur_mmr.insert(kmer.get_2bits(pos));
				pos -= 2;
			}
			CMmer min_mmr(cur_mmr);
			for (uint32 i = sig_len; i < header.kmer_len; ++i)
			{
				cur_mmr.insert(kmer.get_2bits(pos));
				pos -= 2;

				if (cur_mmr < min_mmr)
					min_mmr = cur_mmr;
			}
			uint32 signature = min_mmr.get();

			uint64 map_size = ((1 << 2 * header.signature_len) + 1) * sizeof(uint32);
			my_fseek(prefix_file, 0ULL - (8 + header.header_offset + map_size) + signature * sizeof(uint32), SEEK_END);

			uint32 prefix_array = 0;
			uint32 prefix_arry_size = (1 << 2 * header.lut_prefix_len)*(uint32)sizeof(uint64);
			fread(&prefix_array, sizeof(uint32), 1, prefix_file);

			uint64 prefix_pos = 4 + prefix_arry_size*prefix_array + sizeof(uint64)*prefix;
			my_fseek(prefix_file, prefix_pos, SEEK_SET);

			fread(&lower, sizeof(uint64), 1, prefix_file);
			fread(&upper, sizeof(uint64), 1, prefix_file);
		}
		else
		{
			std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << "\t" << __LINE__ << "\n";
			exit(1);
		}
	}
public:
	CKmerCheck(const CKmerFileHeader& header, const CInputDesc& input_desc) :
		config(CConfig::GetInstance()),
		header(header),
		input_desc(input_desc)
	{
		mask.set_n_1((header.kmer_len - header.lut_prefix_len) * 2);

		std::string file_src = input_desc.file_src;
		prefix_file = fopen((file_src + ".kmc_pre").c_str(), "rb");
		if (!prefix_file)
		{
			std::cerr << "Error: cannot open file : " << (file_src + ".kmc_pre") << "\n";
			exit(1);
		}
		file_src = input_desc.file_src;
		suffix_file = fopen((file_src + ".kmc_suf").c_str(), "rb");
		if (!suffix_file)
		{
			std::cerr << "Error: cannot open file : " << (file_src + ".kmc_suf") << "\n";
			exit(1);
		}

		setvbuf(prefix_file, NULL, _IOFBF, (1ULL << 26));
		setvbuf(suffix_file, NULL, _IOFBF, (1ULL << 26));

		max_prefix = (1 << (2 * header.lut_prefix_len)) - 1;
		record_size = (header.kmer_len - header.lut_prefix_len) / 4 + header.counter_size;
		rec = std::make_unique<uchar[]>(record_size);
	}

	~CKmerCheck()
	{
		fclose(prefix_file);
		fclose(suffix_file);
	}

	uint32 CheckKmer(CKmer<SIZE> kmer)
	{
		uint64 prefix = kmer.remove_suffix((header.kmer_len - header.lut_prefix_len) * 2);
		
		uint64 lower = 0, upper = 0;		
		get_lower_upper(prefix, kmer, lower, upper);

		kmer.mask(mask);
		return check_in_suffix_file(lower, upper, kmer);
	}
	
	bool CheckKmer()
	{
		const std::string& kmer = config.check_params.kmer;
		if (kmer.length() != header.kmer_len)
		{
			std::cerr << "Error: invalid k-mer length\n";
			exit(1);
		}
		char codes[255];
		for (uint32 i = 0; i < 255; ++i)
			codes[i] = -1;
		codes['A'] = codes['a'] = 0;
		codes['C'] = codes['c'] = 1;
		codes['G'] = codes['g'] = 2;
		codes['T'] = codes['t'] = 3;

		
		CKmer<SIZE> _kmer;
		uint64 prefix = 0;
		_kmer.clear();

		for (uint32 i = 0; i < header.lut_prefix_len; ++i)
		{
			char d = codes[(uchar)kmer[i]];
			if (d < 0)
			{
				std::cerr << "Error: invalid k-mer format\n";
				exit(1);
			}
			prefix <<= 2;
			prefix += d;
			_kmer.SHL_insert_2bits(d);
		}
		for (uint32 i = header.lut_prefix_len; i<header.kmer_len; ++i)
		{
			char d = codes[(uchar)kmer[i]];
			if (d < 0)
			{
				std::cerr << "Error: invalid k-mer format\n";
				exit(1);
			}
			_kmer.SHL_insert_2bits(d);
		}
		
		uint64 lower = 0, upper = 0;
		get_lower_upper(prefix, _kmer, lower, upper);
		_kmer.mask(mask);
		uint32 counter = check_in_suffix_file(lower, upper, _kmer);
		std::cout << counter << "\n";
		return true;
	}
};


template<unsigned SIZE>
class CKmerCheckKFF
{
	CConfig& config;
	const CKmerFileHeader& header;
	const CInputDesc& input_desc;
	FILE* file = nullptr;

	std::unique_ptr<uchar[]> rec;

	bool check_in_file(uint64 offset, uint64 lower, uint64 upper, CKmer<SIZE>& kmer, uint32& counter, uint32 record_size, uint32 counter_size)
	{
		if (!upper)
		{
			return false;
		}
		upper--;

		if (upper < lower)
			return 0;

		uint32 cutoff_range = input_desc.cutoff_max - input_desc.cutoff_min;

		auto read_at = [](FILE* file, uint64 pos, uchar*& tmp, uint32 record_size, uint32 counter_size) -> CKmer<SIZE>
		{
			my_fseek(file, pos, SEEK_SET);
			fread(tmp, 1, record_size, file);
			CKmer<SIZE> res;
			res.load(tmp, record_size - counter_size);
			return res;
		};

		while (upper >= lower)
		{
			uint64 middle = (upper + lower) / 2;
			uchar* tmp = rec.get();
			auto middle_kmer = read_at(file, offset + record_size * middle, tmp, record_size, counter_size);
			if (middle_kmer < kmer)
			{
				lower = middle + 1;
			}
			else if (kmer < middle_kmer)
			{
				upper = middle - 1;
			}
			else
			{
				counter = 0;
				for (int32 i = counter_size - 1; i >= 0; --i)
				{
					counter += (((uint32)*tmp++) << (i << 3));
				}
				if (counter - input_desc.cutoff_min < cutoff_range)
				{
					return true;
				}
				break;
			}
		}

		return false;
	}
public:
	CKmerCheckKFF(const CKmerFileHeader& header, const CInputDesc& input_desc) :
		config(CConfig::GetInstance()),
		header(header),
		input_desc(input_desc)
	{
		std::string file_src = input_desc.file_src;

		file = fopen(file_src.c_str(), "rb");
		if (!file)
			file = fopen((file_src + ".kff").c_str(), "rb");

		if (!file)
		{
			std::cerr << "Error: cannot open file " << file_src << "\n";
			exit(1);
		}

		setvbuf(file, NULL, _IOFBF, (1ULL << 26));
	}

	~CKmerCheckKFF()
	{
		fclose(file);
	}

	bool CheckKmer()
	{
		const std::string& kmer = config.check_params.kmer;
		if (kmer.length() != header.kmer_len)
		{
			std::cerr << "Error: invalid k-mer length\n";
			exit(1);
		}

		char codes[255];
		for (uint32 i = 0; i < 255; ++i)
			codes[i] = -1;

		uint8_t encoding = header.GetEncoding();

		uint8_t A_representation = (encoding >> 6) & 3;
		uint8_t C_representation = (encoding >> 4) & 3;
		uint8_t G_representation = (encoding >> 2) & 3;
		uint8_t T_representation =  encoding       & 3;

		codes['A'] = codes['a'] = A_representation;
		codes['C'] = codes['c'] = C_representation;
		codes['G'] = codes['g'] = G_representation;
		codes['T'] = codes['t'] = T_representation;

		CKmer<SIZE> _kmer;

		_kmer.clear();

		for (uint32 i = 0; i < header.kmer_len; ++i)
		{
			char d = codes[(uchar)kmer[i]];
			if (d < 0)
			{
				std::cerr << "Error: invalid k-mer format\n";
				exit(1);
			}
			_kmer.SHL_insert_2bits(d);
		}

		uint32 counter;
		bool any = false;

		uint32 counter_size = header.counter_size;
		uint32 record_size = (header.kmer_len + 3) / 4 + counter_size;
		rec = std::make_unique<uchar[]>(record_size);
		uint32 prev_counter_size = counter_size;
		for (const auto& scope : header.kff_file_struct.scopes)
		{
			counter_size = scope.data_size;
			if(counter_size != prev_counter_size)
			{
				record_size = (header.kmer_len + 3) / 4 + counter_size;
				rec = std::make_unique<uchar[]>(record_size);
				prev_counter_size = counter_size;
			}

			for (const auto& section : scope.data_sections)
			{
				 if(check_in_file(section.data_start_pos, 0, section.nb_blocks, _kmer, counter, record_size, counter_size))
				 {
					 any = true;
					 std::cout << counter << "\n";
				 }
			}
		}
		return any;
	}
};

#endif