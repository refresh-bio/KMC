/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _KFF_RANDOM_ACCESS_H
#define _KFF_RANDOM_ACCESS_H
#include "config.h"
#include "kmer.h"
#include <vector>
#include <string>
#include "kff_db_reader.h"
#include "db_reader_factory.h"
#include "../kmc_api/kmc_file.h"
//TODO KFF: remember to use apropriate encoding! Important!
class CKffAndKMCRandomAccess : public CKMCFile
{
	uint32_t suffix_bits;
	uint64_t CalcLutPrefixLen(uint64_t tot_kmers, uint32_t kmer_len, uint32 counter_size)
	{
		uint32_t best_lut_prefix_len = 0;
		uint64_t best_mem_amount = 1ull << 62;

		for (uint64_t lut_prefix_len = 1; lut_prefix_len < 16; ++lut_prefix_len)
		{
			uint32 suffix_len;
			if (lut_prefix_len > kmer_len)
				suffix_len = 0;
			else
				suffix_len = kmer_len - lut_prefix_len;

			if (suffix_len % 4)
				continue;

			uint64 suf_mem = tot_kmers * (suffix_len / 4 + counter_size);
			uint64 lut_mem = (1ull << (2 * lut_prefix_len)) * sizeof(uint64);

			if (suf_mem + lut_mem < best_mem_amount)
			{
				best_lut_prefix_len = lut_prefix_len;
				best_mem_amount = suf_mem + lut_mem;
			}
		}
		return best_lut_prefix_len;
	}

	template<unsigned SIZE>
	void Store(CKmer<SIZE>& kmer, uint32_t counter, uchar* &ptr)
	{
		uint64_t prefix = kmer.remove_suffix(suffix_bits);
		++prefix_file_buf[prefix];
		kmer.store(ptr, sufix_size);
		for (uint32 j = 0; j < counter_size; ++j)
			*ptr++ = (counter >> (j * 8)) & 0xFF;
	}

public:
	template<unsigned SIZE>
	void OpenKFF(CKmerFileHeader& kmer_file_header, CInputDesc& input_desc)
	{
		auto& config = CConfig::GetInstance();
		uint32_t n_threads = config.input_desc.front().threads;
		total_kmers = 0;
		for (auto& scope : kmer_file_header.kff_file_struct.scopes)
			for (auto& section : scope.data_sections)
				total_kmers += section.nb_blocks;
		std::cerr << "Reading " << total_kmers << " k-mers from KFF file using " << n_threads << " threads\n";
		
		CKFFDbReader<SIZE>* kff_db_reader = new CKFFDbReader<SIZE>(kmer_file_header, input_desc, config.percent_progress, KmerDBOpenMode::sorted);

		kmer_length = kmer_file_header.kmer_len;
		counter_size = kmer_file_header.counter_size;
		lut_prefix_length = CalcLutPrefixLen(total_kmers, kmer_length, counter_size);
		prefix_file_buf_size = (1ull << (2 * lut_prefix_length)) + sizeof(uint64_t);
		prefix_file_buf = new uint64[prefix_file_buf_size];

		std::fill(prefix_file_buf, prefix_file_buf + prefix_file_buf_size, 0ull);

		sufix_size = (kmer_length - lut_prefix_length) / 4;

		suffix_bits = sufix_size * 8;
		
		sufix_rec_size = sufix_size + counter_size;

		sufix_file_buf = new uchar[sufix_rec_size * total_kmers];

		CBundle<SIZE> bundle(kff_db_reader);
		
		uchar* ptr = sufix_file_buf;
		while (!bundle.Finished())
		{			
			Store(bundle.TopKmer(), bundle.TopCounter(), ptr);
			bundle.Pop();
		}
		kmc_version = 0; //set as KMC1
		both_strands = kmer_file_header.kff_file_struct.both_strands;
		uint64_t sum{};
		for (uint64 prefix = 0; prefix < prefix_file_buf_size; ++prefix)
		{
			auto tmp = prefix_file_buf[prefix];
			prefix_file_buf[prefix] = sum;
			sum += tmp;
		}
		assert(prefix_file_buf[prefix_file_buf_size - 1] == total_kmers);

		original_min_count = min_count = kmer_file_header.min_count;
		original_max_count = max_count = kmer_file_header.max_count;
		is_opened = opened_for_RA;
	}	
};

#endif //_KFF_RANDOM_ACCESS_H


