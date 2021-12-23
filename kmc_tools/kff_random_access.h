/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.2.0
  Date   : 2021-12-23
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

class CKffAndKMCRandomAccess : protected CKMCFile
{
	uint32_t suffix_bits;
	bool need_to_encode_reads = false; //if KFF uses different encoding than ACGT -> 0b00011011 the symbols in reads needs to be encoded according to this encoding
	char enocde_reads_map[256];
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

		uint8_t encoding = kmer_file_header.GetEncoding();
		if (encoding != 0b00011011)
		{
			need_to_encode_reads = true;
			for (int i = 0; i < 256; ++i)
				enocde_reads_map[i] = i;

			uint8_t A_representation = (encoding >> 6) & 3;
			uint8_t C_representation = (encoding >> 4) & 3;
			uint8_t G_representation = (encoding >> 2) & 3;
			uint8_t T_representation =  encoding       & 3;

			enocde_reads_map['A'] = "ACGT"[A_representation];
			enocde_reads_map['a'] = "acgt"[A_representation];

			enocde_reads_map['C'] = "ACGT"[C_representation];
			enocde_reads_map['c'] = "acgt"[C_representation];

			enocde_reads_map['G'] = "ACGT"[G_representation];
			enocde_reads_map['g'] = "acgt"[G_representation];

			enocde_reads_map['T'] = "ACGT"[T_representation];
			enocde_reads_map['t'] = "acgt"[T_representation];
		}
	}

	bool GetCountersForRead(const std::string& read, std::vector<uint32>& counters)
	{
		if (need_to_encode_reads)
		{
			std::string encoded_read = read;
			for (auto& c : encoded_read)
				c = enocde_reads_map[(int)c];
			return CKMCFile::GetCountersForRead(encoded_read, counters);
		}
		return CKMCFile::GetCountersForRead(read, counters);
	}

	bool OpenKMC(const std::string& file_name)
	{
		return CKMCFile::OpenForRA(file_name);
	}

	bool SetMinCount(uint32 x)
	{
		return CKMCFile::SetMinCount(x);
	}

	bool SetMaxCount(uint32 x)
	{
		return CKMCFile::SetMaxCount(x);
	}
};

#endif //_KFF_RANDOM_ACCESS_H


