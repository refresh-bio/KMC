/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/
#ifndef _KB_COMPLETER_H
#define _KB_COMPLETER_H

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include "radix.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <stdio.h>
#include "small_k_buf.h"

//************************************************************************************************************
// CKmerBinCompleter - complete the sorted bins and store in a file
//************************************************************************************************************
class CKmerBinCompleter {
	CMemoryMonitor *mm;
	string file_name, kmer_file_name, lut_file_name;
	CKmerQueue *kq;
	CBinDesc *bd;
	CSignatureMapper *s_mapper;
	uint32 *sig_map;
	uint64 _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total;
	uint64 n_recs;

	FILE *out_kmer, *out_lut;
	uint32 lut_pos;
	uint32 sig_map_size;
	uint64 counter_size;

	CMemoryBins *memory_bins;

	bool use_strict_mem;
	
	CBigBinKmerPartQueue* bbkpq;
	CMemoryPool *sm_pmm_merger_lut, *sm_pmm_merger_suff;

	uint32 lut_prefix_len;
	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total;
	uint32 kmer_t_size;
	uint32 cutoff_min, cutoff_max;
	uint32 counter_max;
	int32 kmer_len;
	int32 signature_len;	
	bool both_strands;
	bool without_output;
	bool store_uint(FILE *out, uint64 x, uint32 size);	

public:
	CKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues);
	~CKmerBinCompleter();

	void ProcessBinsFirstStage();
	void ProcessBinsSecondStage();
	void GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total);
	void InitStage2(CKMCParams& Params, CKMCQueues& Queues);
};


//************************************************************************************************************
// CWKmerBinCompleter - wrapper for multithreading purposes
//************************************************************************************************************
class CWKmerBinCompleter {
	CKmerBinCompleter *kbc;

public:
	CWKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues);
	~CWKmerBinCompleter();

	void operator()(bool first_stage);

	void GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total);
	void InitStage2(CKMCParams& Params, CKMCQueues& Queues);
};


//************************************************************************************************************
// SmallKCompleter - completer for small k optimization
//************************************************************************************************************
//TODO: remove quake move implementation to cpp file
class CSmallKCompleter
{
	CMemoryPool *pmm_small_k_completer;
	uint64 n_unique, n_cutoff_min, n_cutoff_max;
	uint32 lut_prefix_len;
	int64 cutoff_max, counter_max;
	uint32 cutoff_min;
	uint32 kmer_len;
	int64 mem_tot_small_k_completer;
	std::string output_file_name;
	bool both_strands;	
	bool without_output;

	inline bool store_uint(FILE *out, uint64 x, uint32 size);
public:
	inline CSmallKCompleter(CKMCParams& Params, CKMCQueues& Queues);

	template<typename COUNTER_TYPE>
	bool Complete(CSmallKBuf<COUNTER_TYPE> results);
	inline void GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max);

};

CSmallKCompleter::CSmallKCompleter(CKMCParams& Params, CKMCQueues& Queues)
{
	pmm_small_k_completer = Queues.pmm_small_k_completer;
	n_unique = n_cutoff_min = n_cutoff_max = 0;
	lut_prefix_len = Params.lut_prefix_len;
	cutoff_max = Params.cutoff_max;
	cutoff_min = Params.cutoff_min;
	counter_max = Params.counter_max;
	both_strands = Params.both_strands;
	kmer_len = (uint32)Params.kmer_len;
	without_output = Params.without_output;

	mem_tot_small_k_completer = Params.mem_tot_small_k_completer;
	output_file_name = Params.output_file_name;
}

bool CSmallKCompleter::store_uint(FILE *out, uint64 x, uint32 size)
{
	for (uint32 i = 0; i < size; ++i)
		putc((x >> (i * 8)) & 0xFF, out);

	return true;
}

template<typename COUNTER_TYPE>
bool CSmallKCompleter::Complete(CSmallKBuf<COUNTER_TYPE> result)
{
	uchar* raw_buffer;
	uint64 counter_size = 0;

	counter_size = min(BYTE_LOG_ULL((uint64)cutoff_max), BYTE_LOG_ULL((uint64)counter_max));
	uint64 kmer_suf_bytes = (kmer_len - lut_prefix_len) / 4;


	pmm_small_k_completer->reserve(raw_buffer);
	uint32 lut_recs = (1 << 2 * lut_prefix_len);
	uint32 lut_buf_recs = (uint32)(MIN(lut_recs * sizeof(uint64), (uint64)mem_tot_small_k_completer / 2) / sizeof(uint64));
	uint32 lut_buf_pos = 0;
	uint32 suf_size = (uint32)(mem_tot_small_k_completer - lut_buf_recs * sizeof(uint64));
	uint32 suf_recs = (uint32)(suf_size / (counter_size + kmer_suf_bytes) * (counter_size + kmer_suf_bytes));
	uint32 suf_pos = 0;

	uint64* lut = (uint64*)raw_buffer;
	uchar* suf = raw_buffer + lut_buf_recs * sizeof(uint64);

	FILE* suf_file = nullptr;
	FILE* pre_file = nullptr;

	string pre_file_name = output_file_name + ".kmc_pre";
	string suf_file_name = output_file_name + ".kmc_suf";

	if (!this->without_output)
	{
		pre_file = fopen(pre_file_name.c_str(), "wb");
		if (!pre_file)
		{
			cerr << "Error: Cannot create " << pre_file_name << "\n";
			exit(1);
			return false;
		}
		suf_file = fopen(suf_file_name.c_str(), "wb");

		if (!suf_file)
		{
			cerr << "Error: Cannot create " << suf_file_name << "\n";
			fclose(pre_file);
			exit(1);
			return false;
		}
	}
	char s_kmc_pre[] = "KMCP";
	char s_kmc_suf[] = "KMCS";

	if (!this->without_output)
	{
		// Markers at the beginning
		fwrite(s_kmc_pre, 1, 4, pre_file);
		fwrite(s_kmc_suf, 1, 4, suf_file);
	}

	CKmer<1> kmer;
	
	uint64 prev_prefix = 0, prefix;

	lut[lut_buf_pos++] = 0;
	uint64 kmer_no = 0;
	for (kmer.data = 0; kmer.data < (1ull << 2 * kmer_len); ++kmer.data) 
	{
		if (!this->without_output)
		{
			prefix = kmer.remove_suffix(2 * (kmer_len - lut_prefix_len));

			if (prefix != prev_prefix) //new prefix
			{
				prev_prefix = prefix;
				lut[lut_buf_pos++] = kmer_no;
				if (lut_buf_pos >= lut_buf_recs)
				{
					fwrite(lut, sizeof(uint64), lut_buf_pos, pre_file);
					lut_buf_pos = 0;
				}
			}
		}

		if (result.buf[kmer.data]) //k-mer exists
		{
			++n_unique;

			if (result.buf[kmer.data] < cutoff_min)
				++n_cutoff_min;
			else if (result.buf[kmer.data] > (uint64)cutoff_max)
				++n_cutoff_max;
			else
			{
				++kmer_no;
				if (!this->without_output)
				{
					if (result.buf[kmer.data] > (uint64)counter_max)
						result.buf[kmer.data] = (COUNTER_TYPE)counter_max;

					for (int32 j = (int32)kmer_suf_bytes - 1; j >= 0; --j)
						suf[suf_pos++] = kmer.get_byte(j);

					result.Store(kmer.data, suf, suf_pos, counter_size);

					if (suf_pos >= suf_recs * (kmer_suf_bytes + counter_size))
					{
						fwrite(suf, 1, suf_pos, suf_file);
						suf_pos = 0;
					}
				}
			}
		}
	}

	if (!this->without_output)
	{
		fwrite(lut, sizeof(uint64), lut_buf_pos, pre_file);
		fwrite(suf, 1, suf_pos, suf_file);
	
		uint32 offset = 0;

		store_uint(pre_file, kmer_len, 4);					offset += 4;
		store_uint(pre_file, (uint32)0, 4);					offset += 4;	// mode: 0 (counting), 1 (Quake-compatibile counting, not supported anymore)
		store_uint(pre_file, counter_size, 4);				offset += 4;
		store_uint(pre_file, lut_prefix_len, 4);			offset += 4;
		store_uint(pre_file, cutoff_min, 4);				offset += 4;
		store_uint(pre_file, cutoff_max, 4);				offset += 4;
		store_uint(pre_file, n_unique - n_cutoff_min - n_cutoff_max, 8);		offset += 8;


		store_uint(pre_file, both_strands ? 0 : 1, 1);			offset++;

		store_uint(pre_file, 0, 1);			offset++;
		store_uint(pre_file, 0, 1);			offset++;
		store_uint(pre_file, 0, 1);			offset++;

		store_uint(pre_file, cutoff_max >> 32, 4);				offset += 4;
		// Space for future use
		for (int32 i = 0; i < 20; ++i)
		{
			store_uint(pre_file, 0, 1);
			offset++;
		}

		store_uint(pre_file, 0x0, 4); //KMC 1.x format
		offset += 4;

		store_uint(pre_file, offset, 4);

		// Markers at the end
		fwrite(s_kmc_pre, 1, 4, pre_file);
		fwrite(s_kmc_suf, 1, 4, suf_file);
		fclose(pre_file);
		fclose(suf_file);
	}
	pmm_small_k_completer->free(raw_buffer);
	

	return true;
}

void CSmallKCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max)
{
	_n_unique = n_unique;
	_n_cutoff_min = n_cutoff_min;
	_n_cutoff_max = n_cutoff_max;
}
#endif

// ***** EOF
