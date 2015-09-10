/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _KB_SORTER_H
#define _KB_SORTER_H

#define DEBUGG_INFO

#include "defs.h"
#include "prob_qual.h"
#include "params.h"
#include "kmer.h"
#include "radix.h"
#include "s_mapper.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <vector>
#include <stdio.h>
#include <cstddef>

#include "kxmer_set.h"
#include "rev_byte.h"


//************************************************************************************************************
template <typename KMER_T, unsigned SIZE> class CKmerBinSorter_Impl;

//************************************************************************************************************
// CKmerBinSorter - sorting of k-mers in a bin
//************************************************************************************************************
template <typename KMER_T, unsigned SIZE> class CKmerBinSorter {
	mutable mutex expander_mtx;
	uint64 input_pos;
	CMemoryMonitor *mm;
	CBinDesc *bd;
	CBinQueue *bq;
	CKmerQueue *kq;
	CMemoryPool *pmm_prob, *pmm_radix_buf, *pmm_expand;
	CMemoryBins *memory_bins;	

	CKXmerSet<KMER_T, SIZE> kxmer_set;	

	int32 n_bins;
	int32 bin_id;
	
	uchar *data;
	uint64 size;
	uint64 n_rec;
	uint64 n_plus_x_recs;
	string desc;
	uint32 buffer_size;
	uint32 kmer_len;
	uint32 max_x;

	uint64 sum_n_rec, sum_n_plus_x_rec;

	int n_omp_threads;

	bool both_strands;
	bool use_quake;
	CSignatureMapper* s_mapper;

	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total;
	int32 cutoff_min, cutoff_max;
	int32 lut_prefix_len;
	int32 counter_max;

	KMER_T *buffer_input, *buffer_tmp, *buffer;
	uint32 *kxmer_counters;

	void Sort();

	friend class CKmerBinSorter_Impl<KMER_T, SIZE>;

public:
	static uint32 PROB_BUF_SIZE;
	CKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, int thread_no);
	~CKmerBinSorter();


	void GetDebugStats(uint64& _sum_n_recs, uint64& _sum_n_plus_x_recs)
	{
		_sum_n_recs = sum_n_rec;
		_sum_n_plus_x_recs = sum_n_plus_x_rec;		
	}

	void set_omp_threads(uint32 _n_omp_threads) { n_omp_threads = _n_omp_threads; }
	void ProcessBins();
};

template <typename KMER_T, unsigned SIZE> uint32 CKmerBinSorter<KMER_T, SIZE>::PROB_BUF_SIZE = 1 << 14;


//************************************************************************************************************
// CKmerBinSorter_Impl - implementation of k-mer type- and size-dependent functions
//************************************************************************************************************
template <typename KMER_T, unsigned SIZE> class CKmerBinSorter_Impl {
public:
	static void Compact(CKmerBinSorter<KMER_T, SIZE> &ptr);
	static void Expand(CKmerBinSorter<KMER_T, SIZE> &ptr, uint64 tmp_size);
	static void ComapctKXmers(CKmerBinSorter<KMER_T, SIZE> &ptr, uint64& compacted_count);
};

template <unsigned SIZE> class CKmerBinSorter_Impl<CKmer<SIZE>, SIZE> {
	static uint64 FindFirstSymbOccur(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb);
	static void InitKXMerSet(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth);
	static void CompactKxmers(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr);
	static void PreCompactKxmers(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64& compacted_count);
	static void CompactKmers(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr);
	static void ExpandKxmersAll(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size);
	static void ExpandKxmersBoth(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size);
	static void ExpandKmersAll(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size);
	static void ExpandKmersBoth(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size);
	static void GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p);
	static void FromChildThread(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, CKmer<SIZE>* thread_buffer, uint64 size);
	static void ExpandKxmerBothParaller(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 start_pos, uint64 end_pos);
	
public:
	static void Compact(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr);
	static void Expand(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 tmp_size);
};

template <unsigned SIZE> class CKmerBinSorter_Impl<CKmerQuake<SIZE>, SIZE> {
public:
	static void Compact(CKmerBinSorter<CKmerQuake<SIZE>, SIZE> &ptr);
	static void Expand(CKmerBinSorter<CKmerQuake<SIZE>, SIZE> &ptr, uint64 tmp_size);
};

//************************************************************************************************************
// CKmerBinSorter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Assign queues and monitors
template <typename KMER_T, unsigned SIZE> CKmerBinSorter<KMER_T, SIZE>::CKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, int thread_no) : kxmer_set(Params.kmer_len)
{
	both_strands = Params.both_strands;
	mm = Queues.mm;
	n_bins = Params.n_bins;
	bd = Queues.bd;
	bq = Queues.bq;
	kq = Queues.kq;


	s_mapper = Queues.s_mapper;

	pmm_radix_buf = Queues.pmm_radix_buf;
	pmm_prob = Queues.pmm_prob;
	pmm_expand = Queues.pmm_expand;	
	
	memory_bins = Queues.memory_bins;

	cutoff_min = Params.cutoff_min;
	cutoff_max = (int32)Params.cutoff_max;
	counter_max = (int32)Params.counter_max;
	max_x = Params.max_x;
	use_quake = Params.use_quake;
	
	lut_prefix_len = Params.lut_prefix_len;

	n_omp_threads = Params.n_omp_threads[thread_no];

	sum_n_rec = sum_n_plus_x_rec = 0;
}

//----------------------------------------------------------------------------------
template <typename KMER_T, unsigned SIZE> CKmerBinSorter<KMER_T, SIZE>::~CKmerBinSorter()
{

}

//----------------------------------------------------------------------------------
// Process the bins
template <typename KMER_T, unsigned SIZE> void CKmerBinSorter<KMER_T, SIZE>::ProcessBins()
{
	uint64 tmp_size;
	uint64 tmp_n_rec;
	CMemDiskFile *file;
	
	SetMemcpyCacheLimit(8);

	// Process bins
	while (!bq->completed())
	{
		// Gat bin data description to sort
		if (!bq->pop(bin_id, data, size, n_rec))
		{
			continue;
		}

		// Get bin data
		bd->read(bin_id, file, desc, tmp_size, tmp_n_rec, n_plus_x_recs, buffer_size, kmer_len);


		// Uncompact the kmers - append truncate prefixes
		//Expand(tmp_size);
		CKmerBinSorter_Impl<KMER_T, SIZE>::Expand(*this, tmp_size);
		memory_bins->free(bin_id, CMemoryBins::mba_input_file);

		// Perfor sorting of kmers in a bin
		Sort();
		// Compact the same kmers (occurring at neighbour positions now)
		CKmerBinSorter_Impl<KMER_T, SIZE>::Compact(*this);
	}

	// Mark all the kmers are already processed
	kq->mark_completed();
}

template <unsigned SIZE> inline void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p)
{
	symb = (data_p[pos] >> byte_shift) & 3;
	if (byte_shift == 0)
	{
		++pos;
		byte_shift = 6;
	}
	else
		byte_shift -= 2;
}

template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::ExpandKmersAll(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size)
{
	uint64 pos = 0;
	ptr.input_pos = 0;
	CKmer<SIZE> kmer;
	uint32 kmer_bytes = (ptr.kmer_len + 3) / 4;

	CKmer<SIZE> kmer_mask;
	kmer_mask.set_n_1(ptr.kmer_len * 2);
	uchar *data_p = ptr.data;
	uchar additional_symbols;
	uint32 kmer_shr = SIZE * 32 - ptr.kmer_len;
	while (pos < tmp_size)
	{
		kmer.clear();
		additional_symbols = data_p[pos++];		
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1; i < kmer_bytes; ++i, --kmer_pos)
		{
			kmer.set_byte(kmer_pos, data_p[pos + i]);
		}
		pos += kmer_bytes;
		uchar byte_shift = 6 - (ptr.kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		ptr.buffer_input[ptr.input_pos++].set(kmer);
		for (int i = 0; i < additional_symbols; ++i)
		{
			uchar symb = (data_p[pos] >> byte_shift) & 3;
			if (byte_shift == 0)
			{
				++pos;
				byte_shift = 6;
			}
			else
				byte_shift -= 2;
			kmer.SHL_insert_2bits(symb);
			kmer.mask(kmer_mask);
			ptr.buffer_input[ptr.input_pos++].set(kmer);
		}
		if (byte_shift != 6)
			++pos;
	}
}
template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::ExpandKmersBoth(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size)
{
	uint64 pos = 0;
	CKmer<SIZE> kmer;
	CKmer<SIZE> rev_kmer;
	CKmer<SIZE> kmer_can;

	uint32 kmer_bytes = (ptr.kmer_len + 3) / 4;
	uint32 kmer_len_shift = (ptr.kmer_len - 1) * 2;
	CKmer<SIZE> kmer_mask;
	kmer_mask.set_n_1(ptr.kmer_len * 2);
	uchar *data_p = ptr.data;
	ptr.input_pos = 0;
	uint32 kmer_shr = SIZE * 32 - ptr.kmer_len;

	uchar additional_symbols;

	uchar symb;
	while (pos < tmp_size)
	{
		kmer.clear();
		rev_kmer.clear();
		additional_symbols = data_p[pos++];

		//building kmer
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1, kmer_rev_pos = 0; i < kmer_bytes; ++i, --kmer_pos, ++kmer_rev_pos)
		{
			kmer.set_byte(kmer_pos, data_p[pos + i]);
			rev_kmer.set_byte(kmer_rev_pos, CRev_byte::lut[data_p[pos + i]]);
		}
		pos += kmer_bytes;
		uchar byte_shift = 6 - (ptr.kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		rev_kmer.mask(kmer_mask);

		kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
		ptr.buffer_input[ptr.input_pos++].set(kmer_can);

		for (int i = 0; i < additional_symbols; ++i)
		{
			symb = (data_p[pos] >> byte_shift) & 3;
			if (byte_shift == 0)
			{
				++pos;
				byte_shift = 6;
			}
			else
				byte_shift -= 2;
			kmer.SHL_insert_2bits(symb);
			kmer.mask(kmer_mask);
			rev_kmer.SHR_insert_2bits(3 - symb, kmer_len_shift);
			kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
			ptr.buffer_input[ptr.input_pos++].set(kmer_can);
		}
		if (byte_shift != 6)
			++pos;
	}
}

template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::FromChildThread(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, CKmer<SIZE>* thread_buffer, uint64 size)
{
	lock_guard<mutex> lcx(ptr.expander_mtx);
	A_memcpy(ptr.buffer_input + ptr.input_pos, thread_buffer, size * sizeof(CKmer<SIZE>));
	ptr.input_pos += size;
}

template<unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::ExpandKxmerBothParaller(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 start_pos, uint64 end_pos)
{
	uchar* _raw_buffer;	
	ptr.pmm_expand->reserve(_raw_buffer);
	CKmer<SIZE>* buffer = (CKmer<SIZE>*)_raw_buffer;

	CKmer<SIZE> kmer, rev_kmer, kmer_mask;
	CKmer<SIZE>  kxmer_mask;
	bool kmer_lower; //true if kmer is lower than its rev. comp
	uint32 x, additional_symbols;
	uchar symb;
	uint32 kmer_bytes = (ptr.kmer_len + 3) / 4;
	uint32 rev_shift = ptr.kmer_len * 2 - 2;
	uchar *data_p = ptr.data;
	kmer_mask.set_n_1(ptr.kmer_len * 2);
	uint32 kmer_shr = SIZE * 32 - ptr.kmer_len;

	kxmer_mask.set_n_1((ptr.kmer_len + ptr.max_x + 1) * 2);

	uint64 buffer_pos = 0;
	uint64 pos = start_pos;

	while (pos < end_pos)
	{
		kmer.clear();
		rev_kmer.clear();
		additional_symbols = data_p[pos++];

		//building kmer
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1, kmer_rev_pos = 0; i < kmer_bytes; ++i, --kmer_pos, ++kmer_rev_pos)
		{
			kmer.set_byte(kmer_pos, data_p[pos + i]);
			rev_kmer.set_byte(kmer_rev_pos, CRev_byte::lut[data_p[pos + i]]);
		}
		pos += kmer_bytes;
		uchar byte_shift = 6 - (ptr.kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		rev_kmer.mask(kmer_mask);

		kmer_lower = kmer < rev_kmer;
		x = 0;
		if (kmer_lower)
			buffer[buffer_pos].set(kmer);
		else
			buffer[buffer_pos].set(rev_kmer);

		uint32 symbols_left = additional_symbols;
		while (symbols_left)
		{
			GetNextSymb(symb, byte_shift, pos, data_p);
			kmer.SHL_insert_2bits(symb);
			kmer.mask(kmer_mask);
			rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
			--symbols_left;

			if (kmer_lower)
			{
				if (kmer < rev_kmer)
				{
					buffer[buffer_pos].SHL_insert_2bits(symb);
					++x;
					if (x == ptr.max_x)
					{
						if (!symbols_left)
							break;

						buffer[buffer_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
						if (buffer_pos >= EXPAND_BUFFER_RECS)
						{
							FromChildThread(ptr, buffer, buffer_pos);
							buffer_pos = 0;
						}
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;

						if (kmer_lower)
							buffer[buffer_pos].set(kmer);
						else
							buffer[buffer_pos].set(rev_kmer);
					}
				}
				else
				{
					buffer[buffer_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
					if (buffer_pos >= EXPAND_BUFFER_RECS)
					{
						FromChildThread(ptr, buffer, buffer_pos);
						buffer_pos = 0;
					}
					x = 0;

					kmer_lower = false;
					buffer[buffer_pos].set(rev_kmer);

				}
			}
			else
			{
				if (!(kmer < rev_kmer))
				{
					buffer[buffer_pos].set_2bits(3 - symb, ptr.kmer_len * 2 + x * 2);
					++x;
					if (x == ptr.max_x)
					{
						if (!symbols_left)
							break;

						buffer[buffer_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
						if (buffer_pos >= EXPAND_BUFFER_RECS)
						{
							FromChildThread(ptr, buffer, buffer_pos);
							buffer_pos = 0;
						}
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;

						if (kmer_lower)
							buffer[buffer_pos].set(kmer);
						else
							buffer[buffer_pos].set(rev_kmer);
					}
				}
				else
				{
					buffer[buffer_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
					if (buffer_pos >= EXPAND_BUFFER_RECS)
					{
						FromChildThread(ptr, buffer, buffer_pos);
						buffer_pos = 0;
					}
					x = 0;

					buffer[buffer_pos].set(kmer);
					kmer_lower = true;
				}
			}
			
		}
		buffer[buffer_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
		if (buffer_pos >= EXPAND_BUFFER_RECS)
		{
			FromChildThread(ptr, buffer, buffer_pos);
			buffer_pos = 0;
		}
		if (byte_shift != 6)
			++pos;
	}
	if (buffer_pos)
	{
		FromChildThread(ptr, buffer, buffer_pos);
		buffer_pos = 0;
	}
	ptr.pmm_expand->free(_raw_buffer);
}


template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::ExpandKxmersBoth(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size)
{	
	ptr.input_pos = 0;
	uint32 threads = ptr.n_omp_threads;

	uint64 bytes_per_thread = (tmp_size + threads - 1) / threads;
	uint32 thread_no = 0;
	vector<thread> exp_threads;
	uint64 start = 0;
	uint64 pos = 0;
	for (; pos < tmp_size; pos += 1 + (ptr.data[pos] + ptr.kmer_len + 3) / 4)
	{
		if ((thread_no + 1) * bytes_per_thread <= pos)
		{
			exp_threads.push_back(thread(ExpandKxmerBothParaller, std::ref(ptr), start, pos));
			start = pos;
			++thread_no;
		}
	}
	if (start < pos)
	{
		exp_threads.push_back(thread(ExpandKxmerBothParaller, std::ref(ptr), start, tmp_size));
	}

	for (auto& p : exp_threads)
		p.join();

	ptr.n_plus_x_recs = ptr.input_pos;// !!!!!!!!		
}

template<unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::ExpandKxmersAll(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size)
{
	ptr.input_pos = 0;
	uint64 pos = 0;
	CKmer<SIZE> kmer_mask;

	CKmer<SIZE> kxmer;
	CKmer<SIZE> kxmer_mask;
	kxmer_mask.set_n_1((ptr.kmer_len + ptr.max_x) * 2);
	uchar *data_p = ptr.data;

	kmer_mask.set_n_1(ptr.kmer_len * 2);
	while (pos < tmp_size)
	{
		kxmer.clear();
		uint32 additional_symbols = data_p[pos++];

		uchar symb;

		uint32 kmer_bytes = (ptr.kmer_len + 3) / 4;
		//building kmer
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1; i < kmer_bytes; ++i, --kmer_pos)
		{
			kxmer.set_byte(kmer_pos, data_p[pos + i]);
		}

		pos += kmer_bytes;
		uchar byte_shift = 6 - (ptr.kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;
		uint32 kmer_shr = SIZE * 32 - ptr.kmer_len;

		if (kmer_shr)
			kxmer.SHR(kmer_shr);

		kxmer.mask(kmer_mask);
		uint32 tmp = MIN(ptr.max_x, additional_symbols);

		for (uint32 i = 0; i < tmp; ++i)
		{
			GetNextSymb(symb, byte_shift, pos, data_p);
			kxmer.SHL_insert_2bits(symb);
		}
		kxmer.set_2bits(tmp, (ptr.kmer_len + ptr.max_x) * 2);

		ptr.buffer_input[ptr.input_pos++].set(kxmer);
		additional_symbols -= tmp;

		uint32 kxmers_count = additional_symbols / (ptr.max_x + 1);
		uint32 kxmer_rest = additional_symbols % (ptr.max_x + 1);

		for (uint32 j = 0; j < kxmers_count; ++j)
		{
			for (uint32 i = 0; i < ptr.max_x + 1; ++i)
			{
				GetNextSymb(symb, byte_shift, pos, data_p);
				kxmer.SHL_insert_2bits(symb);
			}

			kxmer.mask(kxmer_mask);

			kxmer.set_2bits(ptr.max_x, (ptr.kmer_len + ptr.max_x) * 2);

			ptr.buffer_input[ptr.input_pos++].set(kxmer);
		}
		if (kxmer_rest)
		{
			uint32 i = 0;
			GetNextSymb(symb, byte_shift, pos, data_p);
			kxmer.SHL_insert_2bits(symb);
			kxmer.mask(kmer_mask);
			--kxmer_rest;
			for (; i < kxmer_rest; ++i)
			{
				GetNextSymb(symb, byte_shift, pos, data_p);
				kxmer.SHL_insert_2bits(symb);
			}

			kxmer.set_2bits(kxmer_rest, (ptr.kmer_len + ptr.max_x) * 2);
			ptr.buffer_input[ptr.input_pos++].set(kxmer);
		}
		if (byte_shift != 6)
			++pos;
	}
}

//----------------------------------------------------------------------------------
// Uncompact the kmers
template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::Expand(CKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64 tmp_size)
{
	uchar *raw_buffer_input, *raw_buffer_tmp;

	ptr.memory_bins->reserve(ptr.bin_id, raw_buffer_input, CMemoryBins::mba_input_array);
	ptr.memory_bins->reserve(ptr.bin_id, raw_buffer_tmp, CMemoryBins::mba_tmp_array);

	ptr.buffer_input = (CKmer<SIZE> *) raw_buffer_input;
	ptr.buffer_tmp = (CKmer<SIZE> *) raw_buffer_tmp;

	if (ptr.max_x)
	{
		if (ptr.both_strands)
			ExpandKxmersBoth(ptr, tmp_size);
		else
			ExpandKxmersAll(ptr, tmp_size);
	}
	else
	{
		if (ptr.both_strands)
			ExpandKmersBoth(ptr, tmp_size);
		else
			ExpandKmersAll(ptr, tmp_size);
	}
}


//----
template <unsigned SIZE> void CKmerBinSorter_Impl<CKmerQuake<SIZE>, SIZE>::Expand(CKmerBinSorter<CKmerQuake<SIZE>, SIZE>& ptr, uint64 tmp_size)
{
	uchar *data_p = ptr.data;
	uchar *raw_buffer_input, *raw_buffer_tmp;

	ptr.memory_bins->reserve(ptr.bin_id, raw_buffer_input, CMemoryBins::mba_input_array);
	ptr.memory_bins->reserve(ptr.bin_id, raw_buffer_tmp, CMemoryBins::mba_tmp_array);


	ptr.buffer_input = (CKmerQuake<SIZE> *) raw_buffer_input;
	ptr.buffer_tmp = (CKmerQuake<SIZE> *) raw_buffer_tmp;
	CKmerQuake<SIZE> current_kmer;
	CKmerQuake<SIZE> kmer_rev;
	CKmerQuake<SIZE> kmer_can;
	kmer_rev.clear();
	uint32 kmer_len_shift = (ptr.kmer_len - 1) * 2;
	CKmerQuake<SIZE> kmer_mask;
	kmer_mask.set_n_1(ptr.kmer_len * 2);

	ptr.input_pos = 0;
	uint64 pos = 0;

	double *inv_probs;
	ptr.pmm_prob->reserve(inv_probs);
	double kmer_prob;
	uchar qual, symb;
	uint32 inv_probs_pos;
	if (ptr.both_strands)
		while (pos < tmp_size)
		{
			uchar additional_symbols = data_p[pos++];
			inv_probs_pos = 0;
			kmer_prob = 1.0;
			
			for (uint32 i = 0; i < ptr.kmer_len; ++i)
			{
				symb = (data_p[pos] >> 6) & 3;
				qual = data_p[pos++] & 63;

				inv_probs[inv_probs_pos++] = CProbQual::inv_prob_qual[qual];

				current_kmer.SHL_insert_2bits(symb);
				kmer_rev.SHR_insert_2bits(3 - symb, kmer_len_shift);
				kmer_prob *= CProbQual::prob_qual[qual];
			}
			current_kmer.mask(kmer_mask);			
			if (kmer_prob >= CProbQual::MIN_PROB_QUAL_VALUE)
			{
				kmer_can = current_kmer < kmer_rev ? current_kmer : kmer_rev;
				kmer_can.quality = (float)kmer_prob;
				ptr.buffer_input[ptr.input_pos++].set(kmer_can);
			}
			for (uint32 i = 0; i < additional_symbols; ++i)
			{
				symb = (data_p[pos] >> 6) & 3;
				qual = data_p[pos++] & 63;

				current_kmer.SHL_insert_2bits(symb);
				current_kmer.mask(kmer_mask);
				kmer_rev.SHR_insert_2bits(3 - symb, kmer_len_shift);
				
				kmer_prob *= CProbQual::prob_qual[qual] * inv_probs[inv_probs_pos - ptr.kmer_len];
				inv_probs[inv_probs_pos++] = CProbQual::inv_prob_qual[qual];
				if (kmer_prob >= CProbQual::MIN_PROB_QUAL_VALUE)
				{
					kmer_can = current_kmer < kmer_rev ? current_kmer : kmer_rev;
					kmer_can.quality = (float)kmer_prob;
					ptr.buffer_input[ptr.input_pos++].set(kmer_can);
				}
			}
		}
	else
		while (pos < tmp_size)
		{
			uchar additional_symbols = data_p[pos++];
			inv_probs_pos = 0;
			kmer_prob = 1.0;

			for (uint32 i = 0; i < ptr.kmer_len; ++i)
			{
				symb = (data_p[pos] >> 6) & 3;
				qual = data_p[pos++] & 63;

				inv_probs[inv_probs_pos++] = CProbQual::inv_prob_qual[qual];

				current_kmer.SHL_insert_2bits(symb);
				kmer_prob *= CProbQual::prob_qual[qual];
			}
			current_kmer.mask(kmer_mask);
			if (kmer_prob >= CProbQual::MIN_PROB_QUAL_VALUE)
			{
				current_kmer.quality = (float)kmer_prob;
				ptr.buffer_input[ptr.input_pos++].set(current_kmer);
			}
			for (uint32 i = 0; i < additional_symbols; ++i)
			{
				symb = (data_p[pos] >> 6) & 3;
				qual = data_p[pos++] & 63;

				current_kmer.SHL_insert_2bits(symb);
				current_kmer.mask(kmer_mask);

				kmer_prob *= CProbQual::prob_qual[qual] * inv_probs[inv_probs_pos - ptr.kmer_len];
				inv_probs[inv_probs_pos++] = CProbQual::inv_prob_qual[qual];
				if (kmer_prob >= CProbQual::MIN_PROB_QUAL_VALUE)
				{
					current_kmer.quality = (float)kmer_prob;
					ptr.buffer_input[ptr.input_pos++].set(current_kmer);
				}
			}
		}
	ptr.pmm_prob->free(inv_probs);
}


//----------------------------------------------------------------------------------
// Sort the kmers
template <typename KMER_T, unsigned SIZE> void CKmerBinSorter<KMER_T, SIZE>::Sort()
{
	uint32 rec_len;
	uint64 sort_rec;
	if (max_x && !use_quake)
	{
		sort_rec = n_plus_x_recs;
		rec_len = (kmer_len + max_x + 1 + 3) / 4;
	}
	else
	{
		sort_rec = n_rec;
		rec_len = (kmer_len + 3) / 4;
	}
	sum_n_plus_x_rec += n_plus_x_recs;
	sum_n_rec += n_rec;
	
	if (sizeof(KMER_T) == 8)
	{
		uint64 *_buffer_input = (uint64*)buffer_input;
		uint64 *_buffer_tmp = (uint64*)buffer_tmp;

		RadixSort_buffer(pmm_radix_buf, _buffer_input, _buffer_tmp, sort_rec, rec_len, n_omp_threads);

		if (rec_len % 2)
			buffer = (KMER_T*)_buffer_tmp;
		else
			buffer = (KMER_T*)_buffer_input;
	}
	else
	{
		uint32 *_buffer_input = (uint32*)buffer_input;
		uint32 *_buffer_tmp = (uint32*)buffer_tmp;

		RadixSort_uint8(_buffer_input, _buffer_tmp, sort_rec, sizeof(KMER_T), offsetof(KMER_T, data), SIZE*sizeof(typename KMER_T::data_t), rec_len, n_omp_threads);
		if (rec_len % 2)
			buffer = (KMER_T*)_buffer_tmp;
		else
			buffer = (KMER_T*)_buffer_input;
	}
}

//----------------------------------------------------------------------------------
//Binary search position of first occurence of symbol 'symb' in [start_pos,end_pos). Offset defines which symbol in k+x-mer is taken.
template <unsigned SIZE> uint64 CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::FindFirstSymbOccur(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb)
{
	uint32 kxmer_offset = (ptr.kmer_len + ptr.max_x - offset) * 2;
	uint64 middle_pos;
	uchar middle_symb;
	while (start_pos < end_pos)
	{
		middle_pos = (start_pos + end_pos) / 2;
		middle_symb = ptr.buffer[middle_pos].get_2bits(kxmer_offset);
		if (middle_symb < symb)
			start_pos = middle_pos + 1;
		else
			end_pos = middle_pos;
	}
	return end_pos;
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::InitKXMerSet(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth)
{
	if (start_pos == end_pos)
		return;
	uint32 shr = ptr.max_x + 1 - offset;
	ptr.kxmer_set.init_add(start_pos, end_pos, shr);

	--depth;
	if (depth > 0)
	{
		uint64 pos[5];
		pos[0] = start_pos;
		pos[4] = end_pos;
		for (uint32 i = 1; i < 4; ++i)
			pos[i] = FindFirstSymbOccur(ptr, pos[i - 1], end_pos, offset, i);
		for (uint32 i = 1; i < 5; ++i)
			InitKXMerSet(ptr, pos[i - 1], pos[i], offset + 1, depth);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::PreCompactKxmers(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64& compacted_count)
{
	compacted_count = 0;

	CKmer<SIZE> *act_kmer;
	act_kmer = &ptr.buffer[0];
	ptr.kxmer_counters[compacted_count] = 1;

	for (uint32 i = 1; i < ptr.n_plus_x_recs; ++i)
	{
		if (*act_kmer == ptr.buffer[i])
			++ptr.kxmer_counters[compacted_count];
		else
		{
			ptr.buffer[compacted_count++] = *act_kmer;
			ptr.kxmer_counters[compacted_count] = 1;
			act_kmer = &ptr.buffer[i];
		}
	}
	ptr.buffer[compacted_count++] = *act_kmer;
}

//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::CompactKxmers(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr)
{
	ptr.kxmer_set.clear();
	ptr.kxmer_set.set_buffer(ptr.buffer);
	ptr.n_unique = 0;
	ptr.n_cutoff_min = 0;
	ptr.n_cutoff_max = 0;
	ptr.n_total = 0;

	uint32 kmer_symbols = ptr.kmer_len - ptr.lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;
	uint64 lut_recs = 1 << (2 * ptr.lut_prefix_len);
	uint64 lut_size = lut_recs * sizeof(uint64);


	uchar *out_buffer = NULL;
	uchar *raw_lut = NULL;

	ptr.memory_bins->reserve(ptr.bin_id, out_buffer, CMemoryBins::mba_suffix);
	ptr.memory_bins->reserve(ptr.bin_id, raw_lut, CMemoryBins::mba_lut);

	uint64 *lut = (uint64*)raw_lut;
	fill_n(lut, lut_recs, 0);

	uint64 out_pos = 0;

	if (ptr.n_plus_x_recs)
	{
		uchar* raw_kxmer_counters = NULL;
		ptr.memory_bins->reserve(ptr.bin_id, raw_kxmer_counters, CMemoryBins::mba_kxmer_counters);
		ptr.kxmer_counters = (uint32*)raw_kxmer_counters;
		uint64 compacted_count;
		PreCompactKxmers(ptr, compacted_count);

		uint64 pos[5];//pos[symb] is first position where symb occur (at first position of k+x-mer) and pos[symb+1] is first position where symb is not starting symbol of k+x-mer
		pos[0] = 0;
		pos[4] = compacted_count;
		for (uint32 i = 1; i < 4; ++i)
			pos[i] = FindFirstSymbOccur(ptr, pos[i - 1], compacted_count, 0, i);
		for (uint32 i = 1; i < 5; ++i)
			InitKXMerSet(ptr, pos[i - 1], pos[i], ptr.max_x + 2 - i, i);

		
		uint64 counter_pos = 0;

		uint64 counter_size = min(BYTE_LOG(ptr.cutoff_max), BYTE_LOG(ptr.counter_max));

		CKmer<SIZE> kmer, next_kmer;
		kmer.clear();
		next_kmer.clear();
		CKmer<SIZE> kmer_mask;
		kmer_mask.set_n_1(ptr.kmer_len * 2);
		uint32 count;
		//first
		ptr.kxmer_set.get_min(counter_pos, kmer);
		count = ptr.kxmer_counters[counter_pos];
		//rest
		while (ptr.kxmer_set.get_min(counter_pos, next_kmer))
		{
			if (kmer == next_kmer)
				count += ptr.kxmer_counters[counter_pos];
			else
			{
				ptr.n_total += count;
				++ptr.n_unique;
				if (count < (uint32)ptr.cutoff_min)
					ptr.n_cutoff_min++;
				else if (count >(uint32)ptr.cutoff_max)
					ptr.n_cutoff_max++;
				else
				{
					lut[kmer.remove_suffix(2 * kmer_symbols)]++;
					if (count > (uint32)ptr.counter_max)
						count = ptr.counter_max;

					// Store compacted kmer

					for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
						out_buffer[out_pos++] = kmer.get_byte(j);
					for (int32 j = 0; j < (int32)counter_size; ++j)
						out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;
				}
				count = ptr.kxmer_counters[counter_pos];
				kmer = next_kmer;
			}
		}


		//last one
		++ptr.n_unique;
		ptr.n_total += count;
		if (count < (uint32)ptr.cutoff_min)
			ptr.n_cutoff_min++;
		else if (count >(uint32)ptr.cutoff_max)
			ptr.n_cutoff_max++;
		else
		{
			lut[kmer.remove_suffix(2 * kmer_symbols)]++;
			if (count > (uint32)ptr.counter_max)
				count = ptr.counter_max;

			// Store compacted kmer
			for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
				out_buffer[out_pos++] = kmer.get_byte(j);
			for (int32 j = 0; j < (int32)counter_size; ++j)
				out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;
		}


		ptr.memory_bins->free(ptr.bin_id, CMemoryBins::mba_kxmer_counters);
	}


	// Push the sorted and compacted kmer bin to a priority queue in a form ready to be stored to HDD
	ptr.kq->push(ptr.bin_id, out_buffer, out_pos, raw_lut, lut_size, ptr.n_unique, ptr.n_cutoff_min, ptr.n_cutoff_max, ptr.n_total);

	if (ptr.buffer_input)
	{
		ptr.memory_bins->free(ptr.bin_id, CMemoryBins::mba_input_array);
		ptr.memory_bins->free(ptr.bin_id, CMemoryBins::mba_tmp_array);
	}
	ptr.buffer = NULL;
}




//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::CompactKmers(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr)
{
	uint64 i;

	uint32 kmer_symbols = ptr.kmer_len - ptr.lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;
	uint64 lut_recs = 1 << (2 * (ptr.lut_prefix_len));
	uint64 lut_size = lut_recs * sizeof(uint64);

	uint64 counter_size = min(BYTE_LOG(ptr.cutoff_max), BYTE_LOG(ptr.counter_max));

	uchar *out_buffer;
	uchar *raw_lut;

	ptr.memory_bins->reserve(ptr.bin_id, out_buffer, CMemoryBins::mba_suffix);
	ptr.memory_bins->reserve(ptr.bin_id, raw_lut, CMemoryBins::mba_lut);
	uint64 *lut = (uint64*)raw_lut;
	fill_n(lut, lut_recs, 0);

	uint64 out_pos = 0;
	uint32 count;
	CKmer<SIZE> *act_kmer;
	

	ptr.n_unique = 0;
	ptr.n_cutoff_min = 0;
	ptr.n_cutoff_max = 0;
	ptr.n_total = 0;


	if (ptr.n_rec)			// non-empty bin
	{
		act_kmer = &ptr.buffer[0];
		count = 1;

		ptr.n_total = ptr.n_rec;

		for (i = 1; i < ptr.n_rec; ++i)
		{
			if (*act_kmer == ptr.buffer[i])
				count++;
			else
			{
				if (count < (uint32)ptr.cutoff_min)
				{
					act_kmer = &ptr.buffer[i];
					ptr.n_cutoff_min++;
					ptr.n_unique++;
					count = 1;
				}
				else if (count >(uint32) ptr.cutoff_max)
				{
					act_kmer = &ptr.buffer[i];
					ptr.n_cutoff_max++;
					ptr.n_unique++;
					count = 1;
				}
				else
				{
					if (count > (uint32)ptr.counter_max)
						count = ptr.counter_max;

					// Store compacted kmer
					for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
						out_buffer[out_pos++] = act_kmer->get_byte(j);
					for (int32 j = 0; j < (int32)counter_size; ++j)
						out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;

					lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;

					act_kmer = &ptr.buffer[i];
					count = 1;
					ptr.n_unique++;
				}
			}
		}

		if (count < (uint32)ptr.cutoff_min)
		{
			ptr.n_cutoff_min++;
		}
		else if (count >= (uint32)ptr.cutoff_max)
		{
			ptr.n_cutoff_max++;
		}
		else
		{
			if (count >(uint32) ptr.counter_max)
				count = ptr.counter_max;

			for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
				out_buffer[out_pos++] = act_kmer->get_byte(j);
			for (int32 j = 0; j < (int32)counter_size; ++j)
				out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;
			lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;
		}
		ptr.n_unique++;
	}

	// Push the sorted and compacted kmer bin to a priority queue in a form ready to be stored to HDD
	ptr.kq->push(ptr.bin_id, out_buffer, out_pos, raw_lut, lut_size, ptr.n_unique, ptr.n_cutoff_min, ptr.n_cutoff_max, ptr.n_total);

	if (ptr.buffer_input)
	{
		ptr.memory_bins->free(ptr.bin_id, CMemoryBins::mba_input_array);
		ptr.memory_bins->free(ptr.bin_id, CMemoryBins::mba_tmp_array);
	}
	ptr.buffer = NULL;
}




//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::Compact(CKmerBinSorter<CKmer<SIZE>, SIZE> &ptr)
{
	if (ptr.max_x)
		CompactKxmers(ptr);
	else
		CompactKmers(ptr);
}
//----------------------------------------------------------------------------------
// Compact the kmers - the same kmers (at neighbour positions now) are compated to a single kmer and counter
template <unsigned SIZE> void CKmerBinSorter_Impl<CKmerQuake<SIZE>, SIZE>::Compact(CKmerBinSorter<CKmerQuake<SIZE>, SIZE> &ptr)
{
		uint64 i;
	
		uint32 kmer_symbols = ptr.kmer_len - ptr.lut_prefix_len;
		uint64 kmer_bytes = kmer_symbols / 4;
		uint64 lut_recs = 1 << (2 * (ptr.lut_prefix_len));
		uint64 lut_size = lut_recs * sizeof(uint64);

		uchar *out_buffer;
		uchar *raw_lut;
		ptr.memory_bins->reserve(ptr.bin_id, out_buffer, CMemoryBins::mba_suffix);
		ptr.memory_bins->reserve(ptr.bin_id, raw_lut, CMemoryBins::mba_lut);
		uint64 *lut = (uint64*)raw_lut;
		fill_n(lut, lut_recs, 0);

		uint64 out_pos = 0;
		double count;
		CKmerQuake<SIZE> *act_kmer;
	
		ptr.n_unique     = 0;
		ptr.n_cutoff_min = 0;
		ptr.n_cutoff_max = 0;
		ptr.n_total      = 0;
	
		if(ptr.n_rec)			// non-empty bin
		{
			act_kmer = &ptr.buffer[0];
			count = (double)act_kmer->quality;
			ptr.n_total = ptr.n_rec;
			for(i = 1; i < ptr.n_rec; ++i)
			{
				if(*act_kmer == ptr.buffer[i])
					count += ptr.buffer[i].quality;
				else
				{
					if(count < (double) ptr.cutoff_min)
					{
						act_kmer = &ptr.buffer[i];
						++ptr.n_cutoff_min;
						++ptr.n_unique;
						count = act_kmer->quality;
					}
					else if(count > (double) ptr.cutoff_max)
					{
						act_kmer = &ptr.buffer[i];
						++ptr.n_cutoff_max;
						++ptr.n_unique;
						count = act_kmer->quality;
					}
					else
					{
						if(count > (double) ptr.counter_max)
							count = (double) ptr.counter_max;
	
						// Store compacted kmer
						for(int32 j = (int32) kmer_bytes-1; j >= 0; --j)
							out_buffer[out_pos++] = act_kmer->get_byte(j);
						uint32 tmp;
						float f_count = (float) count;
						memcpy(&tmp, &f_count, 4);
						for(int32 j = 0; j < 4; ++j)
							out_buffer[out_pos++] = (tmp >> (j * 8)) & 0xFF;

						lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;

						act_kmer = &ptr.buffer[i];
						count = act_kmer->quality;
						++ptr.n_unique;
					}
				}
			}
	
			if(count < (double) ptr.cutoff_min)
			{
				++ptr.n_cutoff_min;
			}
			else if(count > (double) ptr.cutoff_max)
			{
				++ptr.n_cutoff_max;
			}
			else
			{
				if(count > (double) ptr.counter_max)
					count = (double) ptr.counter_max;
	
				for(int32 j = (int32) kmer_bytes-1; j >= 0; --j)
					out_buffer[out_pos++] = act_kmer->get_byte(j);
	
				uint32 tmp;
				float f_count = (float) count;
				memcpy(&tmp, &f_count, 4);
				for(int32 j = 0; j < 4; ++j)
					out_buffer[out_pos++] = (tmp >> (j * 8)) & 0xFF;
				
				lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;
			}
			++ptr.n_unique;
		}
	
		//// Push the sorted and compacted kmer bin to a priority queue in a form ready to be stored to HDD
		ptr.kq->push(ptr.bin_id, out_buffer, out_pos, raw_lut, lut_size, ptr.n_unique, ptr.n_cutoff_min, ptr.n_cutoff_max, ptr.n_total);
	
		if(ptr.buffer_input)
		{
			ptr.memory_bins->free(ptr.bin_id, CMemoryBins::mba_input_array);
			ptr.memory_bins->free(ptr.bin_id, CMemoryBins::mba_tmp_array);
		}
		ptr.buffer = NULL;
}


//************************************************************************************************************
// CWKmerBinSorter - wrapper for multithreading purposes
//************************************************************************************************************
template <typename KMER_T, unsigned SIZE> class CWKmerBinSorter {
	CKmerBinSorter<KMER_T, SIZE> *kbs;

public:
	CWKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, int thread_no);
	~CWKmerBinSorter();
	void GetDebugStats(uint64& _sum_n_recs, uint64& _sum_n_plus_x_recs)
	{
		kbs->GetDebugStats(_sum_n_recs, _sum_n_plus_x_recs);
	}
	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template <typename KMER_T, unsigned SIZE> CWKmerBinSorter<KMER_T, SIZE>::CWKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, int thread_no)
{
	kbs = new CKmerBinSorter<KMER_T, SIZE>(Params, Queues, thread_no);
}

//----------------------------------------------------------------------------------
// Destructor
template <typename KMER_T, unsigned SIZE> CWKmerBinSorter<KMER_T, SIZE>::~CWKmerBinSorter()
{
	delete kbs;
}

//----------------------------------------------------------------------------------
// Execution
template <typename KMER_T, unsigned SIZE> void CWKmerBinSorter<KMER_T, SIZE>::operator()()
{
	kbs->ProcessBins();
}


#endif
// ***** EOF