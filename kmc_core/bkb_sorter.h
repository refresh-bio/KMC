/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _BKB_SORTER_H
#define _BKB_SORTER_H

#include "radix.h"
#include "raduls.h"
#include "kxmer_set.h"
#include "params.h"

//************************************************************************************************************
// CBigKmerBinSorter - sorter for part of bin, only in strict memory mode
//************************************************************************************************************
template<unsigned SIZE>
class CBigKmerBinSorter
{			
	CBigBinKXmersQueue* bbkq;
	CBigBinDesc* bbd;	
	CBigBinSortedPartQueue* bbspq;
	CMemoryPool *pmm_radix_buf, *sm_pmm_expand, *sm_pmm_sorter_suffixes, *sm_pmm_sorter_lut, *sm_pmm_sort;

	uchar* _raw_kxmers;

	int64 sm_mem_part_suffixes;

	CKXmerSet<SIZE> kxmer_set;


	CKmer<SIZE>* kxmers;
	CKmer<SIZE>* sort_tmp;
	CKmer<SIZE>* sorted_kxmers;
	uint32 *kxmers_counters;
	uint64 kxmers_size;
	uint64 kxmers_pos;
	
	uint32 *kxmer_counters;

	int32 lut_prefix_len;
	int n_sorting_threads;
	int32 bin_id;
	uint32 sub_bin_id;

	uint32 max_x;
	uint32 kmer_len;	
	uint64 sum_n_rec, sum_n_plus_x_rec;

	SortFunction<CKmer<SIZE>> sort_func;

	void PostProcessKmers();
	void PostProcessKxmers();
	void PreCompactKxmers(uint64& compacted_count, uint32* counters);
	uint64 FindFirstSymbOccur(uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb);
	void InitKXMerSet(uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth);
	void PostProcessSort();

	void Sort();	

public:
	CBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues, SortFunction<CKmer<SIZE>> sort_func);
	~CBigKmerBinSorter();
	void Process();
	
};


//************************************************************************************************************
// CBigKmerBinSorter
//************************************************************************************************************

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter<SIZE>::Process()
{
	int32 curr_bin_id = -1;
	bin_id = -1;
	uchar* data = nullptr;
	uint64 size = 0;	
	kxmers_pos = 0;
	sub_bin_id = 0;
	
	while (bbkq->pop(curr_bin_id, data, size))
	{		
		if (bin_id == -1)
			bin_id = curr_bin_id;

		if (curr_bin_id != bin_id) //new bin
		{
			if (kxmers_pos)
			{			
				Sort();			
				PostProcessSort();				
				kxmers_pos = 0;
			}
			bin_id = curr_bin_id;
			sub_bin_id = 0;			
		}

		if (kxmers_pos + size < kxmers_size)
		{			
			memcpy(kxmers + kxmers_pos, data, size * sizeof(CKmer<SIZE>));
			sm_pmm_expand->free(data);
			kxmers_pos += size;
		}
		else
		{
			Sort();
			PostProcessSort();
			++sub_bin_id;
			memcpy(kxmers, data, size * sizeof(CKmer<SIZE>));
			sm_pmm_expand->free(data);
			kxmers_pos = size;
		}
	}
	if (kxmers_pos)
	{
		Sort();
		PostProcessSort();
	}	
	bbspq->mark_completed();
}

//----------------------------------------------------------------------------------
template<unsigned SIZE>
CBigKmerBinSorter<SIZE>::CBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues, SortFunction<CKmer<SIZE>> sort_func) : 
	kxmer_set(Params.kmer_len), 
	sort_func(sort_func)
{	
	sorted_kxmers = nullptr;
	kxmer_counters = nullptr;
	bbkq = Queues.bbkq;	
	bbspq = Queues.bbspq;
	pmm_radix_buf = Queues.pmm_radix_buf;
	sm_pmm_expand = Queues.sm_pmm_expand;
	sm_pmm_sorter_suffixes = Queues.sm_pmm_sorter_suffixes;
	sm_pmm_sorter_lut = Queues.sm_pmm_sorter_lut;
	sm_pmm_sort = Queues.sm_pmm_sort;

	kxmers_size = Params.sm_mem_part_sort / 2 / sizeof(CKmer<SIZE>);

	sm_mem_part_suffixes = Params.sm_mem_part_suffixes;
		
	sm_pmm_sort->reserve(_raw_kxmers);
	kxmers = (CKmer<SIZE>*)_raw_kxmers;
	
	sort_tmp = kxmers + kxmers_size;
	max_x = Params.max_x;
	bbd = Queues.bbd;	
	kmer_len = Params.kmer_len;	

	lut_prefix_len = Params.lut_prefix_len;

	n_sorting_threads = Params.sm_n_sorting_threads;
	
	sum_n_rec = sum_n_plus_x_rec = 0;	
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> CBigKmerBinSorter<SIZE>::~CBigKmerBinSorter()
{	
	sm_pmm_sort->free(_raw_kxmers);
}


//----------------------------------------------------------------------------------
template<unsigned SIZE>
void CBigKmerBinSorter<SIZE>::Sort()
{
	uint32 rec_len;
	uint64 sort_rec = kxmers_pos;	
	if (max_x)
	{
		rec_len = (kmer_len + max_x + 1 + 3) / 4;
	}
	else
	{
		rec_len = (kmer_len + 3) / 4;
	}
	sum_n_plus_x_rec += kxmers_pos;

	sort_func(kxmers, sort_tmp, sort_rec, rec_len - 1, n_sorting_threads, pmm_radix_buf);
	
	if (rec_len % 2)
	{
		kxmers_counters = (uint32*)kxmers;
		sorted_kxmers = sort_tmp;
	}
	else
	{
		kxmers_counters = (uint32*)sort_tmp;
		sorted_kxmers = kxmers;
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter<SIZE>::PostProcessKmers()
{
	uint32 best_lut_prefix_len = 0;
	uint32 local_lut_prefix_len;
	uint64 best_mem_amount = 1ull << 62;

	uint32 counter_size = sizeof(uint32);

	for (local_lut_prefix_len = 2; local_lut_prefix_len < 13; ++local_lut_prefix_len)
	{
		uint32 suffix_len = kmer_len - local_lut_prefix_len;
		if (suffix_len % 4)
			continue;

		uint64 suf_mem = (suffix_len / 4 + counter_size) * kxmers_pos;
		uint64 lut_mem = (1ull << (2 * local_lut_prefix_len)) * sizeof(uint64);
		if (suf_mem + lut_mem < best_mem_amount)
		{
			best_mem_amount = suf_mem + lut_mem;
			best_lut_prefix_len = local_lut_prefix_len;
		}
	}
	local_lut_prefix_len = best_lut_prefix_len;

	uint32 kmer_symbols = kmer_len - local_lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;

	uint32 suffix_rec_bytes = (kmer_len - local_lut_prefix_len) / 4 + counter_size;
	uint64 lut_recs = 1ull << 2 * local_lut_prefix_len;


	uchar* suff_buff;
	sm_pmm_sorter_suffixes->reserve(suff_buff);
	uchar* _raw_lut;
	sm_pmm_sorter_lut->reserve(_raw_lut);
	uint64* lut = (uint64*)_raw_lut;
	fill_n(lut, lut_recs, 0);

	uint64 suff_buff_size = sm_mem_part_suffixes / suffix_rec_bytes * suffix_rec_bytes;

	uint64 suff_buff_pos = 0;
	uint64 n_recs = 0;
	CKmer<SIZE> *act_kmer;
	uint32 count;
	uint64 i;
	act_kmer = &kxmers[0];
	count = 1;	
	for (i = 1; i < kxmers_pos; ++i)
	{
		if (*act_kmer == kxmers[i])
			count++;
		else
		{
			lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;
			for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
				suff_buff[suff_buff_pos++] = act_kmer->get_byte(j);
			for (int32 j = 0; j < (int32)counter_size; ++j)
				suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;
			++n_recs;

			if (suff_buff_pos >= suff_buff_size)
			{
				bbspq->push(bin_id, sub_bin_id, suff_buff, suff_buff_pos, nullptr, 0, false);
				sm_pmm_sorter_suffixes->reserve(suff_buff);
				suff_buff_pos = 0;
			}

			count = 1;
			act_kmer = &kxmers[i];
		}
	}

	lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;

	for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
		suff_buff[suff_buff_pos++] = act_kmer->get_byte(j);
	for (int32 j = 0; j < (int32)counter_size; ++j)
		suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;

	++n_recs;

	bbspq->push(bin_id, sub_bin_id, suff_buff, suff_buff_pos, nullptr, 0, false);
	bbspq->push(bin_id, sub_bin_id, nullptr, 0, lut, lut_recs, true);
	bbd->push(bin_id, sub_bin_id, local_lut_prefix_len, n_recs, nullptr, "", 0);
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter<SIZE>::PreCompactKxmers(uint64& compacted_count, uint32* counters)
{
	compacted_count = 0;

	CKmer<SIZE> *act_kmer;
	act_kmer = &sorted_kxmers[0];
	counters[compacted_count] = 1;

	for (uint32 i = 1; i < kxmers_pos; ++i)
	{
		if (*act_kmer == sorted_kxmers[i])
			++counters[compacted_count];
		else
		{
			sorted_kxmers[compacted_count++] = *act_kmer;
			counters[compacted_count] = 1;
			act_kmer = &sorted_kxmers[i];
		}
	}
	sorted_kxmers[compacted_count++] = *act_kmer;
}

//----------------------------------------------------------------------------------
//Binary search position of first occurrence of symbol 'symb' in [start_pos,end_pos). Offset defines which symbol in k+x-mer is taken.
template <unsigned SIZE> uint64 CBigKmerBinSorter<SIZE>::FindFirstSymbOccur(uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb)
{
	uint32 kxmer_offset = (kmer_len + max_x - offset) * 2;
	uint64 middle_pos;
	uchar middle_symb;
	while (start_pos < end_pos)
	{
		middle_pos = (start_pos + end_pos) / 2;
		middle_symb = sorted_kxmers[middle_pos].get_2bits(kxmer_offset);
		if (middle_symb < symb)
			start_pos = middle_pos + 1;
		else
			end_pos = middle_pos;
	}
	return end_pos;
}


//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter<SIZE>::InitKXMerSet(uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth)
{
	if (start_pos == end_pos)
		return;
	uint32 shr = max_x + 1 - offset;
	kxmer_set.init_add(start_pos, end_pos, shr);

	--depth;
	if (depth > 0)
	{
		uint64 pos[5];
		pos[0] = start_pos;
		pos[4] = end_pos;
		for (uint32 i = 1; i < 4; ++i)
			pos[i] = FindFirstSymbOccur(pos[i - 1], end_pos, offset, i);
		for (uint32 i = 1; i < 5; ++i)
			InitKXMerSet(pos[i - 1], pos[i], offset + 1, depth);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter<SIZE>::PostProcessKxmers()
{
	kxmer_set.clear();
	kxmer_set.set_buffer(sorted_kxmers);

	uint32 best_lut_prefix_len = 0;
	uint32 local_lut_prefix_len;
	uint64 best_mem_amount = 1ull << 62;

	uint32 counter_size = sizeof(uint32);

	for (local_lut_prefix_len = 2; local_lut_prefix_len < 13; ++local_lut_prefix_len) 
	{
		uint32 suffix_len = kmer_len - local_lut_prefix_len;
		if(suffix_len % 4)
			continue;

		uint64 suf_mem = (suffix_len / 4 + counter_size) * kxmers_pos;
		uint64 lut_mem = (1ull << (2 * local_lut_prefix_len)) * sizeof(uint64);
		if (suf_mem + lut_mem < best_mem_amount)
		{
			best_mem_amount = suf_mem + lut_mem;
			best_lut_prefix_len = local_lut_prefix_len;
		}
	}
	local_lut_prefix_len = best_lut_prefix_len;


	uint32 kmer_symbols = kmer_len - local_lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;

	uint32 suffix_rec_bytes = (kmer_len - local_lut_prefix_len) / 4 + counter_size;
	uint64 lut_recs = 1ull << 2 * local_lut_prefix_len;
	

	uchar* suff_buff;
	sm_pmm_sorter_suffixes->reserve(suff_buff);
	uchar* _raw_lut;
	sm_pmm_sorter_lut->reserve(_raw_lut);
	uint64* lut = (uint64*)_raw_lut;
	fill_n(lut, lut_recs, 0);

	uint64 suff_buff_size = sm_mem_part_suffixes / suffix_rec_bytes * suffix_rec_bytes;
	
	uint64 suff_buff_pos = 0;
	uint64 n_recs = 0;

	uint64 compacted_count;
	PreCompactKxmers(compacted_count, kxmers_counters);
	

	uint64 pos[5];
	pos[0] = 0;
	pos[4] = compacted_count;
	for(uint32 i = 1 ; i < 4 ; ++i)
		pos[i] = FindFirstSymbOccur(pos[i - 1], compacted_count, 0, i);
	for (uint32 i = 1; i < 5; ++i)
		InitKXMerSet(pos[i - 1], pos[i], max_x + 2 - i, i);


	uint64 counter_pos = 0;

	CKmer<SIZE> kmer, next_kmer;
	kmer.clear();
	next_kmer.clear();
	CKmer<SIZE> kmer_mask;
	uint32 count;
	kmer_mask.set_n_1(kmer_len * 2);
	kxmer_set.get_min(counter_pos, kmer);
	count = kxmers_counters[counter_pos];

	while (kxmer_set.get_min(counter_pos, next_kmer))
	{
		if (kmer == next_kmer)
			count += kxmers_counters[counter_pos];
		else
		{
			lut[kmer.remove_suffix(2 * kmer_symbols)]++;			
			for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
				suff_buff[suff_buff_pos++] = kmer.get_byte(j);
			for (int32 j = 0; j < (int32)counter_size; ++j)
				suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;
			++n_recs;

			if (suff_buff_pos >= suff_buff_size)
			{				
				bbspq->push(bin_id, sub_bin_id, suff_buff, suff_buff_pos, nullptr, 0, false);
				sm_pmm_sorter_suffixes->reserve(suff_buff);
				suff_buff_pos = 0;
			}
			
			count = kxmers_counters[counter_pos];
			kmer = next_kmer;
		}
	}
			
	lut[kmer.remove_suffix(2 * kmer_symbols)]++;

	for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
		suff_buff[suff_buff_pos++] = kmer.get_byte(j);
	for (int32 j = 0; j < (int32)counter_size; ++j)
		suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;

	++n_recs;

	bbspq->push(bin_id, sub_bin_id, suff_buff, suff_buff_pos, nullptr, 0, false);
	bbspq->push(bin_id, sub_bin_id, nullptr, 0, lut, lut_recs, true);
	bbd->push(bin_id, sub_bin_id, local_lut_prefix_len, n_recs, nullptr, "", 0);
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter<SIZE>::PostProcessSort()
{
	if (max_x)
		PostProcessKxmers();
	else
		PostProcessKmers();
}


//************************************************************************************************************
// CWBigKmerBinSorter - wrapper for multithreading purposes
//************************************************************************************************************
template<unsigned SIZE>
class CWBigKmerBinSorter
{
	CBigKmerBinSorter<SIZE>* bkb_sorter;
public:
	CWBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues, SortFunction<CKmer<SIZE>> sort_func);
	~CWBigKmerBinSorter();
	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template<unsigned SIZE>
CWBigKmerBinSorter<SIZE>::CWBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues, SortFunction<CKmer<SIZE>> sort_func)
{
	bkb_sorter = new CBigKmerBinSorter<SIZE>(Params, Queues, sort_func);
}

//----------------------------------------------------------------------------------
// Destructor
template<unsigned SIZE>
CWBigKmerBinSorter<SIZE>::~CWBigKmerBinSorter()
{
	delete bkb_sorter;
}

//----------------------------------------------------------------------------------
// Execution
template<unsigned SIZE>
void CWBigKmerBinSorter<SIZE>::operator()()
{
	bkb_sorter->Process();
}
#endif  

// ***** EOF 