/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.0
  Date   : 2014-07-04
*/

#ifndef _BKB_SORTER_H
#define _BKB_SORTER_H

#include "radix.h"
#include "kxmer_set.h"
#include "params.h"

//************************************************************************************************************
template<typename KMER_T, unsigned SIZE> class CBigKmerBinSorter_Impl;

//************************************************************************************************************
// CBigKmerBinSorter - sorter for part of bin, only in strict memory mode
//************************************************************************************************************
template<typename KMER_T, unsigned SIZE>
class CBigKmerBinSorter
{			
	CBigBinKXmersQueue* bbkq;
	CBigBinDesc* bbd;	
	CBigBinSortedPartQueue* bbspq;
	CMemoryPool *pmm_radix_buf, *sm_pmm_expand, *sm_pmm_sorter_suffixes, *sm_pmm_sorter_lut, *sm_pmm_sort;

	uchar* _raw_kxmers;

	int64 sm_mem_part_suffixes;

	CKXmerSet<KMER_T, SIZE> kxmer_set;


	KMER_T* kxmers;
	KMER_T* sort_tmp;
	KMER_T* sorted_kxmers;
	uint32 *kxmers_counters;
	uint64 kxmers_size;
	uint64 kxmers_pos;

	int32 cutoff_min, cutoff_max;
	int32 counter_max;
	uint32 *kxmer_counters;

	int32 lut_prefix_len;
	int n_omp_threads;
	int32 bin_id;
	uint32 sub_bin_id;

	uint32 max_x;
	uint32 kmer_len;
	bool use_quake;
	uint64 sum_n_rec, sum_n_plus_x_rec;

	friend class CBigKmerBinSorter_Impl<KMER_T, SIZE>;

	void Sort();	

public:
	CBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues);
	~CBigKmerBinSorter();
	void Process();
	
};

//************************************************************************************************************
// CBigKmerBinSorter_Impl - implementation of k-mer type- and size-dependent functions
//************************************************************************************************************
template<typename KMER_T, unsigned SIZE>
class CBigKmerBinSorter_Impl
{
public:
	static void PostProcessSort(CBigKmerBinSorter<KMER_T, SIZE>& ptr);
};
template<unsigned SIZE>
class CBigKmerBinSorter_Impl<CKmer<SIZE>, SIZE>
{
	static void PostProcessKmers(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr);
	static void PostProcessKxmers(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr);
	static void PreCompactKxmers(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64& compacted_count, uint32* counters);
	static uint64 FindFirstSymbOccur(CBigKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb);
	static void InitKXMerSet(CBigKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth);
public:
	static void PostProcessSort(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr);
};

template<unsigned SIZE>
class CBigKmerBinSorter_Impl<CKmerQuake<SIZE>, SIZE>
{
public:
	static void PostProcessSort(CBigKmerBinSorter<CKmerQuake<SIZE>, SIZE>& ptr);
};



//************************************************************************************************************
// CBigKmerBinSorter
//************************************************************************************************************

//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE> void CBigKmerBinSorter<KMER_T, SIZE>::Process()
{
	int32 curr_bin_id = -1;
	bin_id = -1;
	uchar* data = NULL;
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
				CBigKmerBinSorter_Impl<KMER_T, SIZE>::PostProcessSort(*this);				
				kxmers_pos = 0;
			}
			bin_id = curr_bin_id;
			sub_bin_id = 0;			
		}

		if (kxmers_pos + size < kxmers_size)
		{			
			A_memcpy(kxmers + kxmers_pos, data, size * sizeof(KMER_T));
			sm_pmm_expand->free(data);
			kxmers_pos += size;
		}
		else
		{
			Sort();
			CBigKmerBinSorter_Impl<KMER_T, SIZE>::PostProcessSort(*this);
			++sub_bin_id;
			A_memcpy(kxmers, data, size * sizeof(KMER_T));
			sm_pmm_expand->free(data);
			kxmers_pos = size;
		}
	}
	if (kxmers_pos)
	{
		Sort();
		CBigKmerBinSorter_Impl<KMER_T, SIZE>::PostProcessSort(*this);
	}	
	bbspq->mark_completed();
}

//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE>
CBigKmerBinSorter<KMER_T, SIZE>::CBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues) : kxmer_set(Params.kmer_len)
{	
	sorted_kxmers = NULL;
	kxmer_counters = NULL;
	bbkq = Queues.bbkq;	
	bbspq = Queues.bbspq;
	pmm_radix_buf = Queues.pmm_radix_buf;
	sm_pmm_expand = Queues.sm_pmm_expand;
	sm_pmm_sorter_suffixes = Queues.sm_pmm_sorter_suffixes;
	sm_pmm_sorter_lut = Queues.sm_pmm_sorter_lut;
	sm_pmm_sort = Queues.sm_pmm_sort;

	kxmers_size = Params.sm_mem_part_sort / 2 / sizeof(KMER_T);

	sm_mem_part_suffixes = Params.sm_mem_part_suffixes;
		
	sm_pmm_sort->reserve(_raw_kxmers);
	kxmers = (KMER_T*)_raw_kxmers;
	
	sort_tmp = kxmers + kxmers_size;
	max_x = Params.max_x;
	bbd = Queues.bbd;
	use_quake = Params.use_quake;
	kmer_len = Params.kmer_len;	

	lut_prefix_len = Params.lut_prefix_len;

	n_omp_threads = Params.sm_n_omp_threads;
	
	sum_n_rec = sum_n_plus_x_rec = 0;
	cutoff_max = Params.cutoff_max;
	cutoff_min = Params.cutoff_min;
	counter_max = Params.counter_max;
}

//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE> CBigKmerBinSorter<KMER_T, SIZE>::~CBigKmerBinSorter()
{	
	sm_pmm_sort->free(_raw_kxmers);
}


//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE>
void CBigKmerBinSorter<KMER_T, SIZE>::Sort()
{
	uint32 rec_len;
	uint64 sort_rec = kxmers_pos;	
	if (max_x && !use_quake)
	{
		rec_len = (kmer_len + max_x + 1 + 3) / 4;
	}
	else
	{
		rec_len = (kmer_len + 3) / 4;
	}
	sum_n_plus_x_rec += kxmers_pos;

	if (sizeof(KMER_T) == 8)
	{
		uint64 *_buffer_input = (uint64*)kxmers;
		uint64 *_buffer_tmp = (uint64*)sort_tmp;
	
		RadixSort_buffer(pmm_radix_buf, _buffer_input, _buffer_tmp, sort_rec, rec_len, n_omp_threads);

		if (rec_len % 2)
		{
			kxmers_counters = (uint32*)kxmers;
			sorted_kxmers = (KMER_T*)sort_tmp;
		}
		else
		{
			kxmers_counters = (uint32*)sort_tmp;
			sorted_kxmers = (KMER_T*)kxmers;
		}
	}
	else
	{
		uint32 *_buffer_input = (uint32*)kxmers;
		uint32 *_buffer_tmp = (uint32*)sort_tmp;

		RadixSort_uint8(_buffer_input, _buffer_tmp, sort_rec, sizeof(KMER_T), offsetof(KMER_T, data), SIZE*sizeof(typename KMER_T::data_t), rec_len, n_omp_threads);

		if (rec_len % 2)
		{
			kxmers_counters = (uint32*)_buffer_input;
			sorted_kxmers = (KMER_T*)_buffer_tmp;
		}
		else
		{
			kxmers_counters = (uint32*)_buffer_tmp;
			sorted_kxmers = (KMER_T*)_buffer_input;
		}
	}
}


//************************************************************************************************************
// CBigKmerBinSorter_Impl
//************************************************************************************************************

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::PostProcessKmers(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr)
{
	uint32 best_lut_prefix_len = 0;
	uint32 local_lut_prefix_len;
	uint64 best_mem_amount = 1ull << 62;

	uint32 counter_size = sizeof(uint32);

	for (local_lut_prefix_len = 2; local_lut_prefix_len < 13; ++local_lut_prefix_len)
	{
		uint32 suffix_len = ptr.kmer_len - local_lut_prefix_len;
		if (suffix_len % 4)
			continue;

		uint64 suf_mem = (suffix_len / 4 + counter_size) * ptr.kxmers_pos;
		uint64 lut_mem = (1ull << (2 * local_lut_prefix_len)) * sizeof(uint64);
		if (suf_mem + lut_mem < best_mem_amount)
		{
			best_mem_amount = suf_mem + lut_mem;
			best_lut_prefix_len = local_lut_prefix_len;
		}
	}
	local_lut_prefix_len = best_lut_prefix_len;

	uint32 kmer_symbols = ptr.kmer_len - local_lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;

	uint32 suffix_rec_bytes = (ptr.kmer_len - local_lut_prefix_len) / 4 + counter_size;
	uint64 lut_recs = 1ull << 2 * local_lut_prefix_len;


	uchar* suff_buff;
	ptr.sm_pmm_sorter_suffixes->reserve(suff_buff);
	uchar* _raw_lut;
	ptr.sm_pmm_sorter_lut->reserve(_raw_lut);
	uint64* lut = (uint64*)_raw_lut;
	fill_n(lut, lut_recs, 0);

	uint64 suff_buff_size = ptr.sm_mem_part_suffixes / suffix_rec_bytes * suffix_rec_bytes;

	uint64 suff_buff_pos = 0;
	uint32 n_recs = 0;
	CKmer<SIZE> *act_kmer;
	uint32 count;
	uint64 i;
	act_kmer = &ptr.kxmers[0];
	count = 1;	
	for (i = 1; i < ptr.kxmers_pos; ++i)
	{
		if (*act_kmer == ptr.kxmers[i])
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
				ptr.bbspq->push(ptr.bin_id, ptr.sub_bin_id, suff_buff, suff_buff_pos, NULL, 0, false);
				ptr.sm_pmm_sorter_suffixes->reserve(suff_buff);
				suff_buff_pos = 0;
			}

			count = 1;
			act_kmer = &ptr.kxmers[i];
		}
	}

	lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;

	for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
		suff_buff[suff_buff_pos++] = act_kmer->get_byte(j);
	for (int32 j = 0; j < (int32)counter_size; ++j)
		suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;

	++n_recs;

	ptr.bbspq->push(ptr.bin_id, ptr.sub_bin_id, suff_buff, suff_buff_pos, NULL, 0, false);
	ptr.bbspq->push(ptr.bin_id, ptr.sub_bin_id, NULL, 0, lut, lut_recs, true);
	ptr.bbd->push(ptr.bin_id, ptr.sub_bin_id, local_lut_prefix_len, n_recs, NULL, "", 0);
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::PreCompactKxmers(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr, uint64& compacted_count, uint32* counters)
{
	compacted_count = 0;

	CKmer<SIZE> *act_kmer;
	act_kmer = &ptr.sorted_kxmers[0];
	counters[compacted_count] = 1;

	for (uint32 i = 1; i < ptr.kxmers_pos; ++i)
	{
		if (*act_kmer == ptr.sorted_kxmers[i])
			++counters[compacted_count];
		else
		{
			ptr.sorted_kxmers[compacted_count++] = *act_kmer;
			counters[compacted_count] = 1;
			act_kmer = &ptr.sorted_kxmers[i];
		}
	}
	ptr.sorted_kxmers[compacted_count++] = *act_kmer;
}

//----------------------------------------------------------------------------------
//Binary search position of first occurence of symbol 'symb' in [start_pos,end_pos). Offset defines which symbol in k+x-mer is taken.
template <unsigned SIZE> uint64 CBigKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::FindFirstSymbOccur(CBigKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb)
{
	uint32 kxmer_offset = (ptr.kmer_len + ptr.max_x - offset) * 2;
	uint64 middle_pos;
	uchar middle_symb;
	while (start_pos < end_pos)
	{
		middle_pos = (start_pos + end_pos) / 2;
		middle_symb = ptr.sorted_kxmers[middle_pos].get_2bits(kxmer_offset);
		if (middle_symb < symb)
			start_pos = middle_pos + 1;
		else
			end_pos = middle_pos;
	}
	return end_pos;
}


//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::InitKXMerSet(CBigKmerBinSorter<CKmer<SIZE>, SIZE> &ptr, uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth)
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
template<unsigned SIZE> void CBigKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::PostProcessKxmers(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr)
{
	ptr.kxmer_set.clear();
	ptr.kxmer_set.set_buffer(ptr.sorted_kxmers);

	uint32 best_lut_prefix_len = 0;
	uint32 local_lut_prefix_len;
	uint64 best_mem_amount = 1ull << 62;

	uint32 counter_size = sizeof(uint32);

	for (local_lut_prefix_len = 2; local_lut_prefix_len < 13; ++local_lut_prefix_len) 
	{
		uint32 suffix_len = ptr.kmer_len - local_lut_prefix_len;
		if(suffix_len % 4)
			continue;

		uint64 suf_mem = (suffix_len / 4 + counter_size) * ptr.kxmers_pos;
		uint64 lut_mem = (1ull << (2 * local_lut_prefix_len)) * sizeof(uint64);
		if (suf_mem + lut_mem < best_mem_amount)
		{
			best_mem_amount = suf_mem + lut_mem;
			best_lut_prefix_len = local_lut_prefix_len;
		}
	}
	local_lut_prefix_len = best_lut_prefix_len;


	uint32 kmer_symbols = ptr.kmer_len - local_lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;

	uint32 suffix_rec_bytes = (ptr.kmer_len - local_lut_prefix_len) / 4 + counter_size;
	uint64 lut_recs = 1ull << 2 * local_lut_prefix_len;
	

	uchar* suff_buff;
	ptr.sm_pmm_sorter_suffixes->reserve(suff_buff);
	uchar* _raw_lut;
	ptr.sm_pmm_sorter_lut->reserve(_raw_lut);
	uint64* lut = (uint64*)_raw_lut;
	fill_n(lut, lut_recs, 0);

	uint64 suff_buff_size = ptr.sm_mem_part_suffixes / suffix_rec_bytes * suffix_rec_bytes;
	
	uint64 suff_buff_pos = 0;
	uint32 n_recs = 0;

	uint64 compacted_count;
	PreCompactKxmers(ptr, compacted_count, ptr.kxmers_counters);
	

	uint64 pos[5];
	pos[0] = 0;
	pos[4] = compacted_count;
	for(uint32 i = 1 ; i < 4 ; ++i)
		pos[i] = FindFirstSymbOccur(ptr, pos[i - 1], compacted_count, 0, i);
	for (uint32 i = 1; i < 5; ++i)
		InitKXMerSet(ptr, pos[i - 1], pos[i], ptr.max_x + 2 - i, i);


	uint64 counter_pos = 0;

	CKmer<SIZE> kmer, next_kmer;
	kmer.clear();
	next_kmer.clear();
	CKmer<SIZE> kmer_mask;
	uint32 count;
	kmer_mask.set_n_1(ptr.kmer_len * 2);
	ptr.kxmer_set.get_min(counter_pos, kmer);
	count = ptr.kxmers_counters[counter_pos];

	while (ptr.kxmer_set.get_min(counter_pos, next_kmer))
	{
		if (kmer == next_kmer)
			count += ptr.kxmers_counters[counter_pos];
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
				ptr.bbspq->push(ptr.bin_id, ptr.sub_bin_id, suff_buff, suff_buff_pos, NULL, 0, false);				
				ptr.sm_pmm_sorter_suffixes->reserve(suff_buff);
				suff_buff_pos = 0;
			}
			
			count = ptr.kxmers_counters[counter_pos];
			kmer = next_kmer;
		}
	}
			
	lut[kmer.remove_suffix(2 * kmer_symbols)]++;

	for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
		suff_buff[suff_buff_pos++] = kmer.get_byte(j);
	for (int32 j = 0; j < (int32)counter_size; ++j)
		suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;

	++n_recs;

	ptr.bbspq->push(ptr.bin_id, ptr.sub_bin_id, suff_buff, suff_buff_pos, NULL, 0, false);
	ptr.bbspq->push(ptr.bin_id, ptr.sub_bin_id, NULL, 0, lut, lut_recs, true);
	ptr.bbd->push(ptr.bin_id, ptr.sub_bin_id, local_lut_prefix_len, n_recs, NULL, "", 0);	
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter_Impl<CKmer<SIZE>, SIZE>::PostProcessSort(CBigKmerBinSorter<CKmer<SIZE>, SIZE>& ptr)
{
	if (ptr.max_x)
		PostProcessKxmers(ptr);
	else
		PostProcessKmers(ptr);
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinSorter_Impl<CKmerQuake<SIZE>, SIZE>::PostProcessSort(CBigKmerBinSorter<CKmerQuake<SIZE>, SIZE>& ptr)
{
	//"Not supported in current release"
}




//************************************************************************************************************
// CWBigKmerBinSorter - wrapper for multithreading purposes
//************************************************************************************************************
template<typename KMER_T, unsigned SIZE>
class CWBigKmerBinSorter
{
	CBigKmerBinSorter<KMER_T, SIZE>* bkb_sorter;
public:
	CWBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues);
	~CWBigKmerBinSorter();
	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template<typename KMER_T, unsigned SIZE>
CWBigKmerBinSorter<KMER_T, SIZE>::CWBigKmerBinSorter(CKMCParams& Params, CKMCQueues& Queues)
{
	bkb_sorter = new CBigKmerBinSorter<KMER_T, SIZE>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
template<typename KMER_T, unsigned SIZE>
CWBigKmerBinSorter<KMER_T, SIZE>::~CWBigKmerBinSorter()
{
	delete bkb_sorter;
}

//----------------------------------------------------------------------------------
// Execution
template<typename KMER_T, unsigned SIZE>
void CWBigKmerBinSorter<KMER_T, SIZE>::operator()()
{
	bkb_sorter->Process();
}
#endif  

// ***** EOF 