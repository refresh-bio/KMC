/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _KB_SORTER_H
#define _KB_SORTER_H

#define DEBUGG_INFO

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include "raduls.h"
#include "radix.h"
#include "s_mapper.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <vector>
#include <stdio.h>
#include <functional>
#include <cstddef>
#include <set>


#include "kxmer_set.h"
#include "rev_byte.h"


template<unsigned SIZE> class CExpandThread;

//************************************************************************************************************
// CKmerBinSorter - sorting of k-mers in a bin
//************************************************************************************************************
template <unsigned SIZE> class CKmerBinSorter {	
private:

	mutable mutex expander_mtx;
	uint64 input_pos;

	CMemoryMonitor *mm;
	CBinDesc *bd;
	CExpanderPackDesc *epd;
	CBinQueue *bq;
	CKmerQueue *kq;
	CMemoryPool *pmm_radix_buf;
	CMemoryBins *memory_bins;	
	CSortersManager* sorters_manager;

	CKXmerSet<SIZE> kxmer_set;
	SortFunction<CKmer<SIZE>> sort_func;

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

	int n_sorting_threads;

	bool both_strands;	
	bool without_output;

	CSignatureMapper* s_mapper;

	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total;
	uint32 cutoff_min, cutoff_max;
	int32 lut_prefix_len;
	uint32 counter_max;

	CKmer<SIZE> *buffer_input, *buffer_tmp, *buffer;
	uint32 *kxmer_counters;

	void Sort();

	friend class CExpandThread<SIZE>;

	uint64 FindFirstSymbOccur(uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb);
	void InitKXMerSet(uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth);
	void InitKXMerSetMultithreaded(CKXmerSetMultiThreaded<SIZE>& kxmer_set_multithreaded, uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth);
	void CompactKxmers();
	void PreCompactKxmers(uint64& compacted_count);
	void CompactKmers();
	void ExpandKxmersAll(uint64 tmp_size);
	void ExpandKxmersBoth(uint64 tmp_size);
	void ExpandKmersAll(uint64 tmp_size);
	void ExpandKmersBoth(uint64 tmp_size);
	void GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p);
	void FromChildThread(CKmer<SIZE>* thread_buffer, uint64 size);
	uint64 ExpandKxmerBothParallel(uint64 start_pos, uint64 end_pos, uint64 output_start, uint64 output_end);

	void Compact();
	void Expand(uint64 tmp_size);


public:
	static uint32 PROB_BUF_SIZE;
	CKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, SortFunction<CKmer<SIZE>> sort_func);
	~CKmerBinSorter();

	void GetDebugStats(uint64& _sum_n_recs, uint64& _sum_n_plus_x_recs)
	{
		_sum_n_recs = sum_n_rec;
		_sum_n_plus_x_recs = sum_n_plus_x_rec;		
	}

	void ProcessBins();
};

template <unsigned SIZE> uint32 CKmerBinSorter<SIZE>::PROB_BUF_SIZE = 1 << 14;

// Queue for expander
class CExpanderPackQueue
{
	uint64 _start = 0;
	uint64 _output_start = 0;
	list<pair<uint64, uint64>> l;
	mutex mtx;
public:
	CExpanderPackQueue(list<pair<uint64, uint64>>& expander_pack)
	{
		l = std::move(expander_pack);
	}
	bool Pop(uint64& start, uint64& end, uint64& output_start, uint64& output_end)
	{
		lock_guard<mutex> lck(mtx);
		if (l.empty())
			return false;
		start = _start;
		_start += l.front().first;
		end = _start;
		output_start = _output_start;
		_output_start += l.front().second;
		output_end = _output_start;

		l.pop_front();
		return true;
	}
};


//************************************************************************************************************
// CKmerBinSorter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Assign queues and monitors
template <unsigned SIZE> CKmerBinSorter<SIZE>::CKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, SortFunction<CKmer<SIZE>> sort_func)
:
	kxmer_set(Params.kmer_len),
	sort_func(sort_func)
{
	both_strands = Params.both_strands;
	mm = Queues.mm;
	n_bins = Params.n_bins;
	bd = Queues.bd;
	epd = Queues.epd;
	bq = Queues.bq;
	kq = Queues.kq;
	
	sorters_manager = Queues.sorters_manager;

	s_mapper = Queues.s_mapper;

	pmm_radix_buf = Queues.pmm_radix_buf;	
	
	memory_bins = Queues.memory_bins;

	cutoff_min = Params.cutoff_min;
	cutoff_max = (uint32)Params.cutoff_max;
	counter_max = (uint32)Params.counter_max;
	max_x = Params.max_x;	
	without_output = Params.without_output;
	
	lut_prefix_len = Params.lut_prefix_len;

	n_sorting_threads = 0;
	//n_sorting_threads = Params.n_sorting_threads[thread_no];

	sum_n_rec = sum_n_plus_x_rec = 0;
}

//----------------------------------------------------------------------------------
template <unsigned SIZE> CKmerBinSorter<SIZE>::~CKmerBinSorter()
{

}

//----------------------------------------------------------------------------------
// Process the bins
template <unsigned SIZE> void CKmerBinSorter<SIZE>::ProcessBins()
{
	uint64 tmp_size;
	uint64 tmp_n_rec;
	CMemDiskFile *file;
	
	while (sorters_manager->GetNext(bin_id, data, size, n_rec, n_sorting_threads))
	{
		// Get bin data
		bd->read(bin_id, file, desc, tmp_size, tmp_n_rec, n_plus_x_recs, buffer_size, kmer_len);

		// Uncompact the kmers - append truncate prefixes		

		Expand(tmp_size);
		
		memory_bins->free(bin_id, CMemoryBins::mba_input_file);

		// Perform sorting of kmers in a bin
		Sort();
		
		// Compact the same kmers (occurring at neighbour positions now)
		Compact();
				
		sorters_manager->ReturnThreads(n_sorting_threads, bin_id);
	}

	kq->mark_completed();
}

template <unsigned SIZE> inline void CKmerBinSorter<SIZE>::GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p)
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

template <unsigned SIZE> void CKmerBinSorter<SIZE>::ExpandKmersAll(uint64 tmp_size)
{
	uint64 pos = 0;
	input_pos = 0;
	CKmer<SIZE> kmer;
	uint32 kmer_bytes = (kmer_len + 3) / 4;

	CKmer<SIZE> kmer_mask;
	kmer_mask.set_n_1(kmer_len * 2);
	uchar *data_p = data;
	uchar additional_symbols;
	uint32 kmer_shr = SIZE * 32 - kmer_len;
	while (pos < tmp_size)
	{
		kmer.clear();
		additional_symbols = data_p[pos++];		
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1; i < kmer_bytes; ++i, --kmer_pos)
		{
			kmer.set_byte(kmer_pos, data_p[pos + i]);
		}
		pos += kmer_bytes;
		uchar byte_shift = 6 - (kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		buffer_input[input_pos++].set(kmer);
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
			buffer_input[input_pos++].set(kmer);
		}
		if (byte_shift != 6)
			++pos;
	}
}
template <unsigned SIZE> void CKmerBinSorter<SIZE>::ExpandKmersBoth(uint64 tmp_size)
{
	uint64 pos = 0;
	CKmer<SIZE> kmer;
	CKmer<SIZE> rev_kmer;
	CKmer<SIZE> kmer_can;

	uint32 kmer_bytes = (kmer_len + 3) / 4;
	uint32 kmer_len_shift = (kmer_len - 1) * 2;
	CKmer<SIZE> kmer_mask;
	kmer_mask.set_n_1(kmer_len * 2);
	uchar *data_p = data;
	input_pos = 0;
	uint32 kmer_shr = SIZE * 32 - kmer_len;

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
		uchar byte_shift = 6 - (kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		rev_kmer.mask(kmer_mask);

		kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
		buffer_input[input_pos++].set(kmer_can);

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
			buffer_input[input_pos++].set(kmer_can);
		}
		if (byte_shift != 6)
			++pos;
	}
}

template <unsigned SIZE> void CKmerBinSorter<SIZE>::FromChildThread(CKmer<SIZE>* thread_buffer, uint64 size)
{
	lock_guard<mutex> lcx(expander_mtx);
	memcpy(buffer_input + input_pos, thread_buffer, size * sizeof(CKmer<SIZE>));
	input_pos += size;
}

template<unsigned SIZE> uint64 CKmerBinSorter<SIZE>::ExpandKxmerBothParallel(uint64 start_pos, uint64 end_pos, uint64 output_start, uint64 output_end)
{
	CKmer<SIZE> kmer, rev_kmer, kmer_mask;
	CKmer<SIZE>  kxmer_mask;
	bool kmer_lower; //true if kmer is lower than its rev. comp
	uint32 x, additional_symbols;
	uchar symb;
	uint32 kmer_bytes = (kmer_len + 3) / 4;
	uint32 rev_shift = kmer_len * 2 - 2;
	uchar *data_p = data;
	kmer_mask.set_n_1(kmer_len * 2);
	uint32 kmer_shr = SIZE * 32 - kmer_len;

	kxmer_mask.set_n_1((kmer_len + max_x + 1) * 2);

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
		uchar byte_shift = 6 - (kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		rev_kmer.mask(kmer_mask);

		kmer_lower = kmer < rev_kmer;
		x = 0;
		if (kmer_lower)
			buffer_input[output_start].set(kmer);
		else
			buffer_input[output_start].set(rev_kmer);

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
					buffer_input[output_start].SHL_insert_2bits(symb);
					++x;
					if (x == max_x)
					{
						if (!symbols_left)
							break;

						buffer_input[output_start++].set_2bits(x, kmer_len * 2 + max_x * 2);
						
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;

						if (kmer_lower)
							buffer_input[output_start].set(kmer);
						else
							buffer_input[output_start].set(rev_kmer);
					}
				}
				else
				{
					buffer_input[output_start++].set_2bits(x, kmer_len * 2 + max_x * 2);
					
					x = 0;

					kmer_lower = false;
					buffer_input[output_start].set(rev_kmer);

				}
			}
			else
			{
				if (!(kmer < rev_kmer))
				{
					buffer_input[output_start].set_2bits(3 - symb, kmer_len * 2 + x * 2);
					++x;
					if (x == max_x)
					{
						if (!symbols_left)
							break;

						buffer_input[output_start++].set_2bits(x, kmer_len * 2 + max_x * 2);
						
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;

						if (kmer_lower)
							buffer_input[output_start].set(kmer);
						else
							buffer_input[output_start].set(rev_kmer);
					}
				}
				else
				{
					buffer_input[output_start++].set_2bits(x, kmer_len * 2 + max_x * 2);
					
					x = 0;

					buffer_input[output_start].set(kmer);
					kmer_lower = true;
				}
			}

		}
		buffer_input[output_start++].set_2bits(x, kmer_len * 2 + max_x * 2);

		if (byte_shift != 6)
			++pos;
	}	
	
	uint64 ret = output_end - output_start;
	/*for (; output_start < output_end; ++output_start)
		buffer_input[output_start].fill_T();*/
	return ret;
}


template<unsigned SIZE>
class CExpandThread
{
	CKmerBinSorter<SIZE>& ptr;
	CExpanderPackQueue& q;

	list<pair<uint64, uint64>> filled_regions;

	uint64 fake_recs = 0;

public:
	CExpandThread(CKmerBinSorter<SIZE>& ptr, CExpanderPackQueue& q)
		:
		ptr(ptr),
		q(q)
	{		
	}
	void operator()()
	{
		uint64 start, end, out_start, out_end;
		while (q.Pop(start, end, out_start, out_end))
		{

			auto res = ptr.ExpandKxmerBothParallel(start, end, out_start, out_end);
			fake_recs += res;
			filled_regions.push_back(make_pair(out_start, out_end - res));
		}
	}


	uint64 GetNFakeRecs()
	{
		return fake_recs;
	}
	list<pair<uint64, uint64>>& GetFilledRegions()
	{
		return filled_regions;
	}

};

template <unsigned SIZE> void CKmerBinSorter<SIZE>::ExpandKxmersBoth(uint64 tmp_size)
{	
	if (!tmp_size)
	{		
		return;
	}
	input_pos = 0;
	uint32 threads = n_sorting_threads;

	list<pair<uint64, uint64>> l;
	epd->pop(bin_id, l);

	if (n_sorting_threads > 1)	
	{	
		CExpanderPackQueue q(l);

		vector<thread> ths;
		vector<CExpandThread<SIZE>*> exp;
		for (uint32 i = 0; i < threads; ++i)
		{
			exp.push_back(new CExpandThread<SIZE>(*this, q));
			ths.push_back(thread(std::ref(*exp.back())));
		}
		
		for (auto& t : ths)
			t.join();

		uint64 n_fake_recs_after_expand = 0;
		vector<pair<uint64, uint64>> filled_regions;
		uint64 filled_regions_size = 0;
		for (auto _ptr : exp)
		{
			filled_regions_size += _ptr->GetFilledRegions().size();
		}
		filled_regions.reserve(filled_regions_size);
		
		for (auto _ptr : exp)
		{
			n_fake_recs_after_expand += _ptr->GetNFakeRecs();
			filled_regions.insert(filled_regions.end(), _ptr->GetFilledRegions().begin(), _ptr->GetFilledRegions().end());
			delete _ptr;
		}


		//remove fakes
		sort(filled_regions.begin(), filled_regions.end(), [](const pair<uint64, uint64>& e1, const pair<uint64, uint64>& e2) {return e1.first < e2.first; });
		
		uint64 first_gap_pos = filled_regions.front().second;
		uint64 next_reg_desc_pos = 1;
		uint64 last_elem_pos = filled_regions.back().second - 1;
		int64 back_reg_desc_pos = filled_regions.size() - 1;

		while (true)
		{
			while (next_reg_desc_pos < filled_regions.size() && first_gap_pos >= filled_regions[next_reg_desc_pos].first)
				first_gap_pos = filled_regions[next_reg_desc_pos++].second;

			while (back_reg_desc_pos >= 0 && last_elem_pos < filled_regions[back_reg_desc_pos].first)
				last_elem_pos = filled_regions[--back_reg_desc_pos].second - 1;

			if (first_gap_pos >= last_elem_pos)
				break;

			buffer_input[first_gap_pos++] = buffer_input[last_elem_pos--];
		}
		
		n_plus_x_recs -= n_fake_recs_after_expand;				
	}
	else
	{
		l.clear();
		auto n_fake_recs = ExpandKxmerBothParallel(0, tmp_size, 0, n_plus_x_recs);
		n_plus_x_recs -= n_fake_recs;
		
	}
	
	input_pos = n_plus_x_recs;
}

template<unsigned SIZE> void CKmerBinSorter<SIZE>::ExpandKxmersAll(uint64 tmp_size)
{
	input_pos = 0;
	uint64 pos = 0;
	CKmer<SIZE> kmer_mask;

	CKmer<SIZE> kxmer;
	CKmer<SIZE> kxmer_mask;
	kxmer_mask.set_n_1((kmer_len + max_x) * 2);
	uchar *data_p = data;

	kmer_mask.set_n_1(kmer_len * 2);
	while (pos < tmp_size)
	{
		kxmer.clear();
		uint32 additional_symbols = data_p[pos++];

		uchar symb;

		uint32 kmer_bytes = (kmer_len + 3) / 4;
		//building kmer
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1; i < kmer_bytes; ++i, --kmer_pos)
		{
			kxmer.set_byte(kmer_pos, data_p[pos + i]);
		}

		pos += kmer_bytes;
		uchar byte_shift = 6 - (kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;
		uint32 kmer_shr = SIZE * 32 - kmer_len;

		if (kmer_shr)
			kxmer.SHR(kmer_shr);

		kxmer.mask(kmer_mask);
		uint32 tmp = MIN(max_x, additional_symbols);

		for (uint32 i = 0; i < tmp; ++i)
		{
			GetNextSymb(symb, byte_shift, pos, data_p);
			kxmer.SHL_insert_2bits(symb);
		}
		kxmer.set_2bits(tmp, (kmer_len + max_x) * 2);

		buffer_input[input_pos++].set(kxmer);
		additional_symbols -= tmp;

		uint32 kxmers_count = additional_symbols / (max_x + 1);
		uint32 kxmer_rest = additional_symbols % (max_x + 1);

		for (uint32 j = 0; j < kxmers_count; ++j)
		{
			for (uint32 i = 0; i < max_x + 1; ++i)
			{
				GetNextSymb(symb, byte_shift, pos, data_p);
				kxmer.SHL_insert_2bits(symb);
			}

			kxmer.mask(kxmer_mask);

			kxmer.set_2bits(max_x, (kmer_len + max_x) * 2);

			buffer_input[input_pos++].set(kxmer);
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

			kxmer.set_2bits(kxmer_rest, (kmer_len + max_x) * 2);
			buffer_input[input_pos++].set(kxmer);
		}
		if (byte_shift != 6)
			++pos;
	}
}

//----------------------------------------------------------------------------------
// Uncompact the kmers
template <unsigned SIZE> void CKmerBinSorter<SIZE>::Expand(uint64 tmp_size)
{
	uchar *raw_buffer_input, *raw_buffer_tmp;

	memory_bins->reserve(bin_id, raw_buffer_input, CMemoryBins::mba_input_array);
	memory_bins->reserve(bin_id, raw_buffer_tmp, CMemoryBins::mba_tmp_array);

	buffer_input = (CKmer<SIZE> *) raw_buffer_input;
	buffer_tmp = (CKmer<SIZE> *) raw_buffer_tmp;

	if (max_x)
	{
		if (both_strands)
			ExpandKxmersBoth(tmp_size);
		else
			ExpandKxmersAll(tmp_size);
	}
	else
	{
		if (both_strands)
			ExpandKmersBoth(tmp_size);
		else
			ExpandKmersAll(tmp_size);
	}
}


//----------------------------------------------------------------------------------
// Sort the kmers
template <unsigned SIZE> void CKmerBinSorter<SIZE>::Sort()
{
	uint32 rec_len;
	uint64 sort_rec;
	if (max_x)
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
	
	sort_func(buffer_input, buffer_tmp, sort_rec, rec_len - 1, n_sorting_threads, pmm_radix_buf);
	if (rec_len % 2)
		buffer = buffer_tmp;
	else
		buffer = buffer_input;	
}

//----------------------------------------------------------------------------------
//Binary search position of first occurrence of symbol 'symb' in [start_pos,end_pos). Offset defines which symbol in k+x-mer is taken.
template <unsigned SIZE> uint64 CKmerBinSorter<SIZE>::FindFirstSymbOccur(uint64 start_pos, uint64 end_pos, uint32 offset, uchar symb)
{
	uint32 kxmer_offset = (kmer_len + max_x - offset) * 2;
	uint64 middle_pos;
	uchar middle_symb;
	while (start_pos < end_pos)
	{
		middle_pos = (start_pos + end_pos) / 2;
		middle_symb = buffer[middle_pos].get_2bits(kxmer_offset);
		if (middle_symb < symb)
			start_pos = middle_pos + 1;
		else
			end_pos = middle_pos;
	}
	return end_pos;
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CKmerBinSorter<SIZE>::InitKXMerSet(uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth)
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
template<unsigned SIZE> void CKmerBinSorter<SIZE>::InitKXMerSetMultithreaded(CKXmerSetMultiThreaded<SIZE>& kxmer_set_multithreaded, uint64 start_pos, uint64 end_pos, uint32 offset, uint32 depth)
{
	if (start_pos == end_pos)
		return;
	uint32 shr = max_x + 1 - offset;
	kxmer_set_multithreaded.InitAdd(start_pos, end_pos, shr);

	--depth;
	if (depth > 0)
	{
		uint64 pos[5];
		pos[0] = start_pos;
		pos[4] = end_pos;
		for (uint32 i = 1; i < 4; ++i)
			pos[i] = FindFirstSymbOccur(pos[i - 1], end_pos, offset, i);
		for (uint32 i = 1; i < 5; ++i)
			InitKXMerSetMultithreaded(kxmer_set_multithreaded, pos[i - 1], pos[i], offset + 1, depth);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CKmerBinSorter<SIZE>::PreCompactKxmers(uint64& compacted_count)
{
	uint32 n_threads = n_sorting_threads;
	vector<thread> threads;
	vector<pair<uint64, uint64>> start_end(n_threads);
	uint64 total_recs = n_plus_x_recs;		
	for (uint32 idx = 0; idx < n_threads; ++idx)
	{
		threads.push_back(thread([idx, n_threads, total_recs, &start_end, this]
		{
			uint64 per_thread = total_recs / n_threads;
			uint64 start = idx * per_thread;
			uint64 end = (idx + 1) * per_thread;
			start_end[idx].first = start;
			if (idx == n_threads - 1)
				end = total_recs;

			if (start == end)
			{
				start_end[idx].second = end;
				return;
			}

			uint64 compacted_pos = start;
			CKmer<SIZE>* act_kmer = &buffer[compacted_pos];
			kxmer_counters[compacted_pos] = 1;

			for (uint64 i = start + 1; i < end; ++i)
			{
				if (*act_kmer == buffer[i])
					++kxmer_counters[compacted_pos];
				else
				{
					buffer[compacted_pos++] = *act_kmer;
					kxmer_counters[compacted_pos] = 1;
					act_kmer = &buffer[i];
				}
			}
			buffer[compacted_pos++] = *act_kmer;
			start_end[idx].second = compacted_pos;
		}));
	}

	for (auto& t : threads)
		t.join();

	compacted_count = start_end[0].second;
	for (uint32 i = 1; i < n_threads; ++i)
	{
		if (start_end[i].second > start_end[i].first)
		{
			if (buffer[compacted_count - 1] == buffer[start_end[i].first])
			{
				kxmer_counters[compacted_count - 1] += kxmer_counters[start_end[i].first++];
			}
			uint64 n_elems = start_end[i].second - start_end[i].first;
			if (!n_elems)
				continue;
			auto kmers_dest = buffer + compacted_count;
			auto kmers_src = buffer + start_end[i].first;
			auto counters_dest = kxmer_counters + compacted_count;
			auto counters_src = kxmer_counters + start_end[i].first;

			if (n_threads < 2 || n_elems < 4096) //size limit achieved experimentally
			{
				memmove(kmers_dest, kmers_src, n_elems * sizeof(CKmer<SIZE>));
				memmove(counters_dest, counters_src, n_elems * sizeof(uint32));
			}
			else
			{
				std::thread th1([kmers_dest, kmers_src, n_elems]
				{
					auto _kmers_dest = kmers_dest;
					auto _kmers_src = kmers_src;
					memmove(_kmers_dest, _kmers_src, n_elems * sizeof(CKmer<SIZE>));
				});

				std::thread th2([counters_dest, counters_src, n_elems]
				{
					auto _counters_dest = counters_dest;
					auto _counters_src = counters_src;
					memmove(_counters_dest, _counters_src, n_elems * sizeof(uint32));
				});
				th1.join();
				th2.join();
			}

			compacted_count += n_elems;
		}
	}
}

//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKmerBinSorter<SIZE>::CompactKxmers()
{
	kxmer_set.clear();
	kxmer_set.set_buffer(buffer);
	n_unique = 0;
	n_cutoff_min = 0;
	n_cutoff_max = 0;
	n_total = 0;

	uint32 kmer_symbols = kmer_len - lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;
	uint64 lut_recs = 1ull << (2 * lut_prefix_len);
	uint64 lut_size = lut_recs * sizeof(uint64);


	uchar *out_buffer = nullptr;
	uchar *raw_lut = nullptr;

	memory_bins->reserve(bin_id, out_buffer, CMemoryBins::mba_suffix);
	memory_bins->reserve(bin_id, raw_lut, CMemoryBins::mba_lut);

	uint64 *lut = (uint64*)raw_lut;
	fill_n(lut, lut_recs, 0);

	list<pair<uint64, uint64>> output_packs_desc;
	if (n_plus_x_recs)
	{
		uchar* raw_kxmer_counters = nullptr;
		memory_bins->reserve(bin_id, raw_kxmer_counters, CMemoryBins::mba_kxmer_counters);
		kxmer_counters = (uint32*)raw_kxmer_counters;
		uint64 compacted_count;
		PreCompactKxmers(compacted_count);
		
		uint64 pos[5];//pos[symb] is first position where symb occur (at first position of k+x-mer) and pos[symb+1] is first position where symb is not starting symbol of k+x-mer
		pos[0] = 0;
		pos[4] = compacted_count;
		for (uint32 i = 1; i < 4; ++i)
			pos[i] = FindFirstSymbOccur(pos[i - 1], compacted_count, 0, i);

		if (n_sorting_threads > 1)
		{			
			CKXmerSetMultiThreaded<SIZE> kxmer_set_multithreaded(buffer, kxmer_counters, compacted_count, 
				cutoff_min, cutoff_max, counter_max, kmer_len, lut_prefix_len, lut, out_buffer, n_sorting_threads);

			for (uint32 i = 1; i < 5; ++i)
				InitKXMerSetMultithreaded(kxmer_set_multithreaded, pos[i - 1], pos[i], max_x + 2 - i, i);

			kxmer_set_multithreaded.Process(without_output);
		
			kxmer_set_multithreaded.GetStats(n_unique, n_cutoff_min, n_cutoff_max, n_total);	
			output_packs_desc = std::move(kxmer_set_multithreaded.GetOutputPacksDesc());			
		}
		else
		{
			uint64 out_pos = 0; 
			for (uint32 i = 1; i < 5; ++i)
				InitKXMerSet(pos[i - 1], pos[i], max_x + 2 - i, i);

			uint64 counter_pos = 0;
			uint64 counter_size = min(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));

			CKmer<SIZE> kmer, next_kmer;
			kmer.clear();
			next_kmer.clear();
			CKmer<SIZE> kmer_mask;
			kmer_mask.set_n_1(kmer_len * 2);
			uint32 count;
			//first
			kxmer_set.get_min(counter_pos, kmer);
			count = kxmer_counters[counter_pos];
			//rest
			while (kxmer_set.get_min(counter_pos, next_kmer))
			{
				if (kmer == next_kmer)
					count += kxmer_counters[counter_pos];
				else
				{
					n_total += count;
					++n_unique;
					if (count < cutoff_min)
						n_cutoff_min++;
					else if (count >cutoff_max)
						n_cutoff_max++;
					else
					{
						if (count > counter_max)
							count = counter_max;

						if (!without_output)
						{
							lut[kmer.remove_suffix(2 * kmer_symbols)]++;
							// Store compacted kmer

							for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
								out_buffer[out_pos++] = kmer.get_byte(j);
							for (int32 j = 0; j < (int32)counter_size; ++j)
								out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;
						}
					}
					count = kxmer_counters[counter_pos];
					kmer = next_kmer;
				}
			}

			//last one
			++n_unique;
			n_total += count;
			if (count < cutoff_min)
				n_cutoff_min++;
			else if (count >cutoff_max)
				n_cutoff_max++;
			else
			{
				
				if (count > counter_max)
					count = counter_max;

				if (!without_output)
				{
					lut[kmer.remove_suffix(2 * kmer_symbols)]++;
					// Store compacted kmer

					for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
						out_buffer[out_pos++] = kmer.get_byte(j);
					for (int32 j = 0; j < (int32)counter_size; ++j)
						out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;
				}
			}

			if(!without_output)
				output_packs_desc.emplace_back(0, out_pos);
		}

		memory_bins->free(bin_id, CMemoryBins::mba_kxmer_counters);
	}


	// Push the sorted and compacted kmer bin to a queue in a form ready to be stored to HDD	
	kq->push(bin_id, out_buffer, output_packs_desc, raw_lut, lut_size, n_unique, n_cutoff_min, n_cutoff_max, n_total);

	if (buffer_input)
	{
		memory_bins->free(bin_id, CMemoryBins::mba_input_array);
		memory_bins->free(bin_id, CMemoryBins::mba_tmp_array);
	}
	buffer = nullptr;
}




//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKmerBinSorter<SIZE>::CompactKmers()
{
	uint64 i;

	uint32 kmer_symbols = kmer_len - lut_prefix_len;
	uint64 kmer_bytes = kmer_symbols / 4;
	uint64 lut_recs = 1ull << (2 * (lut_prefix_len));
	uint64 lut_size = lut_recs * sizeof(uint64);

	uint64 counter_size = min(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));

	uchar *out_buffer;
	uchar *raw_lut;

	memory_bins->reserve(bin_id, out_buffer, CMemoryBins::mba_suffix);
	memory_bins->reserve(bin_id, raw_lut, CMemoryBins::mba_lut);
	uint64 *lut = (uint64*)raw_lut;
	fill_n(lut, lut_recs, 0);

	uint64 out_pos = 0;
	uint32 count;
	CKmer<SIZE> *act_kmer;
	
	n_unique = 0;
	n_cutoff_min = 0;
	n_cutoff_max = 0;
	n_total = 0;

	if (n_rec)			// non-empty bin
	{
		act_kmer = &buffer[0];
		count = 1;

		n_total = n_rec;

		for (i = 1; i < n_rec; ++i)
		{
			if (*act_kmer == buffer[i])
				count++;
			else
			{
				if (count < cutoff_min)
				{
					act_kmer = &buffer[i];
					n_cutoff_min++;
					n_unique++;
					count = 1;
				}
				else if (count > cutoff_max)
				{
					act_kmer = &buffer[i];
					n_cutoff_max++;
					n_unique++;
					count = 1;
				}
				else
				{
					if (count > counter_max)
						count = counter_max;

					if (!without_output)
					{
						// Store compacted kmer
						for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
							out_buffer[out_pos++] = act_kmer->get_byte(j);
						for (int32 j = 0; j < (int32)counter_size; ++j)
							out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;

						lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;
					}
					act_kmer = &buffer[i];
					count = 1;
					n_unique++;
				}
			}
		}

		if (count < cutoff_min)
		{
			n_cutoff_min++;
		}
		else if (count >= cutoff_max)
		{
			n_cutoff_max++;
		}
		else
		{
			if (count > counter_max)
				count = counter_max;

			if (!without_output)
			{
				for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
					out_buffer[out_pos++] = act_kmer->get_byte(j);
				for (int32 j = 0; j < (int32)counter_size; ++j)
					out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;
				lut[act_kmer->remove_suffix(2 * kmer_symbols)]++;
			}
		}
		n_unique++;
	}
	list<pair<uint64, uint64>> data_packs;
	if(!without_output)
		data_packs.emplace_back(0, out_pos);
	// Push the sorted and compacted kmer bin to a priority queue in a form ready to be stored to HDD
	kq->push(bin_id, out_buffer, data_packs, raw_lut, lut_size, n_unique, n_cutoff_min, n_cutoff_max, n_total);

	if (buffer_input)
	{
		memory_bins->free(bin_id, CMemoryBins::mba_input_array);
		memory_bins->free(bin_id, CMemoryBins::mba_tmp_array);
	}
	buffer = nullptr;
}




//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKmerBinSorter<SIZE>::Compact()
{
	if (max_x)
		CompactKxmers();
	else
		CompactKmers();
}

//************************************************************************************************************
// CWKmerBinSorter - wrapper for multithreading purposes
//************************************************************************************************************
template <unsigned SIZE> class CWKmerBinSorter {
	CKmerBinSorter<SIZE> *kbs;

public:	
	CWKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, SortFunction<CKmer<SIZE>> sort_func);
	~CWKmerBinSorter();
	void GetDebugStats(uint64& _sum_n_recs, uint64& _sum_n_plus_x_recs)
	{
		kbs->GetDebugStats(_sum_n_recs, _sum_n_plus_x_recs);
	}
	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template <unsigned SIZE> CWKmerBinSorter<SIZE>::CWKmerBinSorter(CKMCParams &Params, CKMCQueues &Queues, SortFunction<CKmer<SIZE>> sort_func)
{
	kbs = new CKmerBinSorter<SIZE>(Params, Queues, sort_func);
}

//----------------------------------------------------------------------------------
// Destructor
template <unsigned SIZE> CWKmerBinSorter<SIZE>::~CWKmerBinSorter()
{
	delete kbs;
}

//----------------------------------------------------------------------------------
// Execution
template <unsigned SIZE> void CWKmerBinSorter<SIZE>::operator()()
{
	kbs->ProcessBins();
}


#endif
// ***** EOF
