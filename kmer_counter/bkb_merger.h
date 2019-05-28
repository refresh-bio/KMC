/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _HBH_MERGER_H
#define _HBH_MERGER_H
#include "bkb_subbin.h"


//************************************************************************************************************
// CBigKmerBinMerger - merger sorted k-mers from number of subbins 
//************************************************************************************************************
template<unsigned SIZE>
class CBigKmerBinMerger
{
	vector<CSubBin<SIZE>*> sub_bins;
	std::vector<std::tuple<CKmer<SIZE>, uint32, uint32>> curr_min;
	CDiskLogger* disk_logger;
	uint32 size;
	CBigBinDesc* bbd;
	CBigBinKmerPartQueue* bbkpq;
	CCompletedBinsCollector* sm_cbc;
	uint32 kmer_len;
	uint32 lut_prefix_len;	
	uint32 cutoff_min, cutoff_max, counter_max;
	CMemoryPool* sm_pmm_merger_suff, *sm_pmm_merger_lut, *sm_pmm_sub_bin_suff, *sm_pmm_sub_bin_lut;
	int64 sm_mem_part_merger_suff, sm_mem_part_merger_lut, sm_mem_part_sub_bin_suff, sm_mem_part_sub_bin_lut;
	uchar *sub_bin_suff_buff, *sub_bin_lut_buff;
public:
	CBigKmerBinMerger(CKMCParams& Params, CKMCQueues& Queues);
	void init(int32 bin_id, uint32 _size);
	bool get_min(CKmer<SIZE>& kmer, uint32& count);
	void Process();
	~CBigKmerBinMerger();
};

//----------------------------------------------------------------------------------
template<unsigned SIZE>
CBigKmerBinMerger<SIZE>::CBigKmerBinMerger(CKMCParams& Params, CKMCQueues& Queues) 
{
	disk_logger = Queues.disk_logger;
	bbd = Queues.bbd;
	bbkpq = Queues.bbkpq;
	sm_cbc = Queues.sm_cbc;
	kmer_len = Params.kmer_len;
	lut_prefix_len = Params.lut_prefix_len;
	cutoff_min = Params.cutoff_min;
	cutoff_max = (uint32)Params.cutoff_max;
	counter_max = (uint32)Params.counter_max;
	sm_pmm_merger_suff = Queues.sm_pmm_merger_suff;
	sm_pmm_merger_lut = Queues.sm_pmm_merger_lut;
	sm_pmm_sub_bin_suff = Queues.sm_pmm_sub_bin_suff;
	sm_pmm_sub_bin_lut = Queues.sm_pmm_sub_bin_lut;
	sm_mem_part_sub_bin_suff = Params.sm_mem_part_sub_bin_suff;
	sm_mem_part_merger_suff = Params.sm_mem_part_merger_suff;
	sm_mem_part_merger_lut = Params.sm_mem_part_merger_lut;
	sm_mem_part_sub_bin_lut = Params.sm_mem_part_sub_bin_lut;

	sm_pmm_sub_bin_lut->reserve(sub_bin_lut_buff);
	sm_pmm_sub_bin_suff->reserve(sub_bin_suff_buff);
}

//----------------------------------------------------------------------------------
template<unsigned SIZE>
CBigKmerBinMerger<SIZE>::~CBigKmerBinMerger()
{
	for (auto p : sub_bins)
		delete p;
	sm_pmm_sub_bin_lut->free(sub_bin_lut_buff);
	sm_pmm_sub_bin_suff->free(sub_bin_suff_buff);
}

//----------------------------------------------------------------------------------
template<unsigned SIZE>
void CBigKmerBinMerger<SIZE>::init(int32 bin_id, uint32 _size)
{
	size = _size;
	uint32 prev_size = (uint32)sub_bins.size();
	int32 sub_bin_id;
	if (size > prev_size)
	{
		sub_bins.resize(size);
		curr_min.resize(size);
		for (uint32 i = prev_size; i < size; ++i)
		{
			sub_bins[i] = new CSubBin<SIZE>(disk_logger);
		}
	}

	uint32 lut_prefix_len = 0;;
	uint64 n_kmers = 0;
	uint64 file_size = 0;
	FILE* file = nullptr;
	string name;
	uint32 per_sub_bin_lut_size = (uint32)(sm_mem_part_sub_bin_lut / size);
	uint32 per_sub_bin_suff_size = (uint32)(sm_mem_part_sub_bin_suff / size);
	for (uint32 i = 0; i < size; ++i)
	{
		bbd->next_sub_bin(bin_id, sub_bin_id, lut_prefix_len, n_kmers, file, name, file_size);
		sub_bins[i]->init(file, file_size, lut_prefix_len, n_kmers, name, kmer_len, sub_bin_lut_buff + i * per_sub_bin_lut_size, per_sub_bin_lut_size, sub_bin_suff_buff + i * per_sub_bin_suff_size, per_sub_bin_suff_size);
		get<2>(curr_min[i]) = i;
		sub_bins[i]->get_min(get<0>(curr_min[i]), get<1>(curr_min[i]));
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE>
bool CBigKmerBinMerger<SIZE>::get_min(CKmer<SIZE>& kmer, uint32& count)
{
	if (!size)
		return false;
	uint32 min = 0;
	for (uint32 i = 1; i < size; ++i)
		if (get<0>(curr_min[i]) < get<0>(curr_min[min]))
			min = i;

	kmer = get<0>(curr_min[min]);
	count = get<1>(curr_min[min]);
	if (sub_bins[get<2>(curr_min[min])]->get_min(get<0>(curr_min[min]), get<1>(curr_min[min])))
		;
	else
		curr_min[min] = curr_min[--size];
	return true;
}

//----------------------------------------------------------------------------------
template<unsigned SIZE>
void CBigKmerBinMerger<SIZE>::Process()
{
	int32 bin_id;
	uint32 size = 0;
	uint32 counter_size = min(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));
	uint32 lut_recs = 1 << 2 * lut_prefix_len;
	uint32 kmer_symbols = (kmer_len - lut_prefix_len);
	uint32 kmer_bytes = kmer_symbols / 4;
	uint32 suff_rec_bytes = kmer_bytes + counter_size;
	uint64 suff_buff_size = sm_mem_part_merger_suff / suff_rec_bytes * suff_rec_bytes;
	uint64 suff_buff_pos = 0;
	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total;
	CKmer<SIZE> kmer, next_kmer;
	kmer.clear();
	next_kmer.clear();
	uint32 count_tmp = 0, count = 0;
	int32 max_in_lut = (int32)(sm_mem_part_merger_lut / sizeof(uint64));

	while (sm_cbc->pop(bin_id))
	{
		bbd->get_n_sub_bins(bin_id, size);
		uchar *raw_lut;
		sm_pmm_merger_lut->reserve(raw_lut);
		uint64 *lut = (uint64*)raw_lut;
		uchar* suff_buff;
		sm_pmm_merger_suff->reserve(suff_buff);
		suff_buff_pos = 0;
		n_unique = n_cutoff_min = n_cutoff_max = n_total = 0;
		fill_n(lut, max_in_lut, 0);
		init(bin_id, size);

		get_min(kmer, count_tmp);
		count = count_tmp;
		uint32 lut_offset = 0;
		uint64 prefix;
		while (get_min(next_kmer, count_tmp))
		{
			if (kmer == next_kmer)
				count += count_tmp;
			else
			{
				++n_unique;
				n_total += count;
				if (count < cutoff_min)
					n_cutoff_min++;
				else if (count > cutoff_max)
					n_cutoff_max++;
				else
				{
					if (count > counter_max)
						count = counter_max;

					//store
					prefix = kmer.remove_suffix(2 * kmer_symbols);
					if (prefix >= max_in_lut + lut_offset)
					{
						bbkpq->push(bin_id, nullptr, 0, raw_lut, max_in_lut * sizeof(uint64), 0, 0, 0, 0, false);
						lut_offset += max_in_lut;
						sm_pmm_merger_lut->reserve(raw_lut);
						lut = (uint64*)raw_lut;
						fill_n(lut, max_in_lut, 0);
					}

					lut[prefix - lut_offset]++;

					for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
						suff_buff[suff_buff_pos++] = kmer.get_byte(j);
					for (int32 j = 0; j < (int32)counter_size; ++j)
						suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;

					if (suff_buff_pos >= suff_buff_size)
					{
						bbkpq->push(bin_id, suff_buff, suff_buff_pos, nullptr, 0, 0, 0, 0, 0, false);
						suff_buff_pos = 0;
						sm_pmm_merger_suff->reserve(suff_buff);
					}
				}
				count = count_tmp;
				kmer = next_kmer;
			}
		}
		++n_unique;
		n_total += count;
		if (count < cutoff_min)
			++n_cutoff_min;
		else if (count > cutoff_max)
			++n_cutoff_max;
		else
		{
			if (count > counter_max)
				count = counter_max;

			//store
			lut[kmer.remove_suffix(2 * kmer_symbols)]++;

			for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
				suff_buff[suff_buff_pos++] = kmer.get_byte(j);
			for (int32 j = 0; j < (int32)counter_size; ++j)
				suff_buff[suff_buff_pos++] = (count >> (j * 8)) & 0xFF;
		}
		bbkpq->push(bin_id, suff_buff, suff_buff_pos, raw_lut, (lut_recs - lut_offset) * sizeof(uint64), n_unique, n_cutoff_min, n_cutoff_max, n_total, true);
	}

	bbkpq->mark_completed();
}


//************************************************************************************************************
// CWBigKmerBinMerger - wrapper for multithreading purposes
//************************************************************************************************************
template<unsigned SIZE>
class CWBigKmerBinMerger
{
	CBigKmerBinMerger<SIZE> *merger;
public:
	CWBigKmerBinMerger(CKMCParams& Params, CKMCQueues& Queues);
	~CWBigKmerBinMerger();
	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template<unsigned SIZE>
CWBigKmerBinMerger<SIZE>::CWBigKmerBinMerger(CKMCParams& Params, CKMCQueues& Queues)
{
	merger = new CBigKmerBinMerger<SIZE>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
template<unsigned SIZE>
CWBigKmerBinMerger<SIZE>::~CWBigKmerBinMerger()
{
	delete merger;
}

#endif

//----------------------------------------------------------------------------------
// Execution
template<unsigned SIZE>
void CWBigKmerBinMerger<SIZE>::operator()()
{
	merger->Process();
}

// ***** EOF 