/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/
#ifndef _KXMER_SET_
#define _KXMER_SET_
#include "defs.h"
#include <tuple>
#include <queue>

using namespace std;

#define KXMER_SET_SIZE 1024 

#define MAX_FOR_X_3 112

template <unsigned SIZE>
class CKXmerSet;

template<bool> struct Range;

template<unsigned SIZE, unsigned PARENT, typename = Range<true> >
struct ParentFinder
{
	static uint32 Execute(CKXmerSet<SIZE>& ptr, const CKmer<SIZE>& kmer);
};

template<unsigned SIZE, unsigned PARENT>
struct ParentFinder<SIZE, PARENT, Range<(PARENT < 56)> >
{
	FORCE_INLINE static uint32 Execute(CKXmerSet<SIZE>& ptr, const CKmer<SIZE>& kmer)
	{
		if (ptr.data[PARENT * 2].first < ptr.data[PARENT * 2 + 1].first)
		{
			if (ptr.data[PARENT * 2].first < kmer)
			{
				ptr.data[PARENT] = ptr.data[PARENT * 2];
				return ParentFinder<SIZE, PARENT * 2>::Execute(ptr, kmer);
			}
		}
		else
		{
			if (ptr.data[PARENT * 2 + 1].first < kmer)
			{
				ptr.data[PARENT] = ptr.data[PARENT * 2 + 1];
				return ParentFinder<SIZE, PARENT * 2 + 1>::Execute(ptr, kmer);
			}
		}
		return PARENT;
	}
};

template<unsigned SIZE, unsigned PARENT>
struct ParentFinder<SIZE, PARENT, Range<(PARENT == 56)> >
{
	FORCE_INLINE static uint32 Execute(CKXmerSet<SIZE>& ptr, const CKmer<SIZE>& kmer)
	{
		if (ptr.data[PARENT * 2].first < kmer)
		{
			ptr.data[PARENT] = ptr.data[PARENT * 2];
			return PARENT * 2;
		}
		return PARENT;
	}
};

template<unsigned SIZE, unsigned PARENT>
struct ParentFinder<SIZE, PARENT, Range<(PARENT > 56)> >
{
	FORCE_INLINE static uint32 Execute(CKXmerSet<SIZE>& /*ptr*/, const CKmer<SIZE>& /*kmer*/)
	{		
		return PARENT;		
	}
};

template <unsigned SIZE>
class CKXmerSet
{
	typedef tuple<uint64, uint64, uint32> elem_desc_t; //start_pos, end_pos, shr
	typedef pair<CKmer<SIZE>, uint32> heap_elem_t; //kxmer val, desc_id
	elem_desc_t data_desc[KXMER_SET_SIZE];
	
public: 
	heap_elem_t data[KXMER_SET_SIZE];
private:
	uint32 pos;
	uint32 desc_pos;
	CKmer<SIZE> mask;
	CKmer<SIZE>* buffer;

	inline void update_heap()
	{
		uint32 desc_id = data[1].second;
		CKmer<SIZE> kmer;
		if (++get<0>(data_desc[desc_id]) < get<1>(data_desc[desc_id]))
		{
			kmer.from_kxmer(buffer[get<0>(data_desc[desc_id])], get<2>(data_desc[desc_id]), mask);
		}
		else
		{
			kmer.set(data[--pos].first);
			desc_id = data[pos].second;
			data[pos].first.fill_T();		
		}

		uint32 parent = ParentFinder<SIZE, 1>::Execute(*this, kmer);

		data[parent] = make_pair(kmer, desc_id);
	}

public:
	CKXmerSet(uint32 kmer_len)
	{
		pos = 1;
		mask.set_n_1(kmer_len * 2);
		desc_pos = 0;		
	}	

	inline void init_add(uint64 start_pos, uint64 end_pos, uint32 shr)
	{		
		data_desc[desc_pos] = make_tuple(start_pos, end_pos, shr);
		data[pos].first.from_kxmer(buffer[start_pos], shr, mask);
		data[pos].second = desc_pos;
		uint32 child_pos = pos++;

		while (child_pos > 1 && data[child_pos].first < data[child_pos / 2].first)
		{
			swap(data[child_pos], data[child_pos / 2]);
			child_pos /= 2;
		}
		++desc_pos;
	}
	inline void set_buffer(CKmer<SIZE>* _buffer)
	{
		buffer = _buffer;
	}
	inline void clear()
	{
		pos = 1;
		desc_pos = 0;
		for (uint32 i = 0; i < KXMER_SET_SIZE; ++i)
		{
			data[i].first.fill_T();
			data[i].second = -1;
		}
	}

	inline bool get_min(uint64& _pos, CKmer<SIZE>& kmer)
	{
		if (pos <= 1)		
			return false;				
		kmer = data[1].first;
		_pos = get<0>(data_desc[data[1].second]);
		update_heap();
		
		return true;
	}
};



struct SubArrayDesc
{
	uint64 start, end;
	uint32 shr;
	uint64 counters_sum;
};

template<unsigned SIZE>
class CSubArrayDescGenerator
{
	uint32 kmer_len;
	uint32 n_parts;
	uint32 parts_left;
	const vector<SubArrayDesc>& sub_array_descs;
	queue<vector<SubArrayDesc>> data;
	CKmer<SIZE>* buffer;
	mutable std::mutex mtx;	
	uint64 out_start = 0;
	uint32 cutoff_min;
	uint32 rec_len;
	vector<uint64> cumsum;
	uint32* kxmer_counters;
	uint64 n_kxmer_counters;
	uint32 n_threads;
public:	
	CSubArrayDescGenerator(uint32 kmer_len, uint32 n_parts, const vector<SubArrayDesc>& sub_array_descs, CKmer<SIZE>* buffer, uint32 cutoff_min, uint32 rec_len, uint32* kxmer_counters, uint64 n_kxmer_counters, uint32 n_threads) :
		kmer_len(kmer_len),
		n_parts(n_parts),
		parts_left(n_parts),
		sub_array_descs(sub_array_descs),
		buffer(buffer),
		cutoff_min(cutoff_min),
		rec_len(rec_len),		
		kxmer_counters(kxmer_counters),
		n_kxmer_counters(n_kxmer_counters),
		n_threads(n_threads)
	{		
		//calculate cumulative sum
		cumsum.resize(n_kxmer_counters / COMPACT_CUMSUM_PART_SIZE + 1);
		vector<thread> threads;
		uint64 per_thread = (n_kxmer_counters / n_threads / COMPACT_CUMSUM_PART_SIZE + 1) * COMPACT_CUMSUM_PART_SIZE;
		for (uint32 th_id = 0; th_id < n_threads; ++th_id)
		{
			threads.emplace_back([&, th_id]()
			{
				uint64 start = th_id * per_thread;
				uint64 end = (th_id + 1) * per_thread;
				if (th_id == n_threads - 1)
					end = n_kxmer_counters;
				if (end > n_kxmer_counters)
					end = n_kxmer_counters;
				uint64 sum = 0;

				for (uint64 i = start; i < end; ++i)
				{
					sum += kxmer_counters[i];
					cumsum[i / COMPACT_CUMSUM_PART_SIZE] = sum;
				}
			});
		}
		for (auto& t : threads)
			t.join();
		uint64 part_size = per_thread / COMPACT_CUMSUM_PART_SIZE;
		for (uint32 i = 1; i < n_threads && i*part_size - 1 < cumsum.size(); ++i)
		{
			uint64 to_add = cumsum[i*part_size - 1];
			for (uint64 j = i * part_size; j < (i + 1)*part_size&& j < cumsum.size(); ++j)
				cumsum[j] += to_add;
		}

		//find biggest
		uint32 biggest_id = 0;
		for (uint32 i = 0; i < sub_array_descs.size(); ++i)
		{
			if (sub_array_descs[i].end - sub_array_descs[i].start > sub_array_descs[biggest_id].end - sub_array_descs[biggest_id].start)
				biggest_id = i;
		}

		uint64 start_in_biggest = 0;

		vector<SubArrayDesc> sub_array_desc_copy(sub_array_descs.begin(), sub_array_descs.end());

		CKmer<SIZE> mask;
		mask.set_n_1(kmer_len * 2);		
		while (parts_left > 1)
		{
			uint64 end_in_biggest = (sub_array_desc_copy[biggest_id].end - sub_array_desc_copy[biggest_id].start - start_in_biggest) / parts_left + sub_array_desc_copy[biggest_id].start;
			CKmer<SIZE> kmer;
			kmer.from_kxmer(buffer[end_in_biggest], sub_array_desc_copy[biggest_id].shr, mask);

			vector<SubArrayDesc> current;
			for (auto& e : sub_array_desc_copy)
			{
				uint32 shr = e.shr;
				uint64 new_end = std::lower_bound(buffer + e.start , buffer + e.end, kmer, [shr, mask](const CKmer<SIZE>& k1, const CKmer<SIZE>& k2)
				{
					CKmer<SIZE> val;
					val.from_kxmer(k1, shr, mask);
					return val < k2;
				}) - buffer;
				current.push_back({ e.start, new_end, shr, 0 });
				e.start = new_end;
			}

			//count exact number of k-mers in each subarray of part						
			for (auto& elem : current)
				elem.counters_sum = GetCumSum(elem.end) - GetCumSum(elem.start);			
			data.push(move(current));
			--parts_left;
		}

		//last
		for (auto& elem : sub_array_desc_copy)
			elem.counters_sum = GetCumSum(elem.end) - GetCumSum(elem.start);		
		data.push(move(sub_array_desc_copy));		
	}	


	uint64 GetCumSum(uint64 pos)
	{
		if (pos == 0)
			return 0;
		--pos;
		uint64 res = 0;
		if (pos / COMPACT_CUMSUM_PART_SIZE > 0)
			res = cumsum[pos / COMPACT_CUMSUM_PART_SIZE - 1];
		for (uint64 i = pos / COMPACT_CUMSUM_PART_SIZE * COMPACT_CUMSUM_PART_SIZE; i <= pos; ++i)
			res += kxmer_counters[i];
		return res;
	}

	bool GetNext(vector<SubArrayDesc>& desc, uint64& _out_start)
	{
		std::lock_guard<std::mutex> lck(mtx);
		if (data.empty())
			return false;
		desc = move(data.front());
		data.pop();		
		_out_start = out_start;
		
		uint64 n_recs = 0;
		for (auto& e : desc)
			n_recs += e.counters_sum;
		out_start += ((n_recs) / MAX(cutoff_min, 1u)) * rec_len;
			
		return true;
	}
};

class CLutUpdater
{
	uint64* lut;
	std::mutex mtx;
public:
	CLutUpdater(uint64* lut) :lut(lut)
	{
	}
	void UpdateLut(uint64 prefix, uint64 val)
	{
		std::lock_guard<std::mutex> lck(mtx);
		lut[prefix] += val;
	}
};

template<unsigned SIZE>
class CKXmerMerger
{
	const vector<SubArrayDesc>& sub_array_descs;
	CSubArrayDescGenerator<SIZE>& sub_array_desc_generator;
	CLutUpdater& lut_updater;
	uint64 n_total = 0;
	uint64 n_unique = 0;
	uint64 n_cutoff_min = 0;
	uint64 n_cutoff_max = 0;

	CKmer<SIZE>* buffer = nullptr;
	uint32* kxmer_counters = nullptr;
	uint32 cutoff_min;
	uint32 cutoff_max;
	uint32 counter_max;
	uint32 kmer_len;	
	CKXmerSet<SIZE> kxmer_set;
	uint64* lut;
	uint32 counter_size;
	int32 lut_prefix_len;	
	uchar* out_buffer;
	bool without_output;

	list<pair<uint64, uint64>> packs;
public:

	void GetStats(uint64& _n_unique, uint64& _n_cutoff_min, uint64& _n_cutoff_max, uint64& _n_total)
	{
		_n_unique = n_unique;
		_n_cutoff_min = n_cutoff_min;
		_n_cutoff_max = n_cutoff_max;
		_n_total = n_total;
	}

	void Process(vector<SubArrayDesc>& desc, /*uint32 prefix, */uint64 out_start)
	{
		kxmer_set.clear();
		kxmer_set.set_buffer(buffer);

		CKmer<SIZE> mask;
		mask.set_n_1(kmer_len * 2);
		
		uint64 last_prefix = 0;
		uint64 last_prefix_n_recs = 0;
		uint64 first_prefix = 1ull << 2 * lut_prefix_len;
		uint64 first_prefix_n_recs = 0;

		uint32 suffix_len_bits = (kmer_len - lut_prefix_len) * 2;
		uint64 kmer_bytes = suffix_len_bits / 8;

		CKmer<SIZE> candidate_min, candidate_max;
		for (auto &d : desc)
		{
			if (d.end > d.start)
			{
				kxmer_set.init_add(d.start, d.end, d.shr);				
				candidate_min.from_kxmer(buffer[d.start], d.shr, mask);
				candidate_max.from_kxmer(buffer[d.end - 1], d.shr, mask);
				
				uint64 candidate_min_prefix = candidate_min.remove_suffix(suffix_len_bits);
				uint64 candidate_max_prefix = candidate_max.remove_suffix(suffix_len_bits);

				if (candidate_max_prefix > last_prefix)
					last_prefix = candidate_max_prefix;

				if (candidate_min_prefix < first_prefix)
					first_prefix = candidate_min_prefix;
			}
		}
		uint64 counter_pos = 0;

		CKmer<SIZE> kmer, next_kmer;
		kmer.clear();
		next_kmer.clear();	
		uint32 count;


		uint64 out_pos = out_start; 

		//first
		if (kxmer_set.get_min(counter_pos, kmer))
		{			
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
					else if (count > cutoff_max)
						n_cutoff_max++;
					else
					{						
						if (count > counter_max)
							count = counter_max;

						if (!without_output)
						{
							uint64 prefix = kmer.remove_suffix(suffix_len_bits);
							if (prefix == last_prefix)
								++last_prefix_n_recs;
							else if (prefix == first_prefix)
								++first_prefix_n_recs;
							else
								++lut[prefix];

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
			else if (count > cutoff_max)
				n_cutoff_max++;
			else
			{
				if (count > counter_max)
					count = counter_max;

				if (!without_output)
				{
					uint64 prefix = kmer.remove_suffix(suffix_len_bits);
					if (prefix == last_prefix)
						++last_prefix_n_recs;
					else if (prefix == first_prefix)
						++first_prefix_n_recs;
					else
						++lut[prefix];

					for (int32 j = (int32)kmer_bytes - 1; j >= 0; --j)
						out_buffer[out_pos++] = kmer.get_byte(j);
					for (int32 j = 0; j < (int32)counter_size; ++j)
						out_buffer[out_pos++] = (count >> (j * 8)) & 0xFF;
				}
			}
			if (!without_output)
			{
				lut_updater.UpdateLut(last_prefix, last_prefix_n_recs);
				lut_updater.UpdateLut(first_prefix, first_prefix_n_recs);
				packs.emplace_back(out_start, out_pos);
			}
		}
	}

	CKXmerMerger(const vector<SubArrayDesc>& sub_array_descs, 
		CSubArrayDescGenerator<SIZE>& sub_array_desc_generator, 
		CLutUpdater& lut_updater, 
		CKmer<SIZE>* buffer,
		uint32* kxmer_counters, 
		uint32 cutoff_min, 
		uint32 cutoff_max, 
		uint32 counter_max, 
		uint32 kmer_len, 
		uint64* lut, 
		uint32 counter_size, 
		int32 lut_prefix_len, 
		uchar* out_buffer,
		bool without_output)
			:
		sub_array_descs(sub_array_descs), 
		sub_array_desc_generator(sub_array_desc_generator), 
		lut_updater(lut_updater),
		buffer(buffer), 
		kxmer_counters(kxmer_counters), 
		cutoff_min(cutoff_min),
		cutoff_max(cutoff_max),
		counter_max(counter_max), 
		kmer_len(kmer_len), 
		kxmer_set(kmer_len), 
		lut(lut),
		counter_size(counter_size), 
		lut_prefix_len(lut_prefix_len), 
		out_buffer(out_buffer),
		without_output(without_output)
	{
		
	}

	void operator()()
	{
		vector<SubArrayDesc> desc;		
		uint64 out_start;		
		while (sub_array_desc_generator.GetNext(desc, out_start))
		{
			Process(desc, out_start);	
		}
	}

	std::list<std::pair<uint64, uint64>>& GetPacks()
	{
		return packs;
	}
};

template<unsigned SIZE>
class CKXmerSetMultiThreaded
{	
	vector<SubArrayDesc> sub_array_descs;
	CKmer<SIZE>* buffer;
	uint32* kxmer_counters;
	uint64 n_kxmer_counters;
	uint32 cutoff_min;
	uint32 cutoff_max;
	uint32 counter_max;
	uint32 kmer_len;	
	int32 lut_prefix_len;
	uint64* lut;
	uchar* out_buffer;
	uint32 n_threads = 0;	

	uint64 n_unique = 0; 
	uint64 n_cutoff_min = 0;
	uint64 n_cutoff_max = 0;
	uint64 n_total = 0;

	list<pair<uint64, uint64>> output_packs_desc;

public:
	CKXmerSetMultiThreaded(CKmer<SIZE>* buffer,
		uint32* kxmer_counters, 
		uint64 n_kxmer_counters,
		uint32 cutoff_min, 
		uint32 cutoff_max, 
		uint32 counter_max, 
		uint32 kmer_len, 
		int32 lut_prefix_len, 
		uint64* lut, 
		uchar* out_buffer, 
		uint32 n_threads)
			:
		buffer(buffer), 
		kxmer_counters(kxmer_counters), 
		n_kxmer_counters(n_kxmer_counters),
		cutoff_min(cutoff_min), 
		cutoff_max(cutoff_max), 
		counter_max(counter_max),		
		kmer_len(kmer_len),
		lut_prefix_len(lut_prefix_len), 
		lut(lut), 
		out_buffer(out_buffer),
		n_threads(n_threads)		
	{		
	}
	void InitAdd(uint64 start, uint64 end, uint32 shr)
	{
		sub_array_descs.emplace_back(SubArrayDesc{ start, end, shr, 0 });
	}

	void Process(bool without_output)
	{
		uint32 n_parts = 8 * n_threads;
		
		CLutUpdater lut_updater(lut);
		vector<thread> threads;
		vector<CKXmerMerger<SIZE>*> mergers;
		uint32 counter_size = min(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));

		uint32 rec_len = (kmer_len - lut_prefix_len) / 4 + counter_size;

		CSubArrayDescGenerator<SIZE> sub_array_desc_generator(kmer_len, n_parts, sub_array_descs, buffer, cutoff_min, rec_len, kxmer_counters, n_kxmer_counters, n_threads);
		for (uint32 i = 0; i < n_threads; ++i)
		{
			mergers.push_back(new CKXmerMerger<SIZE>(sub_array_descs, sub_array_desc_generator, lut_updater, buffer, kxmer_counters, cutoff_min, 
				cutoff_max, counter_max, kmer_len, lut, counter_size, lut_prefix_len, out_buffer, without_output));
			threads.push_back(thread(ref(*mergers.back())));
		}

		for (auto& t : threads)
			t.join();

		uint64 tmp_n_unique = 0;
		uint64 tmp_n_cutoff_min = 0;
		uint64 tmp_n_cutoff_max = 0;
		uint64 tmp_n_total = 0;

		for (auto ptr : mergers)
		{
			ptr->GetStats(tmp_n_unique, tmp_n_cutoff_min, tmp_n_cutoff_max, tmp_n_total);
			output_packs_desc.splice(output_packs_desc.end(), move(ptr->GetPacks()));
			n_unique += tmp_n_unique;
			n_cutoff_min += tmp_n_cutoff_min;
			n_cutoff_max += tmp_n_cutoff_max;
			n_total += tmp_n_total;
		}

		for (auto ptr : mergers)
			delete ptr;
		output_packs_desc.sort([](const pair<uint64, uint64>& e1, const pair<uint64, uint64>& e2){return e1.first < e2.first; });
	}

	void GetStats(uint64& _n_unique, uint64& _n_cutoff_min, uint64& _n_cutoff_max, uint64& _n_total)
	{
		_n_unique = n_unique;
		_n_cutoff_min = n_cutoff_min;
		_n_cutoff_max = n_cutoff_max;
		_n_total = n_total;
	}	

	list<pair<uint64, uint64>>& GetOutputPacksDesc()
	{
		return output_packs_desc;
	}
};

#endif

// ***** EOF