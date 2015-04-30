/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.2.0
  Date   : 2015-04-15
*/
#ifndef _KXMER_SET_
#define _KXMER_SET_
#include "../kmc/definitions.h"
#include <tuple>

using namespace std;

#define KXMER_SET_SIZE 1024 



template <typename KMER_T, unsigned SIZE>
class CKXmerSet
{
	typedef tuple<uint64, uint64, uint32> elem_desc_t; //start_pos, end_pos, shr
	typedef pair<KMER_T, uint32> heap_elem_t; //kxmer val, desc_id
	elem_desc_t data_desc[KXMER_SET_SIZE];
	heap_elem_t data[KXMER_SET_SIZE];
	uint32 pos;
	uint32 desc_pos;
	KMER_T mask;

	KMER_T* buffer;

	inline void update_heap()
	{
		uint32 desc_id = data[1].second;
		KMER_T kmer;
		if (++get<0>(data_desc[desc_id]) < get<1>(data_desc[desc_id]))
		{
			kmer.from_kxmer(buffer[get<0>(data_desc[desc_id])], get<2>(data_desc[desc_id]), mask);
		}
		else
		{
			kmer.set(data[--pos].first);
			desc_id = data[pos].second;
		}
		
		uint32 parent, less;
		parent = less = 1;
		while (true)
		{
			if (parent * 2 >= pos)
				break;
			if (parent * 2 + 1 >= pos)
				less = parent * 2;
			else if (data[parent * 2].first < data[parent * 2 + 1].first)
				less = parent * 2;
			else
				less = parent * 2 + 1;
			if (data[less].first < kmer)
			{
				data[parent] = data[less];
				parent = less;
			}			
			else
				break;
		}
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
	inline void set_buffer(KMER_T* _buffer)
	{
		buffer = _buffer;
	}
	inline void clear()
	{
		pos = 1;
		desc_pos = 0;
	}

	inline bool get_min(uint64& _pos, KMER_T& kmer)
	{
		if (pos <= 1)
			return false;

		kmer = data[1].first;
		_pos = get<0>(data_desc[data[1].second]);
		update_heap();
		
	
		return true;
	}
};



#endif
