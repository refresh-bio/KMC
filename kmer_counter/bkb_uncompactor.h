/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _BKB_UNCOMPACTOR_H
#define _BKB_UNCOMPACTOR_H

#include "params.h"
#include "kmer.h"
#include "rev_byte.h"


//************************************************************************************************************
// CBigKmerBinUncompactor - Unpacking super k-mers to k+x-mers, only in strict memory mode
//************************************************************************************************************
template<unsigned SIZE>
class CBigKmerBinUncompactor
{
	CBigBinPartQueue* bbpq;
	CBigBinKXmersQueue* bbkq;
	CMemoryPool *sm_pmm_expand;
	uint32 max_x;
	bool both_strands;
	uint32 kmer_len;

	CKmer<SIZE>* kxmers;
	int64 sm_mem_part_expand;
	uint32 kxmers_size;
	int32 bin_id;

	uchar* input_data;
	uint64 input_data_size;

	void GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p);
	void Uncompact();
	void ExpandKxmersBoth();
	void ExpandKxmersAll();
	void ExpandKmersBoth();
	void ExpandKmersAll();

	public:
	CBigKmerBinUncompactor(CKMCParams& Params, CKMCQueues& Queues);
	~CBigKmerBinUncompactor();
	void Uncompact(int32 _bin_id, uchar* _data, uint64 _size);
	
};

//************************************************************************************************************
// CBigKmerBinUncompactor
//************************************************************************************************************


//----------------------------------------------------------------------------------
template<unsigned SIZE> CBigKmerBinUncompactor<SIZE>::CBigKmerBinUncompactor(CKMCParams& Params, CKMCQueues& Queues)
{	
	sm_pmm_expand = Queues.sm_pmm_expand;
	bbpq = Queues.bbpq;
	bbkq = Queues.bbkq;
	kmer_len = Params.kmer_len;
	max_x = Params.max_x;
	both_strands = Params.both_strands;
	sm_mem_part_expand = Params.sm_mem_part_expand;
	kxmers_size = (uint32)(sm_mem_part_expand / sizeof(CKmer<SIZE>)); 
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor<SIZE>::Uncompact(int32 _bin_id, uchar* _data, uint64 _size)
{
	bin_id = _bin_id;
	input_data = _data;
	input_data_size = _size;
	Uncompact();
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> CBigKmerBinUncompactor<SIZE>::~CBigKmerBinUncompactor()
{

}


//----------------------------------------------------------------------------------
template <unsigned SIZE> inline void CBigKmerBinUncompactor<SIZE>::GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p)
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

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor<SIZE>::ExpandKxmersBoth()
{
	uchar* _raw_buffer;
	sm_pmm_expand->reserve(_raw_buffer);
	kxmers = (CKmer<SIZE>*)_raw_buffer;

	CKmer<SIZE> kmer, rev_kmer, kmer_mask;
	CKmer<SIZE> kxmer_mask;
	bool kmer_lower;
	uint32 x, additional_symbols;
	uchar symb;
	uint32 kmer_bytes = (kmer_len + 3) / 4;
	uint32 rev_shift = kmer_len * 2 - 2;
	uchar* data_p = input_data;
	kmer_mask.set_n_1(kmer_len * 2);
	uint32 kmer_shr = SIZE * 32 - kmer_len;

	kxmer_mask.set_n_1((kmer_len + max_x + 1) * 2);

	uint64 kxmers_pos = 0;
	uint64 pos = 0;
	while (pos < input_data_size)
	{
		kmer.clear();
		rev_kmer.clear();
		additional_symbols = data_p[pos++];

		//build kmer
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
			kxmers[kxmers_pos].set(kmer);
		else
			kxmers[kxmers_pos].set(rev_kmer);

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
					kxmers[kxmers_pos].SHL_insert_2bits(symb);
					++x;
					if (x == max_x)
					{
						if(!symbols_left)
							break;

						kxmers[kxmers_pos++].set_2bits(x, kmer_len * 2 + max_x * 2);
						if (kxmers_pos >= kxmers_size)
						{
							bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
							kxmers_pos = 0;
							sm_pmm_expand->reserve(_raw_buffer);
							kxmers = (CKmer<SIZE>*)_raw_buffer;
						}
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;
						if (kmer_lower)
							kxmers[kxmers_pos].set(kmer);
						else
							kxmers[kxmers_pos].set(rev_kmer);
					}
				}
				else
				{
					kxmers[kxmers_pos++].set_2bits(x, kmer_len * 2 + max_x * 2);
					if (kxmers_pos >= kxmers_size)
					{
						bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
						kxmers_pos = 0;
						sm_pmm_expand->reserve(_raw_buffer);
						kxmers = (CKmer<SIZE>*)_raw_buffer;
					}
					x = 0;

					kmer_lower = false;
					kxmers[kxmers_pos].set(rev_kmer);
				}
			}
			else
			{
				if (!(kmer < rev_kmer))
				{
					kxmers[kxmers_pos].set_2bits(3 - symb, kmer_len * 2 + x * 2);
					++x;
					if (x == max_x)
					{
						if(!symbols_left)
							break;

						kxmers[kxmers_pos++].set_2bits(x, kmer_len * 2 + max_x * 2);
						if (kxmers_pos >= kxmers_size)
						{
							bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
							kxmers_pos = 0;
							sm_pmm_expand->reserve(_raw_buffer);
							kxmers = (CKmer<SIZE>*)_raw_buffer;
						}
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;

						if (kmer_lower)
							kxmers[kxmers_pos].set(kmer);
						else
							kxmers[kxmers_pos].set(rev_kmer);
					}
				}
				else
				{
					kxmers[kxmers_pos++].set_2bits(x, kmer_len * 2 + max_x * 2);
					if (kxmers_pos >= kxmers_size)
					{
						bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
						kxmers_pos = 0;
						sm_pmm_expand->reserve(_raw_buffer);
						kxmers = (CKmer<SIZE>*)_raw_buffer;
					}
					x = 0;
					
					kxmers[kxmers_pos].set(kmer);
					kmer_lower = true;
				}
			}

		}
		kxmers[kxmers_pos++].set_2bits(x, kmer_len * 2 + max_x * 2);
		if (kxmers_pos >= kxmers_size)
		{
			bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
			kxmers_pos = 0;
			sm_pmm_expand->reserve(_raw_buffer);
			kxmers = (CKmer<SIZE>*)_raw_buffer;
		}
		if (byte_shift != 6)
			++pos;
	}

	if (kxmers_pos)
	{
		bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
	}
	else
	{
		sm_pmm_expand->free(_raw_buffer);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor<SIZE>::ExpandKxmersAll()
{
	uchar* _raw_buffer;
	sm_pmm_expand->reserve(_raw_buffer);
	kxmers = (CKmer<SIZE>*)_raw_buffer;

	uint64 pos = 0;
	CKmer<SIZE> kmer_mask, kxmer, kxmer_mask;
	kxmer_mask.set_n_1((kmer_len + max_x) * 2);
	uchar *data_p = input_data;

	kmer_mask.set_n_1(kmer_len * 2);
	uint64 kxmers_pos = 0;

	while (pos < input_data_size)
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

		kxmers[kxmers_pos++].set(kxmer);
		if (kxmers_pos >= kxmers_size)
		{
			bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
			kxmers_pos = 0;
			sm_pmm_expand->reserve(_raw_buffer);
			kxmers = (CKmer<SIZE>*)_raw_buffer;
		}
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

			kxmers[kxmers_pos++].set(kxmer);
			if (kxmers_pos >= kxmers_size)
			{
				bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
				kxmers_pos = 0;
				sm_pmm_expand->reserve(_raw_buffer);
				kxmers = (CKmer<SIZE>*)_raw_buffer;
			}
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

			kxmers[kxmers_pos++].set(kxmer);
			if (kxmers_pos >= kxmers_size)
			{
				
				bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
				kxmers_pos = 0;
				sm_pmm_expand->reserve(_raw_buffer);
				kxmers = (CKmer<SIZE>*)_raw_buffer;
			}
		}
		if (byte_shift != 6)
			++pos;
	}
	if (kxmers_pos)
	{
		bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
	}
	else
	{
		sm_pmm_expand->free(_raw_buffer);
	}

}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor<SIZE>::ExpandKmersBoth()
{
	uchar* _raw_buffer;
	sm_pmm_expand->reserve(_raw_buffer);
	kxmers = (CKmer<SIZE>*)_raw_buffer;

	CKmer<SIZE> kmer, rev_kmer, kmer_can, kmer_mask;

	uint32 kmer_bytes = (kmer_len + 3) / 4;
	uint32 kmer_len_shift = (kmer_len - 1) * 2;
	kmer_mask.set_n_1(kmer_len * 2);
	uchar *data_p = input_data;

	uint64 kxmers_pos = 0;
	uint64 pos = 0;


	while (pos < input_data_size)
	{
		kmer.clear();
		rev_kmer.clear();
		uint32 additional_symbols = data_p[pos++];
		uchar symb;

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

		uint32 kmer_shr = SIZE * 32 - kmer_len;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		rev_kmer.mask(kmer_mask);

		kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
		kxmers[kxmers_pos++].set(kmer_can);
		if (kxmers_pos >= kxmers_size)
		{
			bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
			kxmers_pos = 0;
			sm_pmm_expand->reserve(_raw_buffer);
			kxmers = (CKmer<SIZE>*)_raw_buffer;
		}

		for (uint32 i = 0; i < additional_symbols; ++i)
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
			kxmers[kxmers_pos++].set(kmer_can);
			if (kxmers_pos >= kxmers_size)
			{
				bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
				kxmers_pos = 0;
				sm_pmm_expand->reserve(_raw_buffer);
				kxmers = (CKmer<SIZE>*)_raw_buffer;
			}
		}
		if (byte_shift != 6)
			++pos;
	}
	if (kxmers_pos)
	{
		bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
	}
	else
	{
		sm_pmm_expand->free(_raw_buffer);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor<SIZE>::ExpandKmersAll()
{	
	uchar* _raw_buffer;
	sm_pmm_expand->reserve(_raw_buffer);
	kxmers = (CKmer<SIZE>*)_raw_buffer;

	uint64 kxmers_pos = 0;
	uint64 pos = 0;
	CKmer<SIZE> kmer;
	uint32 kmer_bytes = (kmer_len + 3) / 4;

	CKmer<SIZE> kmer_mask;
	kmer_mask.set_n_1(kmer_len * 2);
	uchar *data_p = input_data;

	while (pos < input_data_size)
	{
		kmer.clear();
		uint32 additional_symbols = data_p[pos++];		
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1; i < kmer_bytes; ++i, --kmer_pos)
		{
			kmer.set_byte(kmer_pos, data_p[pos + i]);
		}
		pos += kmer_bytes;
		uchar byte_shift = 6 - (kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		uint32 kmer_shr = SIZE * 32 - kmer_len;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		kxmers[kxmers_pos++].set(kmer);
		if (kxmers_pos >= kxmers_size)
		{
			bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
			kxmers_pos = 0;
			sm_pmm_expand->reserve(_raw_buffer);
			kxmers = (CKmer<SIZE>*)_raw_buffer;
		}
		for (uint32 i = 0; i < additional_symbols; ++i)
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
			kxmers[kxmers_pos++].set(kmer);
			if (kxmers_pos >= kxmers_size)
			{
				bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
				kxmers_pos = 0;
				sm_pmm_expand->reserve(_raw_buffer);
				kxmers = (CKmer<SIZE>*)_raw_buffer;
			}
		}
		if (byte_shift != 6)
			++pos;
	}

	if (kxmers_pos)
	{
		bbkq->push(bin_id, (uchar*)kxmers, kxmers_pos);
	}
	else
	{
		sm_pmm_expand->free(_raw_buffer);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor<SIZE>::Uncompact()
{
	if (max_x)
	{
		if (both_strands)
			ExpandKxmersBoth();
		else
			ExpandKxmersAll();
	}
	else
	{
		if (both_strands)
			ExpandKmersBoth();
		else
			ExpandKmersAll();
	}
}


//************************************************************************************************************
// CWBigKmerBinUncompactor - wrapper for multithreading purposes
//************************************************************************************************************
template<unsigned SIZE>
class CWBigKmerBinUncompactor
{
	CBigKmerBinUncompactor<SIZE>* bkb_uncompactor;
	CBigBinPartQueue* bbpq;
	CBigBinKXmersQueue* bbkq;
	CMemoryPool* sm_pmm_input_file;
public:
	CWBigKmerBinUncompactor(CKMCParams& Params, CKMCQueues& Queues);
	~CWBigKmerBinUncompactor();
	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template<unsigned SIZE>
CWBigKmerBinUncompactor<SIZE>::CWBigKmerBinUncompactor(CKMCParams& Params, CKMCQueues& Queues)
{
	bkb_uncompactor = new CBigKmerBinUncompactor<SIZE>(Params, Queues);
	bbpq = Queues.bbpq;
	bbkq = Queues.bbkq;
	sm_pmm_input_file = Queues.sm_pmm_input_file;
}

//----------------------------------------------------------------------------------
// Destructor
template<unsigned SIZE>
CWBigKmerBinUncompactor<SIZE>::~CWBigKmerBinUncompactor()
{
	delete bkb_uncompactor;
}

//----------------------------------------------------------------------------------
// Execution
template<unsigned SIZE>
void CWBigKmerBinUncompactor<SIZE>::operator()()
{
	int32 bin_id;
	uchar* data;
	uint64 size;
	while (bbpq->pop(bin_id, data, size))
	{
		bkb_uncompactor->Uncompact(bin_id, data, size);		
		sm_pmm_input_file->free(data);
	}
	bbkq->mark_completed();
}


#endif 

// ***** EOF