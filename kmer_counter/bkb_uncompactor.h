/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.1
  Date   : 2014-12-18
*/

#ifndef _BKB_UNCOMPACTOR_H
#define _BKB_UNCOMPACTOR_H

#include "params.h"
#include "kmer.h"
#include "rev_byte.h"

//************************************************************************************************************
template<typename KMER_T, unsigned SIZE> class CBigKmerBinUncompactor_Impl;

//************************************************************************************************************
// CBigKmerBinUncompactor - Unpacking super k-mers to k+x-mers, only in strict memory mode
//************************************************************************************************************
template<typename KMER_T, unsigned SIZE>
class CBigKmerBinUncompactor
{
	CBigBinPartQueue* bbpq;
	CBigBinKXmersQueue* bbkq;
	CMemoryPool *sm_pmm_expand;
	uint32 max_x;
	bool both_strands;
	uint32 kmer_len;

	KMER_T* kxmers;
	int64 sm_mem_part_expand;
	uint32 kxmers_size;
	int32 bin_id;

	uchar* input_data;
	uint64 input_data_size;

	friend class CBigKmerBinUncompactor_Impl<KMER_T, SIZE>;
	public:
	CBigKmerBinUncompactor(CKMCParams& Params, CKMCQueues& Queues);
	~CBigKmerBinUncompactor();
	void Uncompact(int32 _bin_id, uchar* _data, uint64 _size);
	
};

//************************************************************************************************************
// CBigKmerBinUncompactor_Impl - implementation of k-mer type- and size-dependent functions
//************************************************************************************************************
template<typename KMER_T, unsigned SIZE>
class CBigKmerBinUncompactor_Impl
{
public:
	static void Uncompact(CBigKmerBinUncompactor<KMER_T, SIZE>& ptr);
};

template<unsigned SIZE>
class CBigKmerBinUncompactor_Impl < CKmer<SIZE>, SIZE >
{
public:
	static void GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p);
	static void Uncompact(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr);
	static void ExpandKxmersBoth(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr);
	static void ExpandKxmersAll(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr);
	static void ExpandKmersBoth(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr);
	static void ExpandKmersAll(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr);
};


template<unsigned SIZE>
class CBigKmerBinUncompactor_Impl < CKmerQuake<SIZE>, SIZE >
{
public:
	static void Uncompact(CBigKmerBinUncompactor<CKmerQuake<SIZE>, SIZE>& ptr);
};

//************************************************************************************************************
// CBigKmerBinUncompactor
//************************************************************************************************************


//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE> CBigKmerBinUncompactor<KMER_T, SIZE>::CBigKmerBinUncompactor(CKMCParams& Params, CKMCQueues& Queues)
{	
	sm_pmm_expand = Queues.sm_pmm_expand;
	bbpq = Queues.bbpq;
	bbkq = Queues.bbkq;
	kmer_len = Params.kmer_len;
	max_x = Params.max_x;
	both_strands = Params.both_strands;
	sm_mem_part_expand = Params.sm_mem_part_expand;
	kxmers_size = (uint32)(sm_mem_part_expand / sizeof(KMER_T)); 
}

//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE> void CBigKmerBinUncompactor<KMER_T, SIZE>::Uncompact(int32 _bin_id, uchar* _data, uint64 _size)
{
	bin_id = _bin_id;
	input_data = _data;
	input_data_size = _size;
	CBigKmerBinUncompactor_Impl<KMER_T, SIZE>::Uncompact(*this);
}

//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE> CBigKmerBinUncompactor<KMER_T, SIZE>::~CBigKmerBinUncompactor()
{

}

//************************************************************************************************************
// CBigKmerBinUncompactor_Impl
//************************************************************************************************************

//----------------------------------------------------------------------------------
template <unsigned SIZE> inline void CBigKmerBinUncompactor_Impl<CKmer<SIZE>, SIZE>::GetNextSymb(uchar& symb, uchar& byte_shift, uint64& pos, uchar* data_p)
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
template<unsigned SIZE> void CBigKmerBinUncompactor_Impl<CKmer<SIZE>, SIZE>::ExpandKxmersBoth(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr)
{
	uchar* _raw_buffer;
	ptr.sm_pmm_expand->reserve(_raw_buffer);
	ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;

	CKmer<SIZE> kmer, rev_kmer, kmer_mask;
	CKmer<SIZE> kxmer_mask;
	bool kmer_lower;
	uint32 x, additional_symbols;
	uchar symb;
	uint32 kmer_bytes = (ptr.kmer_len + 3) / 4;
	uint32 rev_shift = ptr.kmer_len * 2 - 2;
	uchar* data_p = ptr.input_data;
	kmer_mask.set_n_1(ptr.kmer_len * 2);
	uint32 kmer_shr = SIZE * 32 - ptr.kmer_len;

	kxmer_mask.set_n_1((ptr.kmer_len + ptr.max_x + 1) * 2);

	uint64 kxmers_pos = 0;
	uint64 pos = 0;
	while (pos < ptr.input_data_size)
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
			ptr.kxmers[kxmers_pos].set(kmer);
		else
			ptr.kxmers[kxmers_pos].set(rev_kmer);

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
					ptr.kxmers[kxmers_pos].SHL_insert_2bits(symb);
					++x;
					if (x == ptr.max_x)
					{
						if(!symbols_left)
							break;

						ptr.kxmers[kxmers_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
						if (kxmers_pos >= ptr.kxmers_size)
						{
							ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
							kxmers_pos = 0;
							ptr.sm_pmm_expand->reserve(_raw_buffer);
							ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
						}
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;
						if (kmer_lower)
							ptr.kxmers[kxmers_pos].set(kmer);
						else
							ptr.kxmers[kxmers_pos].set(rev_kmer);
					}
				}
				else
				{
					ptr.kxmers[kxmers_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
					if (kxmers_pos >= ptr.kxmers_size)
					{
						ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
						kxmers_pos = 0;
						ptr.sm_pmm_expand->reserve(_raw_buffer);
						ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
					}
					x = 0;

					kmer_lower = false;
					ptr.kxmers[kxmers_pos].set(rev_kmer);
				}
			}
			else
			{
				if (!(kmer < rev_kmer))
				{
					ptr.kxmers[kxmers_pos].set_2bits(3 - symb, ptr.kmer_len * 2 + x * 2);
					++x;
					if (x == ptr.max_x)
					{
						if(!symbols_left)
							break;

						ptr.kxmers[kxmers_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
						if (kxmers_pos >= ptr.kxmers_size)
						{
							ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
							kxmers_pos = 0;
							ptr.sm_pmm_expand->reserve(_raw_buffer);
							ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
						}
						x = 0;

						GetNextSymb(symb, byte_shift, pos, data_p);
						kmer.SHL_insert_2bits(symb);
						kmer.mask(kmer_mask);
						rev_kmer.SHR_insert_2bits(3 - symb, rev_shift);
						--symbols_left;

						kmer_lower = kmer < rev_kmer;

						if (kmer_lower)
							ptr.kxmers[kxmers_pos].set(kmer);
						else
							ptr.kxmers[kxmers_pos].set(rev_kmer);
					}
				}
				else
				{
					ptr.kxmers[kxmers_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
					if (kxmers_pos >= ptr.kxmers_size)
					{
						ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
						kxmers_pos = 0;
						ptr.sm_pmm_expand->reserve(_raw_buffer);
						ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
					}
					x = 0;
					
					ptr.kxmers[kxmers_pos].set(kmer);
					kmer_lower = true;
				}
			}

		}
		ptr.kxmers[kxmers_pos++].set_2bits(x, ptr.kmer_len * 2 + ptr.max_x * 2);
		if (kxmers_pos >= ptr.kxmers_size)
		{
			ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
			kxmers_pos = 0;
			ptr.sm_pmm_expand->reserve(_raw_buffer);
			ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
		}
		if (byte_shift != 6)
			++pos;
	}

	if (kxmers_pos)
	{
		ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
	}
	else
	{
		ptr.sm_pmm_expand->free(_raw_buffer);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor_Impl<CKmer<SIZE>, SIZE>::ExpandKxmersAll(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr)
{
	uchar* _raw_buffer;
	ptr.sm_pmm_expand->reserve(_raw_buffer);
	ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;

	uint64 pos = 0;
	CKmer<SIZE> kmer_mask, kxmer, kxmer_mask;
	kxmer_mask.set_n_1((ptr.kmer_len + ptr.max_x) * 2);
	uchar *data_p = ptr.input_data;

	kmer_mask.set_n_1(ptr.kmer_len * 2);
	uint64 kxmers_pos = 0;

	while (pos < ptr.input_data_size)
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

		ptr.kxmers[kxmers_pos++].set(kxmer);
		if (kxmers_pos >= ptr.kxmers_size)
		{
			ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
			kxmers_pos = 0;
			ptr.sm_pmm_expand->reserve(_raw_buffer);
			ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
		}
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

			ptr.kxmers[kxmers_pos++].set(kxmer);
			if (kxmers_pos >= ptr.kxmers_size)
			{
				ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
				kxmers_pos = 0;
				ptr.sm_pmm_expand->reserve(_raw_buffer);
				ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
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

			kxmer.set_2bits(kxmer_rest, (ptr.kmer_len + ptr.max_x) * 2);

			ptr.kxmers[kxmers_pos++].set(kxmer);
			if (kxmers_pos >= ptr.kxmers_size)
			{
				
				ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
				kxmers_pos = 0;
				ptr.sm_pmm_expand->reserve(_raw_buffer);
				ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
			}
		}
		if (byte_shift != 6)
			++pos;
	}
	if (kxmers_pos)
	{
		ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
	}
	else
	{
		ptr.sm_pmm_expand->free(_raw_buffer);
	}

}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor_Impl<CKmer<SIZE>, SIZE>::ExpandKmersBoth(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr)
{
	uchar* _raw_buffer;
	ptr.sm_pmm_expand->reserve(_raw_buffer);
	ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;

	CKmer<SIZE> kmer, rev_kmer, kmer_can, kmer_mask;

	uint32 kmer_bytes = (ptr.kmer_len + 3) / 4;
	uint32 kmer_len_shift = (ptr.kmer_len - 1) * 2;
	kmer_mask.set_n_1(ptr.kmer_len * 2);
	uchar *data_p = ptr.input_data;

	uint64 kxmers_pos = 0;
	uint64 pos = 0;


	while (pos < ptr.input_data_size)
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
		uchar byte_shift = 6 - (ptr.kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		uint32 kmer_shr = SIZE * 32 - ptr.kmer_len;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		rev_kmer.mask(kmer_mask);

		kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
		ptr.kxmers[kxmers_pos++].set(kmer_can);
		if (kxmers_pos >= ptr.kxmers_size)
		{
			ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
			kxmers_pos = 0;
			ptr.sm_pmm_expand->reserve(_raw_buffer);
			ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
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
			ptr.kxmers[kxmers_pos++].set(kmer_can);
			if (kxmers_pos >= ptr.kxmers_size)
			{
				ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
				kxmers_pos = 0;
				ptr.sm_pmm_expand->reserve(_raw_buffer);
				ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
			}
		}
		if (byte_shift != 6)
			++pos;
	}
	if (kxmers_pos)
	{
		ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
	}
	else
	{
		ptr.sm_pmm_expand->free(_raw_buffer);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor_Impl<CKmer<SIZE>, SIZE>::ExpandKmersAll(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr)
{	
	uchar* _raw_buffer;
	ptr.sm_pmm_expand->reserve(_raw_buffer);
	ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;

	uint64 kxmers_pos = 0;
	uint64 pos = 0;
	CKmer<SIZE> kmer;
	uint32 kmer_bytes = (ptr.kmer_len + 3) / 4;

	CKmer<SIZE> kmer_mask;
	kmer_mask.set_n_1(ptr.kmer_len * 2);
	uchar *data_p = ptr.input_data;

	while (pos < ptr.input_data_size)
	{
		kmer.clear();
		uint32 additional_symbols = data_p[pos++];		
		for (uint32 i = 0, kmer_pos = 8 * SIZE - 1; i < kmer_bytes; ++i, --kmer_pos)
		{
			kmer.set_byte(kmer_pos, data_p[pos + i]);
		}
		pos += kmer_bytes;
		uchar byte_shift = 6 - (ptr.kmer_len % 4) * 2;
		if (byte_shift != 6)
			--pos;

		uint32 kmer_shr = SIZE * 32 - ptr.kmer_len;

		if (kmer_shr)
			kmer.SHR(kmer_shr);

		kmer.mask(kmer_mask);
		ptr.kxmers[kxmers_pos++].set(kmer);
		if (kxmers_pos >= ptr.kxmers_size)
		{
			ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
			kxmers_pos = 0;
			ptr.sm_pmm_expand->reserve(_raw_buffer);
			ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
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
			ptr.kxmers[kxmers_pos++].set(kmer);
			if (kxmers_pos >= ptr.kxmers_size)
			{
				ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
				kxmers_pos = 0;
				ptr.sm_pmm_expand->reserve(_raw_buffer);
				ptr.kxmers = (CKmer<SIZE>*)_raw_buffer;
			}
		}
		if (byte_shift != 6)
			++pos;
	}

	if (kxmers_pos)
	{
		ptr.bbkq->push(ptr.bin_id, (uchar*)ptr.kxmers, kxmers_pos);
	}
	else
	{
		ptr.sm_pmm_expand->free(_raw_buffer);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor_Impl<CKmer<SIZE>, SIZE>::Uncompact(CBigKmerBinUncompactor<CKmer<SIZE>, SIZE>& ptr)
{
	if (ptr.max_x)
	{
		if (ptr.both_strands)
			ExpandKxmersBoth(ptr);
		else
			ExpandKxmersAll(ptr);
	}
	else
	{
		if (ptr.both_strands)
			ExpandKmersBoth(ptr);
		else
			ExpandKmersAll(ptr);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CBigKmerBinUncompactor_Impl<CKmerQuake<SIZE>, SIZE>::Uncompact(CBigKmerBinUncompactor<CKmerQuake<SIZE>, SIZE>& ptr)
{
	//"Not supported in current release"
}


//************************************************************************************************************
// CWBigKmerBinUncompactor - wrapper for multithreading purposes
//************************************************************************************************************
template<typename KMER_T, unsigned SIZE>
class CWBigKmerBinUncompactor
{
	CBigKmerBinUncompactor<KMER_T, SIZE>* bkb_uncompactor;
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
template<typename KMER_T, unsigned SIZE>
CWBigKmerBinUncompactor<KMER_T, SIZE>::CWBigKmerBinUncompactor(CKMCParams& Params, CKMCQueues& Queues)
{
	bkb_uncompactor = new CBigKmerBinUncompactor<KMER_T, SIZE>(Params, Queues);
	bbpq = Queues.bbpq;
	bbkq = Queues.bbkq;
	sm_pmm_input_file = Queues.sm_pmm_input_file;
}

//----------------------------------------------------------------------------------
// Destructor
template<typename KMER_T, unsigned SIZE>
CWBigKmerBinUncompactor<KMER_T, SIZE>::~CWBigKmerBinUncompactor()
{
	delete bkb_uncompactor;
}

//----------------------------------------------------------------------------------
// Execution
template<typename KMER_T, unsigned SIZE>
void CWBigKmerBinUncompactor<KMER_T, SIZE>::operator()()
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