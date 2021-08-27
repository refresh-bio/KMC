/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef RADULS_H
#define RADULS_H

#include <functional>
#include "kmer.h"

#define MAGIC_NUMBER 8

template<typename KMER_T>
using SortFunction = function<void(KMER_T*, KMER_T*, uint64, uint32, uint32, CMemoryPool*)>;

namespace RadulsSort
{
	template<typename KMER_T>
	void RadixSortMSD_SSE2(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);

	template<typename KMER_T>
	void RadixSortMSD_SSE41(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);

	template<typename KMER_T>
	void RadixSortMSD_AVX(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);

	template<typename KMER_T>
	void RadixSortMSD_AVX2(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);
}

#endif // RADULS_H

// ***** EOF