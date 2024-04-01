/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
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
#ifndef __aarch64__
#ifdef HAVE_SSE2_INSTRUCTIONS
	template<typename KMER_T>
	void RadixSortMSD_SSE2(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);
#endif
#ifdef HAVE_SSE4_1_INSTRUCTIONS
	template<typename KMER_T>
	void RadixSortMSD_SSE41(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);
#endif
#ifdef HAVE_AVX_INSTRUCTIONS
	template<typename KMER_T>
	void RadixSortMSD_AVX(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);
#endif
#ifdef HAVE_AVX2_INSTRUCTIONS
	template<typename KMER_T>
	void RadixSortMSD_AVX2(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);
#endif
#else
	template<typename KMER_T>
	void RadixSortMSD_NEON(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf);
#endif
}

#endif // RADULS_H

// ***** EOF
