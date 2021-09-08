/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _DEFS_H
#define _DEFS_H

#include <cinttypes>

#define KMC_VER		"3.1.1"
#define KMC_DATE	"2019-05-19"

#define _CRT_SECURE_NO_WARNINGS

#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define MAX(x,y)	((x) > (y) ? (x) : (y))
#define NORM(x, lower, upper)	((x) < (lower) ? (lower) : (x) > (upper) ? (upper) : (x))

#define uchar	unsigned char

#include <time.h>

//#define DEBUG_MODE

//#define USE_META_PROG

#define COMPACT_CUMSUM_PART_SIZE (1<<10)

#define KMER_X		3

#define STATS_FASTQ_SIZE (1 << 28)

#define EXPAND_BUFFER_RECS (1 << 16)

#define MIN_N_BINS 64
#define MAX_N_BINS 2000

#ifndef MAX_K
#define MAX_K		256
#endif

#define MIN_K		1

#define MIN_MEM		2

// Range of number of FASTQ/FASTA reading threads
#define MIN_SF		1
#define MAX_SF		32

// Range of number of signature length
#define MIN_SL		5
#define MAX_SL		11

// Range of number of splitting threads
#define MIN_SP		1
#define MAX_SP		64

// Range of number of threads for 2nd stage
#define MIN_SR		1
#define MAX_SR		128

//Range of number of sorter threads pre sorter in strict memory mode
#define MIN_SMSO	1
#define MAX_SMSO	16

//Range of number of uncompactor threads in strict memory mode
#define MIN_SMUN	1
#define MAX_SMUN	16

//Range of number of merger threads in strict memory mode
#define MIN_SMME	1
#define MAX_SMME	16



typedef float	count_t;

#define KMER_WORDS		((MAX_K + 31) / 32)


#ifdef WIN32
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
using int32 = int32_t;
using uint32 = uint32_t;
using int64 = int64_t;
using uint64 = uint64_t;

#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell

using int32 = int32_t;
using uint32 = uint32_t;
using int64 = int64_t;
using uint64 = uint64_t;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <stdio.h>
#include <ext/algorithm>
using __gnu_cxx::copy_n;

#endif


const int32 MAX_STR_LEN = 32768;
#define ALIGNMENT 0x100

#define BYTE_LOG(x) (((x) < (1 << 8)) ? 1 : ((x) < (1 << 16)) ? 2 : ((x) < (1 << 24)) ? 3 : 4)

#define BYTE_LOG_ULL(x) (((x) < (1ull << 8)) ? 1 : ((x) < (1ull << 16)) ? 2 : ((x) < (1ull << 24)) ? 3 : ((x) < (1ull << 32)) ? 4 : ((x) < (1ull << 40) ? 5 : ((x) < (1ull << 48) ? 6 : ((x) < (1ull << 56)) ? 7 : 8)))

#ifdef WIN32
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif

//for radix
#define WIN_ALIGNMENT 64
//designated experimentally
constexpr int32 BUFFER_WIDTHS[] = { -1, 32, 16, 16, 8, 8, 4, 8, 4 };
constexpr int32_t BUFFER_WIDTHS_ABOVE_CACHE_LINE_SIZE[] = { -1, 8, 4, 8, 2, 8, 4, 8, 1 };
constexpr uint32_t GetBufferWidth(uint32_t index)
{
	return (index <= 8) ? BUFFER_WIDTHS[index] : BUFFER_WIDTHS_ABOVE_CACHE_LINE_SIZE[(index - 1) % 8 + 1];
}

#ifdef WIN32
#define ALIGN_ARRAY __declspec(align(WIN_ALIGNMENT))
#else
#define ALIGN_ARRAY __attribute__((aligned(ALIGNMENT)))
#endif

#ifdef WIN32
typedef unsigned __int8 uint8_t;
#else
#include <stdint.h>
#endif


inline uint32 calc_counter_size(int64 cutoff_max, int64 counter_max)
{
	if (counter_max == 1)
		return 0;
	return MIN(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));
}

inline uint32 calc_counter_size_ull(int64 cutoff_max, int64 counter_max)
{
	if (counter_max == 1)
		return 0;
	return MIN(BYTE_LOG_ULL((uint64)cutoff_max), BYTE_LOG_ULL((uint64)counter_max));
}



#endif

// ***** EOF
