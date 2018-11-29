/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _DEFS_H
#define _DEFS_H

#include <stdio.h>
#include <stdlib.h>

using uint32 = unsigned int;
using uint64 = unsigned long long;
using int32 = int;
using int64 = long long;
using uchar = unsigned char;

#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define MAX(x,y)	((x) > (y) ? (x) : (y))
#define NORM(x, lower, upper)	((x) < (lower) ? (lower) : (x) > (upper) ? (upper) : (x))

#define BYTE_LOG(x) (((x) < (1 << 8)) ? 1 : ((x) < (1 << 16)) ? 2 : ((x) < (1 << 24)) ? 3 : 4)


//#define ENABLE_DEBUG
//#define ENABLE_LOGGER

#define KMC_VER		"3.1.0"
#define KMC_DATE	"2018-05-10"

#define DEFAULT_CIRCULAL_QUEUE_CAPACITY (4)

#define SUFFIX_WRITE_QUEUE_CAPACITY (10)


#define KMC1_DB_READER_PREFIX_BUFF_BYTES (1 << 24)
#define KMC1_DB_READER_SUFFIX_BUFF_BYTES  (1 << 24)

#define KMC2_DB_READER_PREFIX_BUFF_BYTES (1 << 24)
#define KMC2_DB_READER_SUFFIX_BUFF_BYTES  (1 << 24)

#define KMC1_DB_WRITER_PREFIX_BUFF_BYTES (1 << 24)
#define KMC1_DB_WRITER_SUFFIX_BUFF_BYTES  (1 << 24)

#define HISTOGRAM_MAX_COUNTER_DEFAULT 10000

#define DUMP_BUF_SIZE (1 << 24)

//Increasing this value will lead to more memory consumption, but from preliminary observations it has no performance(is sense of time) impact, so it is recommended to not change this value
#define BUNDLE_CAPACITY (1 << 12) //in kmers, for kmers and counters. 

#define KMC2_DB_READER_BUNDLE_CAPACITY (1 << 22)

//this value has high impact to used memory, max value of memory is = 2 * SINGLE_BIN_BUFF_SIZE_FOR_DB2_READER * number_of_kmc2_input_dbs * number_of_bins_per_in_db
//increasing this value can have positive performance impact when running on HDD
#define SINGLE_BIN_BUFF_SIZE_FOR_DB2_READER (1 << 21) //if less is needed less will be allocated

//default values
#define CUTOFF_MIN 2 
#define CUTOFF_MAX 1000000000
#define COUNTER_MAX 255

#define MAX_K		256
#define KMER_WORDS		((MAX_K + 31) / 32)


#define USE_META_PROG

#ifdef WIN32
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell
#endif


#ifdef _MSC_VER
	#define _bswap_uint64(X) _byteswap_uint64(X)
	#define _bswap_uint32(X) _byteswap_ulong(X)
#elif defined(__GNUC__)
	#define _bswap_uint64(X) __builtin_bswap64(X)
	#define _bswap_uint32(X) __builtin_bswap32(X)
#else //unknown. Use the fastest "standard" way I've found
	#define _bswap_uint64(val) \
		val = ((val << 8) & 0xFF00FF00FF00FF00ULL) + ((val >> 8) & 0x00FF00FF00FF00FFULL); \
		val = ((val << 16) & 0xFFFF0000FFFF0000ULL) + ((val >> 16) & 0x0000FFFF0000FFFFULL); \
		val = (val << 32) + (val >> 32);
	#define _bswap_uint32(val) \
		val = (val<<24) | ((val<<8) & 0x00ff0000) | ((val >> 8) & 0x0000ff00) | (val >> 24);
#endif

#endif


// ***** EOF