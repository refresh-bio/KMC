/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
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


//#define DISABLE_ASMLIB

//#define ENABLE_DEBUG
//#define ENABLE_LOGGER

#define KMC_VER		"2.3.0"
#define KMC_DATE	"2015-08-21"




#define DEFAULT_CIRCULAL_QUEUE_CAPACITY (4)

#define SUFIX_WRITE_QUEUE_CAPACITY (10)


#define KMC1_DB_READER_PREFIX_BUFF_BYTES (1 << 24)
#define KMC1_DB_READER_SUFIX_BUFF_BYTES  (1 << 24)

#define KMC2_DB_READER_PREFIX_BUFF_BYTES (1 << 24)
#define KMC2_DB_READER_SUFIX_BUFF_BYTES  (1 << 24)

#define KMC1_DB_WRITER_PREFIX_BUFF_BYTES (1 << 24)
#define KMC1_DB_WRITER_SUFIX_BUFF_BYTES  (1 << 24)

#define HISTOGRAM_MAX_COUNTER_DEFAULT 10000

#define DUMP_BUF_SIZE (1 << 24)

//Increasing this value will lead to more memory consumption, but from preliminary observations it has no performance(is sense of time) impact, so it is recommended to not change this value
#define BUNDLE_CAPACITY (1 << 12) //in kmers, for kmers and counters. 

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

#endif


// ***** EOF