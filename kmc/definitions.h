/*
    This file is a part of KMC software distributed under GNU GPL 3 licence.
    The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

    Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

    Version: 220
    Date   : 20150415
*/

#ifndef _DEFS_H
#define _DEFS_H

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
typedef float			count_t;
typedef unsigned char 		uchar;
typedef int 			int32;
typedef unsigned int 		uint32;
typedef long long 		int64;
typedef unsigned long long 	uint64;

const int32 MAX_STR_LEN 	= 32768;

#define ALIGNMENT		0x100
#define MIN(x,y)		((x) < (y) ? (x) : (y))
#define MAX(x,y)		((x) > (y) ? (x) : (y))
#define NORM(x, lower, upper)	((x) < (lower) ? (lower) : (x) > (upper) ? (upper) : (x))
#define BYTE_LOG(x) 		(((x) < (1 << 8)) ? 1 : ((x) < (1 << 16)) ? 2 : ((x) < (1 << 24)) ? 3 : 4)
#define STATS_FASTQ_SIZE 	(1 << 28)
#define EXPAND_BUFFER_RECS 	(1 << 16)
#define KMER_WORDS		((MAX_K + 31) / 32)

#endif
// ***** EOF
