/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz and Agnieszka Debudaj-Grabysz

  Version: 3.1.1
  Date   : 2019-05-19
*/


#ifndef _KMER_DEFS_H
#define _KMER_DEFS_H

#define KMC_VER		"3.1.1"
#define KMC_DATE	"2019-05-19"

#ifndef MIN
#define MIN(x,y)	((x) < (y) ? (x) : (y))
#endif

#ifndef WIN32
	#include <stdint.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <cmath>
	#include <string.h>

	#define _TCHAR	char
	#define _tmain	main

	#define my_fopen    fopen
	#define my_fseek    fseek
	#define my_ftell    ftell


	#include <stdio.h>
	#include <iostream>

#else
	#define my_fopen    fopen
	#define my_fseek    _fseeki64
	#define my_ftell    _ftelli64
#endif
	//typedef unsigned char uchar;

	typedef int int32;
	typedef unsigned int uint32;
	typedef long long int64;
	typedef unsigned long long uint64;
	typedef unsigned char uchar;
#endif

// ***** EOF
