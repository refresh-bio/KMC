/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2017-01-28
*/

#ifndef _SMALL_K_BUF
#define _SMALL_K_BUF

#include "defs.h"

template<typename COUNTER_TYPE>
struct CSmallKBuf
{
	COUNTER_TYPE* buf;
	void Store(uint64 index, uchar* _buf, uint32& buf_pos, uint64 counter_size)
	{
		for (uint64 j = 0; j < counter_size; ++j)
			_buf[buf_pos++] = (buf[index] >> (j * 8)) & 0xFF;
	}
};

#endif

// ***** EOF