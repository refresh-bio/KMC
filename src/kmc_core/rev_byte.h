/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/

#ifndef _REV_BYTE_H
#define _REV_BYTE_H

#include "defs.h"
struct CRev_byte
{
	static uchar lut[256];
	struct _si
	{
		_si()
		{
			for (uint32 i = 0; i < 256; ++i)
				lut[i] = ((3 - (i & 3)) << 6) + ((3 - ((i >> 2) & 3)) << 4) + ((3 - ((i >> 4) & 3)) << 2) + (3 - ((i >> 6) & 3));
		}

	}static _init;
};

#endif

// ***** EOF