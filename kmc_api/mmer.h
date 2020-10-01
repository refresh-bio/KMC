/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _MMER_H
#define _MMER_H
#include "kmer_defs.h"

// *************************************************************************
// *************************************************************************


class CMmer
{
	uint32 str;
	uint32 mask;
	uint32 current_val;
	uint32* norm;
	uint32 len;
	static uint32 norm5[1 << 10];
	static uint32 norm6[1 << 12];
	static uint32 norm7[1 << 14];	
	static uint32 norm8[1 << 16];
	static uint32 norm9[1 << 18];
	static uint32 norm10[1 << 20];
	static uint32 norm11[1 << 22];

	static bool is_allowed(uint32 mmer, uint32 len)
	{
		if ((mmer & 0x3f) == 0x3f)            // TTT suffix
			return false;
		if ((mmer & 0x3f) == 0x3b)            // TGT suffix
			return false;
		if ((mmer & 0x3c) == 0x3c)            // TG* suffix !!!! consider issue #152
			return false;

		for (uint32 j = 0; j < len - 3; ++j)
		if ((mmer & 0xf) == 0)                // AA inside
			return false;
		else
			mmer >>= 2;

		if (mmer == 0)            // AAA prefix
			return false;
		if (mmer == 0x04)        // ACA prefix
			return false;
		if ((mmer & 0xf) == 0)    // *AA prefix
			return false;
	
		return true;
	}

	friend class CSignatureMapper;
	struct _si
	{			
		static uint32 get_rev(uint32 mmer, uint32 len)
		{
			uint32 rev = 0;
			uint32 shift = len*2 - 2;
			for(uint32 i = 0 ; i < len ; ++i)
			{
				rev += (3 - (mmer & 3)) << shift;
				mmer >>= 2;
				shift -= 2;
			}
			return rev;
		}

		

		static void init_norm(uint32* norm, uint32 len)
		{
			uint32 special = 1 << len * 2;
			for(uint32 i = 0 ; i < special ; ++i)
			{				
				uint32 rev = get_rev(i, len);
				uint32 str_val = is_allowed(i, len) ? i : special;
				uint32 rev_val = is_allowed(rev, len) ? rev : special;
				norm[i] = MIN(str_val, rev_val);				
			}
		}

		_si()
		{
			init_norm(norm5, 5);
			init_norm(norm6, 6);
			init_norm(norm7, 7);
			init_norm(norm8, 8);
			init_norm(norm9, 9);
			init_norm(norm10, 10);
			init_norm(norm11, 11);
		}

	}static _init;
public:
	CMmer(uint32 _len);
	inline void insert(uchar symb);
	inline uint32 get() const;
	inline bool operator==(const CMmer& x);
	inline bool operator<(const CMmer& x);
	inline void clear();
	inline bool operator<=(const CMmer& x);
	inline void set(const CMmer& x);
	inline void insert(const char* seq);
	
};



//--------------------------------------------------------------------------
inline void CMmer::insert(uchar symb)
{
	str <<= 2;
	str += symb;
	str &= mask;

	current_val = norm[str];
}

//--------------------------------------------------------------------------
inline uint32 CMmer::get() const
{
	return current_val;
}

//--------------------------------------------------------------------------
inline bool CMmer::operator==(const CMmer& x)
{
	return current_val == x.current_val;
}

//--------------------------------------------------------------------------
inline bool CMmer::operator<(const CMmer& x)
{
	return current_val < x.current_val;
}

//--------------------------------------------------------------------------
inline void CMmer::clear()
{
	str = 0;
}

//--------------------------------------------------------------------------
inline bool CMmer::operator<=(const CMmer& x)
{
	return current_val <= x.current_val;
}

//--------------------------------------------------------------------------
inline void CMmer::set(const CMmer& x)
{
	str = x.str;
	current_val = x.current_val;
}

//--------------------------------------------------------------------------
inline void CMmer::insert(const char* seq)
{
	switch (len)
	{
	case 5: 
		str = (seq[0] << 8) + (seq[1] << 6) + (seq[2] << 4) + (seq[3] << 2) + (seq[4]);
		break;
	case 6:
		str = (seq[0] << 10) + (seq[1] << 8) + (seq[2] << 6) + (seq[3] << 4) + (seq[4] << 2) + (seq[5]);
		break;
	case 7:
		str = (seq[0] << 12) + (seq[1] << 10) + (seq[2] << 8) + (seq[3] << 6) + (seq[4] << 4 ) + (seq[5] << 2) + (seq[6]);
		break;
	case 8:
		str = (seq[0] << 14) + (seq[1] << 12) + (seq[2] << 10) + (seq[3] << 8) + (seq[4] << 6) + (seq[5] << 4) + (seq[6] << 2) + (seq[7]);
		break;
	case 9:
		str = (seq[0] << 16) + (seq[1] << 14) + (seq[2] << 12) + (seq[3] << 10) + (seq[4] << 8) + (seq[5] << 6) + (seq[6] << 4) + (seq[7] << 2) + (seq[8]);
		break;
	case 10:
		str = (seq[0] << 18) + (seq[1] << 16) + (seq[2] << 14) + (seq[3] << 12) + (seq[4] << 10) + (seq[5] << 8) + (seq[6] << 6) + (seq[7] << 4) + (seq[8] << 2) + (seq[9]);
		break;
	case 11:
		str = (seq[0] << 20) + (seq[1] << 18) + (seq[2] << 16) + (seq[3] << 14) + (seq[4] << 12) + (seq[5] << 10) + (seq[6] << 8) + (seq[7] << 6) + (seq[8] << 4) + (seq[9] << 2) + (seq[10]);
		break;
	default:
		break;
	}

	current_val = norm[str];
}


#endif