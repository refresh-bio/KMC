/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  This file demonstrates the example usage of kmc_api software. 
  It reads kmer_counter's output and prints kmers to an output file.

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.4
  Date   : 2024-02-09
*/

#include <string>
#include <cmath>
#include "../kmc_api/kmer_defs.h"

#ifndef _NC_UTILS_H
#define _NC_UTILS_H
class CNumericConversions {
public:
    static uchar digits[100000*5];
	static int powOf10[30];
    struct _si {
        _si()
        {
            for(int i = 0; i < 100000; ++i)
            {
                int dig = i;

                digits[i*5+4] = '0' + (dig % 10);
                dig /= 10;
                digits[i*5+3] = '0' + (dig % 10);
                dig /= 10;
                digits[i*5+2] = '0' + (dig % 10);
                dig /= 10;
                digits[i*5+1] = '0' + (dig % 10);
                dig /= 10;
                digits[i*5+0] = '0' + dig;
            }
			powOf10[0] = 1;
			for(int i = 1 ; i < 30 ; ++i)
			{
				powOf10[i] = powOf10[i-1]*10;
			}
        }
    } static _init;

    static int NDigits(uint64 val)
    {
        if(val >= 10000)
            return 5;
        else if(val >= 1000)
            return 4;
        else if(val >= 100)
            return 3;
        else if(val >= 10)
            return 2;
        else
            return 1;
    }

    static int Int2PChar(uint64 val, uchar *str)
    {
        if(val >= 1000000000000000ull)
        {
            uint64 dig1 = val / 1000000000000000ull;
            val -= dig1 * 1000000000000000ull;
            uint64 dig2 = val / 10000000000ull;
            val -= dig2 * 10000000000ull;
            uint64 dig3 = val / 100000ull;
            uint64 dig4 = val - dig3 * 100000ull;

            int ndig = NDigits(dig1);

            memcpy(str, digits+dig1*5+(5-ndig), ndig);
            memcpy(str+ndig, digits+dig2*5, 5);
            memcpy(str+ndig+5, digits+dig3*5, 5);
            memcpy(str+ndig+10, digits+dig4*5, 5);

            return ndig+15;
        }
        else if(val >= 10000000000ull)
        {
            uint64 dig1 = val / 10000000000ull;
            val -= dig1 * 10000000000ull;
            uint64 dig2 = val / 100000ull;
            uint64 dig3 = val - dig2 * 100000ull;

            int ndig = NDigits(dig1);

            memcpy(str, digits+dig1*5+(5-ndig), ndig);
            memcpy(str+ndig, digits+dig2*5, 5);
            memcpy(str+ndig+5, digits+dig3*5, 5);

            return ndig+10;
        }
        else if(val >= 100000ull)
        {
            uint64 dig1 = val / 100000ull;
            uint64 dig2 = val - dig1 * 100000ull;

            int ndig = NDigits(dig1);

            memcpy(str, digits+dig1*5+(5-ndig), ndig);
            memcpy(str+ndig, digits+dig2*5, 5);

            return ndig+5;
        }
        else
        {
            int ndig = NDigits(val);

            memcpy(str, digits+val*5+(5-ndig), ndig);

            return ndig;
        }
	}

	static int Double2PChar(double val, int prec, uchar *str)
	{
		double corrector = .5 / powOf10[prec];
		val += corrector;
		double ipart;
		double fractPart = std::modf(val, &ipart);
		uint32 intPart = (uint32)ipart;
		uint32 len = Int2PChar(intPart, str);
		uint32 pos = len;
		str[pos++] = '.';
		for(int i = 0 ; i < prec ; ++i)
		{
			fractPart *= 10;
			str[pos++] = '0' + (uint32)fractPart  % 10 ;
		}
		return len + prec + 1;
	}
};

#endif
