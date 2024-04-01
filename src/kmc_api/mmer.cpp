/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/

#include "../kmc_api/mmer.h"


uint32_t CMmer::norm5[];
uint32_t CMmer::norm6[];
uint32_t CMmer::norm7[];
uint32_t CMmer::norm8[];
uint32_t CMmer::norm9[];
uint32_t CMmer::norm10[];
uint32_t CMmer::norm11[];

CMmer::_si CMmer::_init;


//--------------------------------------------------------------------------
CMmer::CMmer(uint32_t _len)
{
	switch (_len)
	{
	case 5:
		norm = norm5;
		break;
	case 6:
		norm = norm6;
		break;
	case 7:
		norm = norm7;
		break;
	case 8:
		norm = norm8;
		break;
	case 9:
		norm = norm9;
		break;
	case 10:
		norm = norm10;
		break;
	case 11:
		norm = norm11;
		break;
	default:
		break;
	}
	len = _len;
	mask = (1 << _len * 2) - 1;
	str = 0;
}

//--------------------------------------------------------------------------

