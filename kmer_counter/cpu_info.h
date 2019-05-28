/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _CPU_INFO_
#define _CPU_INFO_
#include <string>

class CCpuInfo
{

public:
	static const std::string& GetVendor();
	static const std::string& GetBrand();

	static bool SSE_Enabled();
	static bool SSE2_Enabled();
	static bool SSE3_Enabled();
	static bool SSE41_Enabled();
	static bool SSE42_Enabled();
	static bool AVX_Enabled();
	static bool AVX2_Enabled();

};

#endif

// ***** EOF