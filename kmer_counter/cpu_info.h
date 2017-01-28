/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#ifndef _CPU_INFO_
#define _CPU_INFO_
#include <string>

class CCpuInfo
{
	static void cpuid(int output[4], int functionnumber);
public:
	static std::string GetVendor();
	static std::string GetBrand();
};

#endif