/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include "cpu_info.h"
#include <vector>
#include <array>
#include <cstring>
using namespace std;

//Implementation taken from Agner Fog vectorclass library (http://www.agner.org/optimize/#vectorclass)
void CCpuInfo::cpuid(int output[4], int functionnumber)
{
#if defined (_MSC_VER) || defined (__INTEL_COMPILER)       // Microsoft or Intel compiler, intrin.h included

	__cpuidex(output, functionnumber, 0);                  // intrinsic function for CPUID

#elif defined(__GNUC__) || defined(__clang__)              // use inline assembly, Gnu/AT&T syntax

	int a, b, c, d;
	__asm("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(functionnumber), "c"(0) : );
	output[0] = a;
	output[1] = b;
	output[2] = c;
	output[3] = d;

#else                                                      // unknown platform. try inline assembly with masm/intel syntax

	__asm {
		mov eax, functionnumber
			xor ecx, ecx
			cpuid;
		mov esi, output
			mov[esi], eax
			mov[esi + 4], ebx
			mov[esi + 8], ecx
			mov[esi + 12], edx
	}

#endif
}

std::string CCpuInfo::GetVendor()
{
	array<int, 4> cpui = { -1 };
	cpuid(cpui.data(), 0);
	int nIds_ = cpui[0];
	std::vector<std::array<int, 4>> data_;
	for (int i = 0; i <= nIds_; ++i)
	{
		cpuid(cpui.data(), i);
		data_.push_back(cpui);
	}
	char vendor[0x20];
	memset(vendor, 0, sizeof(vendor));
	memcpy(vendor, &data_[0][1], sizeof(int));
	memcpy(vendor + 4, &data_[0][3], sizeof(int));
	memcpy(vendor + 8, &data_[0][2], sizeof(int));

	//*reinterpret_cast<int*>(vendor) = data_[0][1];
	//*reinterpret_cast<int*>(vendor + 4) = data_[0][3];
	//*reinterpret_cast<int*>(vendor + 8) = data_[0][2];
	return vendor;
}

std::string CCpuInfo::GetBrand()
{
	array<int, 4> cpui = { -1 };
	cpuid(cpui.data(), 0x80000000);
	int nExIds_ = cpui[0];
	char brand[0x40];
	memset(brand, 0, sizeof(brand));
	std::vector<std::array<int, 4>> extdata_;
	for (int i = 0x80000000; i <= nExIds_; ++i)
	{
		cpuid(cpui.data(), i);
		extdata_.push_back(cpui);
	}
	string brand_ = "";
	if ((unsigned int)nExIds_ >= 0x80000004)
	{
		memcpy(brand, extdata_[2].data(), sizeof(cpui)); memcpy(brand, extdata_[2].data(), sizeof(cpui));
		memcpy(brand + 16, extdata_[3].data(), sizeof(cpui)); memcpy(brand + 16, extdata_[3].data(), sizeof(cpui));
		memcpy(brand + 32, extdata_[4].data(), sizeof(cpui)); memcpy(brand + 32, extdata_[4].data(), sizeof(cpui));
		brand_ = brand;
	}
	return brand_;
}