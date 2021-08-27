/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#include "stdafx.h"
#include "cpu_info.h"
#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <vector>
#include <array>
#include <cstring>
#include <bitset>
using std::array;
using std::vector;
using std::string;
using std::bitset;

// This code mostly bases on https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex?view=vs-2017

//TODO: probably I should also check if extensions are enabled in the OS
static struct CpuInfoImpl {

	bool sse = false;
	bool sse2 = false;
	bool sse3 = false;
	bool sse4_1 = false;
	bool sse4_2 = false;
	bool avx = false;
	bool avx2 = false;

	string vendor, brand;
	void cpuid(int *result, int function_id) const
	{
#ifdef _MSC_VER
		__cpuidex(result, function_id, 0);

		//it seems clang defined __GNUC__ so __clang__ is checked first, althought in fact it seems __GNUC__ code seems to work also for clang

#elif defined(__clang__)							//basing on https://clang.llvm.org/doxygen/cpuid_8h_source.html
		__asm("xchgq  %%rbx,%q1\n"
		"cpuid\n"
			"xchgq  %%rbx,%q1"
			:"=a"(result[0]), "=r" (result[1]), "=c"(result[2]), "=d"(result[3])
			: "0"(function_id), "c"(0));
#elif defined(__GNUC__)								//basing on https://github.com/gcc-mirror/gcc/blob/master/gcc/config/i386/cpuid.h#L187		
		__asm__("cpuid\n\t"
			: "=a" (result[0]), "=b" (result[1]), "=c" (result[2]), "=d" (result[3]) : "0" (function_id), "c"(0));
#endif  
	}

	CpuInfoImpl()
	{
		array<int, 4> cpui = { -1 };
		cpuid(cpui.data(), 0);
		int nIds_ = cpui[0];
		vector<array<int, 4>> data_;
		for (int i = 0; i <= nIds_; ++i)
		{
			cpuid(cpui.data(), i);
			data_.push_back(cpui);
		}
		char _vendor[0x20]{};

		memcpy(_vendor, &data_[0][1], sizeof(int));
		memcpy(_vendor + 4, &data_[0][3], sizeof(int));
		memcpy(_vendor + 8, &data_[0][2], sizeof(int));

		vendor = _vendor;
		if (nIds_ > 0)
		{
			std::bitset<32> ECX = data_[1][2];
			std::bitset<32> EDX = data_[1][3];
			sse = EDX[25];
			sse2 = EDX[26];
			sse3 = ECX[0];
			sse4_1 = ECX[19];
			sse4_2 = ECX[20];
			avx = ECX[28];
		}

		if (nIds_ > 6)
		{
			std::bitset<32> EBX = data_[7][1];
			avx2 = EBX[5];
		}
	}

	const string& GetVendor() const
	{
		return vendor;
	}

	const string& GetBrand()
	{
		static bool computed = false;
		if (computed)
			return brand;
		array<int, 4> cpui = { -1 };
		cpuid(cpui.data(), 0x80000000);
		int nExIds_ = cpui[0];
		char _brand[0x40];
		memset(_brand, 0, sizeof(_brand));
		vector<array<int, 4>> extdata_;
		for (int i = 0x80000000; i <= nExIds_; ++i)
		{
			cpuid(cpui.data(), i);
			extdata_.push_back(cpui);
		}

		if ((unsigned int)nExIds_ >= 0x80000004)
		{
			memcpy(_brand, extdata_[2].data(), sizeof(cpui)); memcpy(_brand, extdata_[2].data(), sizeof(cpui));
			memcpy(_brand + 16, extdata_[3].data(), sizeof(cpui)); memcpy(_brand + 16, extdata_[3].data(), sizeof(cpui));
			memcpy(_brand + 32, extdata_[4].data(), sizeof(cpui)); memcpy(_brand + 32, extdata_[4].data(), sizeof(cpui));
			brand = _brand;
		}
		computed = true;
		return brand;
	}
} cpu_info_impl;

const string& CCpuInfo::GetVendor()
{
	return cpu_info_impl.GetVendor();
}

const string& CCpuInfo::GetBrand()
{
	return cpu_info_impl.GetBrand();
}


bool CCpuInfo::SSE_Enabled() { return cpu_info_impl.sse; }
bool CCpuInfo::SSE2_Enabled() { return cpu_info_impl.sse2; }
bool CCpuInfo::SSE3_Enabled() { return cpu_info_impl.sse3; }
bool CCpuInfo::SSE41_Enabled() { return cpu_info_impl.sse4_1; }
bool CCpuInfo::SSE42_Enabled() { return cpu_info_impl.sse4_2; }
bool CCpuInfo::AVX_Enabled() { return cpu_info_impl.avx; }
bool CCpuInfo::AVX2_Enabled() { return cpu_info_impl.avx2; }

// ***** EOF