/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.1
Date   : 2019-05-19
*/
#ifndef _SMALL_SORT_H
#define _SMALL_SORT_H

#include <cstdint>
#include <chrono>
#include <stdlib.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <functional>
#include <array>
#include <string>
#include <array>
#include "defs.h"
#include "kmer.h"

using namespace std;


template<unsigned SIZE>
class CSmallSort {
	static uint32 ArraySize;
	
	static vector<function<void(CKmer<SIZE> *, uint32)>> sorters;
	static vector<function<void(CKmer<SIZE> *, uint32)>> algorithms;

	static CKmer<SIZE> *arr, *arr_orig;
	static vector<vector<double>> sorter_times;

	static void std_sort(CKmer<SIZE> *ptr, uint32 size);
	static void ins_sort_loop(CKmer<SIZE> *ptr, uint32 size);
	static void ins_sort_hybrid(CKmer<SIZE> *ptr, uint32 size);
	static void ins_sort_macro(CKmer<SIZE> *ptr, uint32 size);
	static void shell_sort_1_7(CKmer<SIZE> *ptr, uint32 size);
	static void shell_sort_1_8(CKmer<SIZE> *ptr, uint32 size);
	static void shell_sort_1_10(CKmer<SIZE> *ptr, uint32 size);

	static void PrepareArray(void)
	{
		arr = new CKmer<SIZE>[ArraySize];
		arr_orig = new CKmer<SIZE>[ArraySize];

		mt19937_64 mt;

		sorter_times.clear();

		for (uint32 i = 0; i < ArraySize; ++i)
			for (uint32 j = 0; j < SIZE; ++j)
				arr_orig[i].random_init(j, mt());

	}

	static void ReleaseArray(void)
	{
		delete[] arr;
		delete[] arr_orig;
	}

	static void EvaluateAlgorithms(uint32 max_small_size)
	{
		algorithms.clear();

		algorithms.push_back(CSmallSort::std_sort);
		algorithms.push_back(CSmallSort::ins_sort_loop);
		algorithms.push_back(CSmallSort::ins_sort_hybrid);
		algorithms.push_back(CSmallSort::shell_sort_1_7);
		algorithms.push_back(CSmallSort::shell_sort_1_8);
		algorithms.push_back(CSmallSort::shell_sort_1_10);

		sorter_times.resize(max_small_size + 1);

		sorter_times.front().resize(algorithms.size());				// empty values for 0 elements

		for (uint32 part_size = 1; part_size <= max_small_size; ++part_size)
			for (uint32 j = 0; j < algorithms.size(); ++j)
			{
				for (int64_t i = 0; i < ArraySize; ++i)
					arr[i] = arr_orig[i];

				// ****************
				auto start = std::chrono::high_resolution_clock::now();

				for (uint64_t start = 0; start + part_size < ArraySize; start += part_size)
					algorithms[j](arr + start, part_size);

				auto end = std::chrono::high_resolution_clock::now();

				std::chrono::duration<double> diff = end - start;

				diff /= (ArraySize / part_size);

				sorter_times[part_size].push_back(diff.count());
			}
	}

	static void SmoothTimes(void)
	{
		auto n_sorters = sorter_times.front().size();
		auto n_values = sorter_times.size();

		for (uint32 i = 0; i < n_sorters; ++i)
		{
			for (uint32 j = n_values - 1; j > 0; --j)
				if (sorter_times[j - 1][i] > sorter_times[j][i])
					sorter_times[j - 1][i] = sorter_times[j][i];
		}
	}
	
	static void SelectBestSorters(void)
	{
		auto n_sorters = sorter_times.front().size();
		auto n_values = sorter_times.size();

		sorters.clear();

		for (uint32 i = 0; i < n_values; ++i)
		{
			uint32 best_id = 0;
			double best_time = sorter_times[i].front();

//			cout << i << " : " << best_time << " ";

			for (uint32 j = 1; j < n_sorters; ++j)
			{
				if (sorter_times[i][j] < best_time)
				{
					best_time = sorter_times[i][j];
					best_id = j;
				}
//				cout << sorter_times[i][j] << " ";
			}

//			cout << endl;
			
			sorters.push_back(algorithms[best_id]);

//			cout << best_id << endl;
		}
	}

public:
	CSmallSort();
	~CSmallSort();

	static void Adjust(uint32 max_small_size = 256)
	{
		PrepareArray();
		EvaluateAlgorithms(max_small_size);
		SmoothTimes();
		SelectBestSorters();
		ReleaseArray();
	}

	static void Sort(CKmer<SIZE> *ptr, uint32 size)
	{
		sorters[size](ptr, size);
	}
};


template<unsigned SIZE>
uint32 CSmallSort<SIZE>::ArraySize = 20 << 10;

template<unsigned SIZE>
vector<function<void(CKmer<SIZE> *, uint32)>> CSmallSort<SIZE>::sorters;

template<unsigned SIZE>
vector<function<void(CKmer<SIZE> *, uint32)>> CSmallSort<SIZE>::algorithms;

template<unsigned SIZE>
CKmer<SIZE> *CSmallSort<SIZE>::arr;

template<unsigned SIZE>
CKmer<SIZE> *CSmallSort<SIZE>::arr_orig;

template<unsigned SIZE>
vector<vector<double>> CSmallSort<SIZE>::sorter_times;

// ********** Sorters
template<unsigned SIZE>
void CSmallSort<SIZE>::std_sort(CKmer<SIZE> *ptr, uint32 size)
{
	std::sort(ptr, ptr + size);
}

template<unsigned SIZE>
void CSmallSort<SIZE>::ins_sort_loop(CKmer<SIZE> *ptr, uint32 size)
{
	int64_t i, j;
	CKmer<SIZE> x;
	int64_t n = size;

	CKmer<SIZE> *A = ptr;

	for (i = 1; i < n; i++)
	{
		x = A[i];
		j = i - 1;
		while ((j >= 0) && (x < A[j]))
		{
			A[j + 1] = A[j];
			j--;
		}
		A[j + 1] = x;
	}
}



#define INS_SORT_INNER_LOOP_BEGIN	\
	pA = ++qA;						\
	x = *pA;						\
	if (x < *(pA - 1))				\
		{								\
		*pA = *(pA - 1);			\
		--pA;

#define INS_SORT_INNER_LOOP_END		\
	*pA = x;						\
		}

#define INS_SORT_INNER_LOOP_IF_BEGIN	\
	if (x < *(pA - 1))					\
		{									\
		*pA = *(pA - 1);				\
		--pA;						

#define INS_SORT_INNER_LOOP_IF_END		\
		}

#define INS_SORT_INNER_LOOP_1		\
	if (A[1] < A[0])				\
		{								\
		CKmer<SIZE>  tmp = A[0];				\
		A[0] = A[1];				\
		A[1] = tmp;					\
		}


#define INS_SORT_INNER_LOOP_IF_1	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_2	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_1		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_3	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_2		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_4	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_3		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_5	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_4		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_6	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_5		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_7	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_6		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_8	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_7		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_9	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_8		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_10	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_9		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_11	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_10		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_12	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_11		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_13	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_12		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_14	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_13		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_15	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_14		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_16	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_15		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_17	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_16		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_18	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_17		INS_SORT_INNER_LOOP_IF_END
#define INS_SORT_INNER_LOOP_IF_19	\
	INS_SORT_INNER_LOOP_IF_BEGIN		INS_SORT_INNER_LOOP_IF_18		INS_SORT_INNER_LOOP_IF_END


#define INS_SORT_BODY_2					\
	CKmer<SIZE> x;								\
	CKmer<SIZE> *pA;								\
	CKmer<SIZE> *qA = A + 1;						\
	INS_SORT_INNER_LOOP_1

#define INS_SORT_BODY_3					\
	INS_SORT_BODY_2						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_1		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_4					\
	INS_SORT_BODY_3						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_2		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_5					\
	INS_SORT_BODY_4						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_3		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_6					\
	INS_SORT_BODY_5						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_4		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_7					\
	INS_SORT_BODY_6						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_5		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_8					\
	INS_SORT_BODY_7						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_6		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_9					\
	INS_SORT_BODY_8						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_7		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_10				\
	INS_SORT_BODY_9						INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_8		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_11				\
	INS_SORT_BODY_10					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_9		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_12				\
	INS_SORT_BODY_11					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_10		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_13				\
	INS_SORT_BODY_12					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_11		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_14				\
	INS_SORT_BODY_13					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_12		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_15				\
	INS_SORT_BODY_14					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_13		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_16				\
	INS_SORT_BODY_15					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_14		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_17				\
	INS_SORT_BODY_16					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_15		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_18				\
	INS_SORT_BODY_17					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_16		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_19				\
	INS_SORT_BODY_18					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_17		INS_SORT_INNER_LOOP_END
#define INS_SORT_BODY_20				\
	INS_SORT_BODY_19					INS_SORT_INNER_LOOP_BEGIN		INS_SORT_INNER_LOOP_IF_18		INS_SORT_INNER_LOOP_END



//---------------------------------------------------------------------------
// Sortowanie przez proste wstawianie bez wartownika
template<unsigned SIZE>
void CSmallSort<SIZE>::ins_sort_macro(CKmer<SIZE> *A, uint32 size)
{
	switch (size)
	{
	case 0:
	case 1:
		break;
	case 2:
		if (A[1] < A[0])				
		{								
			CKmer<SIZE> tmp = A[0];	
			A[0] = A[1];				
			A[1] = tmp;					
		}
		break;
	case 3:	 {INS_SORT_BODY_3};		break;
	case 4:  {INS_SORT_BODY_4};		break;
	case 5:  {INS_SORT_BODY_5};		break;
	case 6:  {INS_SORT_BODY_6};		break;
	case 7:  {INS_SORT_BODY_7};		break;
	case 8:  {INS_SORT_BODY_8};		break;
		/*	case 9:  {INS_SORT_BODY_9};		break;
		case 10: {INS_SORT_BODY_10};	break;
		case 11: {INS_SORT_BODY_11};	break;
		case 12: {INS_SORT_BODY_12};	break;
		case 13: {INS_SORT_BODY_13};	break;
		case 14: {INS_SORT_BODY_14};	break;
		case 15: {INS_SORT_BODY_15};	break;
		case 16: {INS_SORT_BODY_16};	break;
		case 17: {INS_SORT_BODY_17};	break;
		case 18: {INS_SORT_BODY_18};	break;
		case 19: {INS_SORT_BODY_19};	break;
		case 20: {INS_SORT_BODY_20};	break;*/
	}

	/*	int n = end - A;
	if (n < 16)
	if (n < 8)
	if (n < 4)
	if (n < 2)
	;
	else
	if (n == 2)		{INS_SORT_BODY_2}
	else			{ INS_SORT_BODY_3 }
	else
	if(n < 6)
	if(n == 4)		{ INS_SORT_BODY_4 }
	else { INS_SORT_BODY_5 }
	else
	if(n == 6)		{ INS_SORT_BODY_6}
	else { INS_SORT_BODY_7 }
	else
	if (n < 12)
	if(n < 10)
	if(n == 8)		{ INS_SORT_BODY_8 }
	else			{ INS_SORT_BODY_9}
	else
	if (n == 10) { INS_SORT_BODY_10 }
	else { INS_SORT_BODY_11 }
	else
	if(n < 14)
	if (n == 12) { INS_SORT_BODY_12 }
	else			{ INS_SORT_BODY_13 }
	else
	if (n == 14) { INS_SORT_BODY_14 }
	else { INS_SORT_BODY_15 }
	else
	if(n < 24)
	if(n < 20)
	if(n < 18)
	if (n == 16) { INS_SORT_BODY_16 }
	else { INS_SORT_BODY_17 }
	else
	if(n == 18) { INS_SORT_BODY_18 }
	else { INS_SORT_BODY_19 }
	else
	if(n < 22)
	if(n == 20) { INS_SORT_BODY_20 }

	*/
}

//---------------------------------------------------------------------------
// Sortowanie przez proste wstawianie bez wartownika
template<unsigned SIZE>
void CSmallSort<SIZE>::ins_sort_hybrid(CKmer<SIZE>* ptr, uint32 size)
{
	int64_t n = size;

	if (n <= 8)
	{
		ins_sort_macro(ptr, n);
		return;
	}

	int64_t i, j;
	CKmer<SIZE> x;

	ins_sort_macro(ptr, 8);

	CKmer<SIZE> *A = ptr;

	for (i = 8; i < n; i++)
	{
		x = A[i];
		j = i - 1;
		while ((j >= 0) && (x < A[j]))
		{
			A[j + 1] = A[j];
			j--;
		}
		A[j + 1] = x;
	}
}

template<unsigned SIZE>
void CSmallSort<SIZE>::shell_sort_1_7(CKmer<SIZE>* ptr, uint32 size)
{
	int i, j;
	CKmer<SIZE> x;
	int n = size;
	
	for (i = 7; i < n; i++)
	{
		j = i;
		x = ptr[i];
		while (j >= 7 && x < ptr[j - 7])
		{
			ptr[j] = ptr[j - 7];
			j -= 7;
		}
		ptr[j] = x;
	}

	for (i = 1; i < n; i++)
	{
		x = ptr[i];
		j = i - 1;
		while (j >= 0 && x < ptr[j])
		{
			ptr[j + 1] = ptr[j];
			j--;
		}
		ptr[j + 1] = x;
	}
}


template<unsigned SIZE>
void CSmallSort<SIZE>::shell_sort_1_8(CKmer<SIZE>* ptr, uint32 size)
{
	int i, j;
	CKmer<SIZE> x;
	int n = size;

	for (i = 8; i < n; i++)
	{
		j = i;
		x = ptr[i];
		while (j >= 8 && x < ptr[j - 8])
		{
			ptr[j] = ptr[j - 8];
			j -= 8;
		}
		ptr[j] = x;
	}

	for (i = 1; i < n; i++)
	{
		x = ptr[i];
		j = i - 1;
		while (j >= 0 && x < ptr[j])
		{
			ptr[j + 1] = ptr[j];
			j--;
		}
		ptr[j + 1] = x;
	}
}

template<unsigned SIZE>
void CSmallSort<SIZE>::shell_sort_1_10(CKmer<SIZE>* ptr, uint32 size)
{
	int i, j;
	CKmer<SIZE> x;
	int n = size;

	for (i = 10; i < n; i++)
	{
		j = i;
		x = ptr[i];
		while (j >= 10 && x < ptr[j - 10])
		{
			ptr[j] = ptr[j - 10];
			j -= 10;
		}
		ptr[j] = x;
	}

	for (i = 1; i < n; i++)
	{
		x = ptr[i];
		j = i - 1;
		while (j >= 0 && x < ptr[j])
		{
			ptr[j + 1] = ptr[j];
			j--;
		}
		ptr[j + 1] = x;
	}
}



#endif
	
// ***** EOF