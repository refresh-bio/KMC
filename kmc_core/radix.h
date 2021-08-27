/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.
Date   : 2019-05-19
*/
#ifndef _RADIX_H
#define _RADIX_H

#include <cassert>
#include <queue>
#include <condition_variable>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "defs.h"
#include "timer.h"
#include <thread>
#include "small_sort.h"
#include "intr_copy.h"

namespace RadixSort
{
	// Thresholds chosen experimentally. Must be extended if MAX_K > 512 !!!
	const uint64 insertion_sort_thresholds[] = { 32, 32, 32, 25, 54, 42, 42, 32, 32, 32, 32, 32, 32, 32, 32, 32 };
	const uint64 shell_sort_thresholds[] = { 32, 180, 180, 256, 134, 165, 87, 103, 103, 103, 103, 103, 103, 103, 103, 103 };
	const uint64 std_sort_thresholds[] = { 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384 };
	const uint64 small_sort_thresholds[] = { 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384 };

	static_assert(sizeof(small_sort_thresholds) / sizeof(uint64) >= MAX_K / 32, "Extend small_sort_threshold and 3 similar arrays");

	template<typename KMER_T, unsigned SIZE>
	inline void SmallSortDispatch(KMER_T* kmers, uint64 size)
	{
		CSmallSort<SIZE>::Sort(kmers, size);


		/*
		if (size <= insertion_sort_thresholds[KMER_T::KMER_SIZE])
		InsertionSortDispatch(kmers, (int)size);
		else if (size <= shell_sort_thresholds[KMER_T::KMER_SIZE])
		ShellSortDispatch(kmers, (int)size);
		else if (size <= std_sort_thresholds[KMER_T::KMER_SIZE])
		StdSortDispatch(kmers, size);*/
	}


	template<typename KMER_T>
	/*inline void ShellSortDispatch(KMER_T* kmers, int size)
	{
	int j, h;
	KMER_T x;
	int h_tab[] = { 1, 8 };
	int h_pos = 1;

	while (h_pos >= 0)
	{
	h = h_tab[h_pos--];

	for (int i = h; i < size; i++)
	{
	j = i;
	x = kmers[i];
	while (j >= h && x < kmers[j - h])
	{
	kmers[j] = kmers[j - h];
	j -= h;
	}
	kmers[j] = x;
	}
	}
	}*/

	inline void ShellSortDispatch(KMER_T* kmers, int size)
	{
		int i, j;
		KMER_T x;

		for (i = 8; i < size; i++)
		{
			j = i;
			x = kmers[i];
			while (j >= 8 && x < kmers[j - 8])
			{
				kmers[j] = kmers[j - 8];
				j -= 8;
			}
			kmers[j] = x;
		}

		for (i = 1; i < size; i++)
		{
			x = kmers[i];
			j = i - 1;
			while (j >= 0 && x < kmers[j])
			{
				kmers[j + 1] = kmers[j];
				j--;
			}
			kmers[j + 1] = x;
		}
	}


	template<typename KMER_T>
	inline void InsertionSortDispatch(KMER_T* kmers, int size)
	{
		int i, j;
		KMER_T x;

		for (i = 1; i < size; i++)
		{
			x = kmers[i];
			j = i - 1;
			while (j >= 0 && x < kmers[j])
			{
				kmers[j + 1] = kmers[j];
				j--;
			}
			kmers[j + 1] = x;
		}
	}

	template<typename KMER_T>
	inline void StdSortDispatch(KMER_T* kmers, uint64 size)
	{
		std::sort(kmers, kmers + size);
	}

	template<typename KMER_T>
	struct CRadixMSDTaskskDesc
	{
		KMER_T* kmers;
		KMER_T* tmp;
		uint64_t n_recs;
		uint32_t byte;

		CRadixMSDTaskskDesc(KMER_T* kmers, KMER_T* tmp, uint64_t n_recs, uint32_t byte) :
			kmers(kmers), tmp(tmp), n_recs(n_recs), byte(byte)
		{
		}
		bool operator<(const CRadixMSDTaskskDesc<KMER_T> &x) const
		{
			return this->n_recs < x.n_recs;
		}
	};


	template<typename KMER_T>
	class CRadixMSDTasksQueue
	{
		std::priority_queue<CRadixMSDTaskskDesc<KMER_T>> tasks;
		//	std::queue<CRadixMSDTaskskDesc<KMER_T>> tasks;
		std::condition_variable cv_pop;
		mutable std::mutex mtx;
		uint64_t tasks_in_progress = 0;
	public:
		void push(KMER_T* kmers, KMER_T* tmp, uint64_t n, uint32_t byte)
		{
			std::lock_guard<std::mutex> lck(mtx);
			tasks_in_progress++;
			tasks.emplace(kmers, tmp, n, byte);

			if (tasks.size() == 1) // was empty
				cv_pop.notify_all();
		}

		bool pop(KMER_T* &kmers, KMER_T* &tmp, uint64_t& n, uint32_t& byte)
		{
			std::unique_lock<std::mutex> lck(mtx);
			cv_pop.wait(lck, [this] { return tasks.size() || !tasks_in_progress; });
			if (!tasks_in_progress)
				return false;
			/*		kmers = tasks.front().kmers;
			tmp = tasks.front().tmp;
			n = tasks.front().n_recs;
			byte = tasks.front().byte;*/
			kmers = tasks.top().kmers;
			tmp = tasks.top().tmp;
			n = tasks.top().n_recs;
			byte = tasks.top().byte;

			tasks.pop();

			return true;
		}

		void notify_task_finished()
		{
			std::lock_guard<std::mutex> lck(mtx);
			--tasks_in_progress;
			if (!tasks_in_progress)
				cv_pop.notify_all();
		}
	};

	template <typename KMER_T, typename COUNTER_TYPE, unsigned SIZE>
	class CRadixSorterMSD
	{
		CRadixMSDTasksQueue<KMER_T>& tasks_queue;
		CMemoryPool* pmm_radix_buf;
		uint64 use_queue_min_recs = 0;

		void Sort(KMER_T* kmers, KMER_T* tmp, uint64_t n_recs, uint32_t byte)
		{
			uint8_t* ptr = (uint8_t*)kmers + byte;
			ALIGN_ARRAY COUNTER_TYPE globalHisto[256] = {};
			ALIGN_ARRAY COUNTER_TYPE copy_globalHisto[257] = {};

			/*for (uint64_t i = 0; i < n_recs; i++)
			{
			globalHisto[*ptr]++;
			ptr += sizeof(KMER_T);
			}*/
			switch (n_recs % 4)
			{
			case 3:
				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);
			case 2:
				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);
			case 1:
				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);
			}
			for (uint64 i = 0; i < n_recs / 4; ++i)
			{
				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);

				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);

				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);

				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);
			}

			COUNTER_TYPE prevSum = 0;
			for (int i = 0; i < 256; ++i)
			{
				COUNTER_TYPE temp = globalHisto[i];
				globalHisto[i] = prevSum;
				copy_globalHisto[i] = prevSum;
				prevSum += temp;
			}
			copy_globalHisto[256] = (COUNTER_TYPE)n_recs;

			KMER_T* src = kmers;
			ptr = (uint8_t*)kmers + byte;

			if (n_recs * sizeof(KMER_T) < (1ull << 17))
			{
				/*for (uint64_t i = 0; i < n_recs; ++i)
				{
				tmp[globalHisto[*ptr]] = src[i];
				globalHisto[*ptr]++;
				ptr += sizeof(KMER_T);
				}*/


				switch (n_recs % 4)
				{
				case 3:
					tmp[globalHisto[*ptr]] = src[(n_recs % 4) - 3];
					globalHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				case 2:
					tmp[globalHisto[*ptr]] = src[(n_recs % 4) - 2];
					globalHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				case 1:
					tmp[globalHisto[*ptr]] = src[(n_recs % 4) - 1];
					globalHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				}

				for (uint64 i = n_recs % 4; i < n_recs; i += 4)
				{
					tmp[globalHisto[*ptr]] = src[i];
					globalHisto[*ptr]++;
					ptr += sizeof(KMER_T);

					tmp[globalHisto[*ptr]] = src[i + 1];
					globalHisto[*ptr]++;
					ptr += sizeof(KMER_T);

					tmp[globalHisto[*ptr]] = src[i + 2];
					globalHisto[*ptr]++;
					ptr += sizeof(KMER_T);

					tmp[globalHisto[*ptr]] = src[i + 3];
					globalHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				}

			}
			else
			{
				//const int32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(KMER_T) / 8];

				constexpr uint32_t BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(KMER_T) / 8];
				constexpr uint32_t BUFFER_WIDTH_IN_128BIT_WORDS = BUFFER_WIDTH * sizeof(KMER_T) / 16;
				constexpr uint32_t BUFFER_16B_ALIGNED = sizeof(KMER_T) % 16 == 0;


				uchar* raw_buffer;
				pmm_radix_buf->reserve(raw_buffer);
				uchar* buffer = raw_buffer;
				while ((uint64_t)buffer % ALIGNMENT)
					++buffer;
				KMER_T *Buffer = (KMER_T*)buffer;

				uint8_t byteValue = 0;
				int index_x = 0;
				//for (uint64_t i = 0; i < n_recs; ++i)
				//{
				//	byteValue = *ptr;

				//	index_x = globalHisto[byteValue] % BUFFER_WIDTH;

				//	Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i];

				//	globalHisto[byteValue]++;

				//	if (index_x == (BUFFER_WIDTH - 1))
				//		memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));

				//	ptr += sizeof(KMER_T);
				//} //end_for


				switch (n_recs % 4)
				{
				case 3:
					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n_recs % 4) - 3];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				case 2:
					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n_recs % 4) - 2];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				case 1:
					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n_recs % 4) - 1];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				}

				for (uint64_t i = n_recs % 4; i < n_recs; i += 4)
				{
					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 1];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 2];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 3];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))												
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

				}



				int64_t elemInBuffer;
				int64_t index_stop;
				int64_t index_start;
				int64_t elemWrittenIntoBuffer;

				for (uint32_t private_i = 0; private_i < 256; private_i++)
				{
					index_stop = globalHisto[private_i] % BUFFER_WIDTH;
					index_start = copy_globalHisto[private_i] % BUFFER_WIDTH;
					elemWrittenIntoBuffer = globalHisto[private_i] - copy_globalHisto[private_i];

					if ((index_stop - elemWrittenIntoBuffer) <= 0)
						elemInBuffer = index_stop;
					else
						elemInBuffer = index_stop - index_start;

					if (elemInBuffer != 0)						
						IntrCopy64fun(&tmp[globalHisto[private_i] - elemInBuffer],
							&Buffer[private_i * BUFFER_WIDTH + (globalHisto[private_i] - elemInBuffer) % BUFFER_WIDTH], elemInBuffer * sizeof(KMER_T) / 8);
				}

				pmm_radix_buf->free(raw_buffer);
			}


			if (byte > 0)
			{
				for (int i = 0; i < 256; i++)
				{
					uint64_t new_n = copy_globalHisto[i + 1] - copy_globalHisto[i];

					if (new_n <= small_sort_thresholds[KMER_T::KMER_SIZE])
					{
						SmallSortDispatch<KMER_T, SIZE>(tmp + copy_globalHisto[i], new_n);
						if (byte % 2 != 0)
						{
							for (COUNTER_TYPE j = copy_globalHisto[i]; j < copy_globalHisto[i] + (COUNTER_TYPE)new_n; ++j)
								kmers[j] = tmp[j];
						}
					}
					else
					{
						if (new_n >= use_queue_min_recs)
							tasks_queue.push(tmp + copy_globalHisto[i], kmers + copy_globalHisto[i], new_n, byte - 1);
						else
							Sort(tmp + copy_globalHisto[i], kmers + copy_globalHisto[i], new_n, byte - 1);
					}

				}
			}

		}
	public:
		CRadixSorterMSD(CRadixMSDTasksQueue<KMER_T>& tasks_queue, CMemoryPool* pmm_radix_buf, uint64 use_queue_min_recs) :
			tasks_queue(tasks_queue),
			pmm_radix_buf(pmm_radix_buf),
			use_queue_min_recs(use_queue_min_recs)
		{
		}
		void operator()()
		{
			KMER_T* kmers;
			KMER_T* tmp;
			uint64_t n_recs;
			uint32 byte;
			while (tasks_queue.pop(kmers, tmp, n_recs, byte))
			{
				Sort(kmers, tmp, n_recs, byte);
				tasks_queue.notify_task_finished();
			}
		}
	};

	template<typename KMER_T, typename COUNTER_TYPE, unsigned SIZE>
	void RadixSortMSD_impl(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf, bool is_first_level,
		uint64 is_big_threshold, uint64 n_total_recs)
	{
		if (n_recs <= small_sort_thresholds[KMER_T::KMER_SIZE])
		{
			SmallSortDispatch<KMER_T, SIZE>(kmers, n_recs);
			if (byte % 2 == 0)
			{
				for (uint64 j = 0; j < n_recs; ++j)
					tmp[j] = kmers[j];
			}
			return;
		}

		//	cout << n_recs << " : " << n_threads << "  :  " << byte << endl;

		uint64 per_thread = n_recs / n_threads;

		std::vector<std::thread> threads;
		std::vector<std::array<COUNTER_TYPE, 256>> histos(n_threads);
		ALIGN_ARRAY COUNTER_TYPE globalHisto[256] = {};
		for (uint32_t th_id = 0; th_id < n_threads; ++th_id)
		{
			threads.push_back(std::thread([th_id, kmers, n_recs, n_threads, per_thread, &histos, byte]
			{
				ALIGN_ARRAY COUNTER_TYPE myHisto[256] = { 0 };

				uint8_t* ptr = (uint8_t*)(kmers + th_id*per_thread) + byte;
				uint64_t n = per_thread;
				if (th_id == n_threads - 1)
					n = n_recs - th_id*per_thread;

				/*for (uint64_t i = 0; i < n; i++)
				{
				myHisto[*ptr]++;
				ptr += sizeof(KMER_T);
				}*/

				switch (n % 4)
				{
				case 3:
					myHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				case 2:
					myHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				case 1:
					myHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				}
				for (uint64 i = 0; i < n / 4; ++i)
				{
					myHisto[*ptr]++;
					ptr += sizeof(KMER_T);

					myHisto[*ptr]++;
					ptr += sizeof(KMER_T);

					myHisto[*ptr]++;
					ptr += sizeof(KMER_T);

					myHisto[*ptr]++;
					ptr += sizeof(KMER_T);
				}

				for (uint32_t i = 0; i < 256; ++i)
				{
					histos[th_id][i] = myHisto[i];
				}
			}));
		}

		for (auto& t : threads)
			t.join();
		threads.clear();


		// ***** collecting counters
		for (int i = 0; i < 256; ++i)
		{
			COUNTER_TYPE prevSum = 0;
			for (uint32_t n = 0; n < n_threads; n++)
			{
				COUNTER_TYPE temp = histos[n][i];
				histos[n][i] = prevSum;
				prevSum += temp;
			}
			globalHisto[i] = prevSum;
		}

		COUNTER_TYPE prevSum = 0;
		for (int i = 0; i < 256; ++i)
		{
			COUNTER_TYPE temp = globalHisto[i];
			globalHisto[i] = prevSum;
			prevSum += temp;
		}

		for (uint32_t n = 0; n < n_threads; ++n)
		{
			for (int i = 0; i < 256; i++)
			{
				histos[n][i] += globalHisto[i];
			}
		}

		std::vector<uchar*> _raw_buffers(n_threads);
		std::vector<std::array<COUNTER_TYPE, 256>> threads_histos(n_threads);

		for (uint32_t th_id = 0; th_id < n_threads; ++th_id)
		{
			threads.push_back(std::thread([th_id, kmers, tmp, n_recs, n_threads, per_thread, byte, &histos, &_raw_buffers, &threads_histos, pmm_radix_buf]
			{
				ALIGN_ARRAY COUNTER_TYPE myHisto[256];

				for (int i = 0; i < 256; ++i)
					myHisto[i] = histos[th_id][i];

				uint8_t* ptr = (uint8_t*)(kmers + th_id*per_thread) + byte;
				uint64_t n = per_thread;
				if (th_id == n_threads - 1)
					n = n_recs - th_id*per_thread;

				KMER_T* src = kmers + th_id*per_thread;

				//const int32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(KMER_T) / 8];
				constexpr uint32_t BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(KMER_T) / 8];
				constexpr uint32_t BUFFER_WIDTH_IN_128BIT_WORDS = BUFFER_WIDTH * sizeof(KMER_T) / 16;
				constexpr uint32_t BUFFER_16B_ALIGNED = sizeof(KMER_T) % 16 == 0;

				uchar* raw_buffer;
				pmm_radix_buf->reserve(raw_buffer);
				_raw_buffers[th_id] = raw_buffer;
				uchar* buffer = raw_buffer;
				while ((uint64_t)buffer % ALIGNMENT)
					++buffer;
				KMER_T *Buffer = (KMER_T*)buffer;

				uint8_t byteValue = 0;
				int index_x = 0;
				//for (uint64_t i = 0; i < n; ++i)
				//{
				//	byteValue = *ptr;

				//	index_x = myHisto[byteValue] % BUFFER_WIDTH;

				//	Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i];

				//	myHisto[byteValue]++;

				//	if (index_x == (BUFFER_WIDTH - 1))
				//		memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));

				//	ptr += sizeof(KMER_T);
				//} //end_for

				switch (n % 4)
				{
				case 3:
					byteValue = *ptr;
					index_x = myHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 3];
					myHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				case 2:
					byteValue = *ptr;
					index_x = myHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 2];
					myHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				case 1:
					byteValue = *ptr;
					index_x = myHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 1];
					myHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				}

				for (uint64 i = n % 4; i < n; i += 4)
				{
					byteValue = *ptr;
					index_x = myHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i];
					myHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = myHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 1];
					myHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = myHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 2];
					myHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = myHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 3];
					myHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))						
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				}


				for (uint32 i = 0; i < 256; ++i)
					threads_histos[th_id][i] = myHisto[i];
			}));
		}
		for (auto& t : threads)
			t.join();
		threads.clear();


		for (uint32_t th_id = 0; th_id < n_threads; ++th_id)
		{
			threads.push_back(std::thread([th_id, kmers, tmp, n_recs, n_threads, per_thread, byte, &histos, &_raw_buffers, &threads_histos, pmm_radix_buf]
			{
				ALIGN_ARRAY COUNTER_TYPE myHisto[256];
				for (int i = 0; i < 256; ++i)
					myHisto[i] = threads_histos[th_id][i];
				uchar* raw_buffer = _raw_buffers[th_id];
				uchar* buffer = raw_buffer;
				while ((uint64_t)buffer % ALIGNMENT)
					++buffer;
				KMER_T *Buffer = (KMER_T*)buffer;

				int64_t elemInBuffer;
				int64_t index_stop;
				int64_t index_start;
				int64_t elemWrittenIntoBuffer;

				const uint32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(KMER_T) / 8];

				for (uint32_t private_i = 0; private_i < 256; private_i++)
				{
					index_stop = myHisto[private_i] % BUFFER_WIDTH;
					index_start = histos[th_id][private_i] % BUFFER_WIDTH;
					elemWrittenIntoBuffer = myHisto[private_i] - histos[th_id][private_i];

					if ((index_stop - elemWrittenIntoBuffer) <= 0)
						elemInBuffer = index_stop;
					else
						elemInBuffer = index_stop - index_start;

					if (elemInBuffer != 0)						
						IntrCopy64fun(&tmp[myHisto[private_i] - elemInBuffer],
							&Buffer[private_i * BUFFER_WIDTH + (myHisto[private_i] - elemInBuffer) % BUFFER_WIDTH], elemInBuffer * sizeof(KMER_T) / 8);
				}
				pmm_radix_buf->free(raw_buffer);
			}));
		}
		for (auto& t : threads)
			t.join();
		threads.clear();

		if (byte > 0)
		{
			CRadixMSDTasksQueue<KMER_T> tasks_queue;

			KMER_T* kmers_ptr = kmers;
			KMER_T* ptr = tmp;

			/*		if (n_threads <= 4)
			is_big_threshold = n_recs; //for 4 or less threads do not extract big bins
			*/

			std::vector<std::tuple<KMER_T*, KMER_T*, uint64>> big_bins;
			uint64_t n_rec_in_big_bins = 0;

			for (uint32_t i = 1; i < 256; ++i)
			{
				uint64_t n = globalHisto[i] - globalHisto[i - 1];
				if (n > 0)
				{
					if (n > is_big_threshold)
						if (!is_first_level)
							RadixSortMSD_impl<KMER_T, COUNTER_TYPE, SIZE>(ptr, kmers_ptr, n, byte - 1, n_threads, pmm_radix_buf, false, is_big_threshold, n_total_recs);
						else
						{
							big_bins.push_back(std::make_tuple(ptr, kmers_ptr, n));
							n_rec_in_big_bins += n;
						}
					else
						tasks_queue.push(ptr, kmers_ptr, n, byte - 1);
				}
				ptr += n;
				kmers_ptr += n;
			}
			uint64_t n = n_recs - globalHisto[255];
			if (n > 0)
			{
				if (n > is_big_threshold)
					if (!is_first_level)
						RadixSortMSD_impl<KMER_T, COUNTER_TYPE, SIZE>(ptr, kmers_ptr, n, byte - 1, n_threads, pmm_radix_buf, false, is_big_threshold, n_total_recs);
					else
					{
						big_bins.push_back(std::make_tuple(ptr, kmers_ptr, n));
						n_rec_in_big_bins += n;
					}
				else
					tasks_queue.push(ptr, kmers_ptr, n, byte - 1);
			}
			ptr += n;
			kmers_ptr += n;

			sort(big_bins.begin(), big_bins.end(), [](std::tuple<KMER_T*, KMER_T*, uint64> x, std::tuple<KMER_T*, KMER_T*, uint64> y) {return get<2>(x) > get<2>(y); });

			//		uint32 n_threads_for_big_bins = 2 * n_threads / 3;
			uint32 n_threads_for_big_bins = uint32(ceil(n_threads * n_rec_in_big_bins * 5.0 / (4 * n_total_recs)));
			if (n_threads_for_big_bins > n_threads)
				n_threads_for_big_bins = n_threads;

			uint32 n_threads_for_small_bins_running;
			uint32 n_threads_for_small_bins = n_threads - n_threads_for_big_bins;
			//		cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;

			std::vector<CRadixSorterMSD<KMER_T, COUNTER_TYPE, SIZE>*> sorters;
			for (n_threads_for_small_bins_running = 0; n_threads_for_small_bins_running < n_threads_for_small_bins; ++n_threads_for_small_bins_running)
			{
				sorters.push_back(new CRadixSorterMSD<KMER_T, COUNTER_TYPE, SIZE>(tasks_queue, pmm_radix_buf, n_recs / 4096));
				threads.push_back(std::thread(std::ref(*sorters.back())));
			}
			//		cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;

			//process big bins (only in first radix pass, for later big_bins.size() equals 0)
			for (auto& big_bin : big_bins)
			{
				RadixSortMSD_impl<KMER_T, COUNTER_TYPE, SIZE>(get<0>(big_bin), get<1>(big_bin), get<2>(big_bin), byte - 1, n_threads_for_big_bins, pmm_radix_buf, false,
					is_big_threshold, n_total_recs);

				/*			n_rec_in_big_bins -= get<2>(big_bin);
				n_threads_for_big_bins = n_threads * n_rec_in_big_bins * 5 / (4 * n_total_recs);
				if (n_threads_for_big_bins > n_threads)
				n_threads_for_big_bins = n_threads;
				n_threads_for_small_bins = n_threads - n_threads_for_big_bins;
				cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;

				for (; n_threads_for_small_bins_running < n_threads_for_small_bins; ++n_threads_for_small_bins_running)
				{
				sorters.push_back(new CRadixSorterMSD<KMER_T, COUNTER_TYPE>(tasks_queue, pmm_radix_buf, n_recs / 4096));
				threads.push_back(std::thread(std::ref(*sorters.back())));
				}
				cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;*/
			}
			//		cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;
			//now i can use threads left after processing big bins to process small ones
			for (; n_threads_for_small_bins_running < n_threads; ++n_threads_for_small_bins_running)
			{
				sorters.push_back(new CRadixSorterMSD<KMER_T, COUNTER_TYPE, SIZE>(tasks_queue, pmm_radix_buf, n_recs / 4096));
				threads.push_back(std::thread(std::ref(*sorters.back())));
			}

			for (auto& t : threads)
				t.join();
			for (auto s : sorters)
				delete s;
		}
	}


	template<typename KMER_T, unsigned SIZE>
	void RadixSortMSD(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf)
	{
		if (n_recs >= (1ull << 31))
			RadixSortMSD_impl<KMER_T, int64, SIZE>(kmers, tmp, n_recs, byte, n_threads, pmm_radix_buf, true, 2 * n_recs / n_threads, n_recs);
		else
			RadixSortMSD_impl<KMER_T, int32, SIZE>(kmers, tmp, n_recs, byte, n_threads, pmm_radix_buf, true, 2 * n_recs / n_threads, n_recs);
	}

}
#endif

// ***** EOF
