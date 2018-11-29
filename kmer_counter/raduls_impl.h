/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _RADULS_IMPL_H
#define _RADULS_IMPL_H

#include <cassert>
#include <queue>
#include <condition_variable>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "defs.h"
#include "kmer.h"
#include "timer.h"
#include <thread>
#include "first_dispatch.h"
#include "intr_copy.h"
#include "raduls.h"

#define IS_NARROW(x, y)	((x) < (y) * 16)

//#define USE_TIMERS
#include "stdafx.h"
#include "raduls_impl.h"

namespace RadulsSort
{
	// Thresholds chosen experimentally. Must be extended if MAX_K > 512 !!!
	const uint64 insertion_sort_thresholds[] = { 32, 32, 32, 25, 54, 42, 42, 32, 32, 32, 32, 32, 32, 32, 32, 32 };
	const uint64 shell_sort_thresholds[] = { 32, 180, 180, 256, 134, 165, 87, 103, 103, 103, 103, 103, 103, 103, 103, 103 };
	const uint64 std_sort_thresholds[] = { 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384 };
	const uint64 small_sort_thresholds[] = { 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384 };
	const uint64 wide_small_sort_thresholds[] = { 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32 };

	static_assert(sizeof(small_sort_thresholds) / sizeof(uint64) >= MAX_K / 32, "Extend small_sort_threshold and 3 similar arrays");



	template<typename KMER_T>
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
	inline void SmallSortDispatch(KMER_T* kmers, uint64 size)
	{
		if (size <= insertion_sort_thresholds[KMER_T::KMER_SIZE])
			InsertionSortDispatch(kmers, (int)size);
		else if (size <= shell_sort_thresholds[KMER_T::KMER_SIZE])
			ShellSortDispatch(kmers, (int)size);
		else if (size <= std_sort_thresholds[KMER_T::KMER_SIZE])
			StdSortDispatch(kmers, size);
	}


	template<typename KMER_T>
	struct CRadixMSDTaskskDesc
	{
		KMER_T* kmers;
		KMER_T* tmp;
		uint64_t n_recs;
		uint32_t byte;
		bool is_narrow;

		CRadixMSDTaskskDesc(KMER_T* kmers, KMER_T* tmp, uint64_t n_recs, uint32_t byte, bool is_narrow) :
			kmers(kmers), tmp(tmp), n_recs(n_recs), byte(byte), is_narrow(is_narrow)
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
		std::condition_variable cv_pop;
		mutable std::mutex mtx;
		uint64_t tasks_in_progress = 0;
	public:
		void push(KMER_T* kmers, KMER_T* tmp, uint64_t n, uint32_t byte, bool is_narrow)
		{
			std::lock_guard<std::mutex> lck(mtx);
			tasks_in_progress++;
			tasks.emplace(kmers, tmp, n, byte, is_narrow);

			if (tasks.size() == 1) // was empty
				cv_pop.notify_all();
		}

		bool pop(KMER_T* &kmers, KMER_T* &tmp, uint64_t& n, uint32_t& byte, bool &is_narrow)
		{
			std::unique_lock<std::mutex> lck(mtx);
			cv_pop.wait(lck, [this] { return tasks.size() || !tasks_in_progress; });
			if (!tasks_in_progress)
				return false;

			kmers = tasks.top().kmers;
			tmp = tasks.top().tmp;
			n = tasks.top().n_recs;
			byte = tasks.top().byte;
			is_narrow = tasks.top().is_narrow;

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

	template <typename KMER_T, typename COUNTER_TYPE>
	class CRadixSorterMSD
	{
		CRadixMSDTasksQueue<KMER_T>& tasks_queue;
		CMemoryPool* pmm_radix_buf;
		uint64 use_queue_min_recs = 0;

		void Sort(KMER_T* kmers, KMER_T* tmp, uint64_t n_recs, uint32_t byte, bool is_narrow)
		{
#ifdef MEASURE_TIMES
			CThreadWatch tw;
			tw.startTimer();
#endif

			uint8_t* ptr = (uint8_t*)kmers + byte;
			ALIGN_ARRAY COUNTER_TYPE globalHisto[256] = {};
			ALIGN_ARRAY COUNTER_TYPE copy_globalHisto[257] = {};

			//-------------------------------------------------------------------------------
			//the following loop is unrolled below:
			//-------------------------------------------------------------------------------
			//for (uint64_t i = 0; i < n_recs; i++)
			//{
			//globalHisto[*ptr]++;
			//ptr += sizeof(KMER_T);
			//}
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
			//-------------------------------------------------------------------------------
			//the end of the unrolled loop
			//-------------------------------------------------------------------------------

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

#ifdef MEASURE_TIMES
			tw.stopTimer();

			times_byte_total[byte] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
			if (n_recs * sizeof(KMER_T) >= (1ull << 17))
				times_satish_stages[byte][0] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
#endif

			if (n_recs * sizeof(KMER_T) < (1ull << 17))
			{
#ifdef MEASURE_TIMES
				tw.startTimer();
#endif
				//-------------------------------------------------------------------------------
				//the following loop is unrolled below:
				//-------------------------------------------------------------------------------

				//for (uint64_t i = 0; i < n_recs; ++i)
				//{
				//tmp[globalHisto[*ptr]] = src[i];
				//globalHisto[*ptr]++;
				//ptr += sizeof(KMER_T);
				//}


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
				//-------------------------------------------------------------------------------
				//the end of the unrolled loop
				//-------------------------------------------------------------------------------
#ifdef MEASURE_TIMES
				tw.stopTimer();

				times_byte_total[byte] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
#endif
			}
			else
			{
#ifdef MEASURE_TIMES
				tw.startTimer();
#endif
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

				//-------------------------------------------------------------------------------
				//the following loop is unrolled below
				//-------------------------------------------------------------------------------

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
						//					memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				case 2:
					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n_recs % 4) - 2];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))
						//					memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				case 1:
					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n_recs % 4) - 1];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))
						//					memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
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
						//					memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 1];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))
						//					memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 2];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))
						//					memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);

					byteValue = *ptr;
					index_x = globalHisto[byteValue] % BUFFER_WIDTH;
					Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 3];
					globalHisto[byteValue]++;
					if (index_x == (BUFFER_WIDTH - 1))
						//					memcpy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
						IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[globalHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);
					ptr += sizeof(KMER_T);
				}
#ifdef MEASURE_TIMES
				tw.stopTimer();

				times_byte_total[byte] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
				times_satish_stages[byte][1] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
#endif
				//-------------------------------------------------------------------------------
				//the end of the unrolled loop
				//-------------------------------------------------------------------------------
#ifdef MEASURE_TIMES
				tw.startTimer();
#endif
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
						//					memcpy(&tmp[globalHisto[private_i] - elemInBuffer], &Buffer[private_i * BUFFER_WIDTH + (globalHisto[private_i] - elemInBuffer) % BUFFER_WIDTH], (elemInBuffer)*sizeof(KMER_T));
						IntrCopy64fun(&tmp[globalHisto[private_i] - elemInBuffer],
							&Buffer[private_i * BUFFER_WIDTH + (globalHisto[private_i] - elemInBuffer) % BUFFER_WIDTH], elemInBuffer * sizeof(KMER_T) / 8);
				}

				pmm_radix_buf->free(raw_buffer);
#ifdef MEASURE_TIMES
				tw.stopTimer();

				times_byte_total[byte] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
				times_satish_stages[byte][2] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
#endif
			}

			if (byte > 0)
			{
				bool must_copy_tmp = byte % 2 != 0;

				uint64_t narrow_small_sort_threshold = small_sort_thresholds[KMER_T::KMER_SIZE];
				uint64_t wide_small_sort_threshold = wide_small_sort_thresholds[KMER_T::KMER_SIZE];

				for (int i = 0; i < 256; i++)
				{
					uint64_t new_n = copy_globalHisto[i + 1] - copy_globalHisto[i];

					if (new_n <= wide_small_sort_threshold || (is_narrow && new_n <= narrow_small_sort_threshold))
					{
						SmallSortDispatch(tmp + copy_globalHisto[i], new_n);
						if (must_copy_tmp)
							for (COUNTER_TYPE j = copy_globalHisto[i]; j < copy_globalHisto[i] + (COUNTER_TYPE)new_n; ++j)
								kmers[j] = tmp[j];
					}
					else
					{
						if (new_n >= use_queue_min_recs)
							tasks_queue.push(tmp + copy_globalHisto[i], kmers + copy_globalHisto[i], new_n, byte - 1, IS_NARROW(n_recs, new_n));
						else
							Sort(tmp + copy_globalHisto[i], kmers + copy_globalHisto[i], new_n, byte - 1, IS_NARROW(n_recs, new_n));
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
			bool is_narrow;

			while (tasks_queue.pop(kmers, tmp, n_recs, byte, is_narrow))
			{
				Sort(kmers, tmp, n_recs, byte, is_narrow);
				tasks_queue.notify_task_finished();
			}
		}
	};


	template<typename KMER_T, typename COUNTER_TYPE>
	void RadixSortMSD_impl(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf, bool is_first_level,
		uint64 is_big_threshold, uint64 n_total_recs)
	{
		uint64_t current_small_sort_threshold = small_sort_thresholds[KMER_T::KMER_SIZE];

		if (n_recs <= current_small_sort_threshold)
		{
			SmallSortDispatch(kmers, n_recs);
			if (byte % 2 == 0)
			{
				for (uint64 j = 0; j < n_recs; ++j)
					tmp[j] = kmers[j];
			}
			return;
		}

		CRangeQueue my_buffer(MAGIC_NUMBER * n_threads, n_recs);

		uint64 per_thread = n_recs / n_threads;

#ifdef USE_TIMERS
		CStopWatch sw;
		sw.startTimer();
#endif

		std::vector<std::thread> threads;
		//	std::vector<ALIGN_ARRAY COUNTER_TYPE[256]> histos(MAGIC_NUMBER * n_threads);
		std::vector<std::array<COUNTER_TYPE, 256>> histos(MAGIC_NUMBER * n_threads);
		ALIGN_ARRAY COUNTER_TYPE globalHisto[256] = {};
		for (uint32_t th_id = 0; th_id < n_threads; ++th_id)
		{
			threads.push_back(std::thread(pierwsze_kolko_etap1<KMER_T, COUNTER_TYPE>, th_id, kmers, n_recs, n_threads, per_thread,
				ref(histos), byte, std::ref(my_buffer)));
		}

		for (auto& t : threads)
			t.join();
		threads.clear();

#ifdef USE_TIMERS
		sw.stopTimer();
		cout << "1: " << sw.getElapsedTime() << endl; 	fflush(stdout);
		sw.startTimer();
#endif

		// ***** collecting counters
		for (int i = 0; i < 256; ++i)
		{
			COUNTER_TYPE prevSum = 0;
			for (uint32_t n = 0; n < MAGIC_NUMBER * n_threads; n++) //<<<---------------
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

		for (uint32_t n = 0; n < MAGIC_NUMBER * n_threads; ++n) //<<<------------------
		{
			for (int i = 0; i < 256; i++)
			{
				histos[n][i] += globalHisto[i];
			}
		}

		my_buffer.reset_indices();

		std::vector<uchar*> _raw_buffers(MAGIC_NUMBER * n_threads);
		//	std::vector<ALIGN_ARRAY COUNTER_TYPE[256]> threads_histos(MAGIC_NUMBER * n_threads);
		std::vector<std::array<COUNTER_TYPE, 256>> threads_histos(MAGIC_NUMBER * n_threads);

		for (uint32_t th_id = 0; th_id < n_threads; ++th_id)
		{
			threads.push_back(std::thread(pierwsze_kolko_etap2<KMER_T, COUNTER_TYPE>, th_id, kmers, tmp,
				n_recs, n_threads, per_thread, byte,
				ref(histos), ref(_raw_buffers), ref(threads_histos), pmm_radix_buf, ref(my_buffer)));
		}
		for (auto& t : threads)
			t.join();
		threads.clear();

#ifdef USE_TIMERS
		sw.stopTimer();
		cout << "2: " << sw.getElapsedTime() << endl; 	fflush(stdout);
		sw.startTimer();
#endif

		my_buffer.reset_indices();
		for (uint32_t th_id = 0; th_id < n_threads; ++th_id)
		{
			threads.push_back(std::thread(pierwsze_kolko_etap3<KMER_T, COUNTER_TYPE>, th_id, kmers, tmp,
				n_recs, n_threads, per_thread, byte,
				ref(histos), ref(_raw_buffers), ref(threads_histos), pmm_radix_buf, ref(my_buffer)));
		}
		for (auto& t : threads)
			t.join();
		threads.clear();

#ifdef USE_TIMERS
		sw.stopTimer();
		cout << "3: " << sw.getElapsedTime() << endl; 	fflush(stdout);
#endif

		if (byte > 0)
		{
			//---------------------------

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
							RadixSortMSD_impl<KMER_T, COUNTER_TYPE>(ptr, kmers_ptr, n, byte - 1, n_threads, pmm_radix_buf, false, is_big_threshold, n_total_recs);
						else
						{
							big_bins.push_back(std::make_tuple(ptr, kmers_ptr, n));
							n_rec_in_big_bins += n;
						}
					else
						tasks_queue.push(ptr, kmers_ptr, n, byte - 1, IS_NARROW(n_recs, n));
					//					tasks_queue.push(ptr, kmers_ptr, n, byte - 1, 0);
				}
				ptr += n;
				kmers_ptr += n;
			}
			uint64_t n = n_recs - globalHisto[255];
			if (n > 0)
			{
				if (n > is_big_threshold)
					if (!is_first_level)
						RadixSortMSD_impl<KMER_T, COUNTER_TYPE>(ptr, kmers_ptr, n, byte - 1, n_threads, pmm_radix_buf, false, is_big_threshold, n_total_recs);
					else
					{
						big_bins.push_back(std::make_tuple(ptr, kmers_ptr, n));
						n_rec_in_big_bins += n;
					}
				else
					tasks_queue.push(ptr, kmers_ptr, n, byte - 1, IS_NARROW(n, n_recs));
				//				tasks_queue.push(ptr, kmers_ptr, n, byte - 1, 0);
			}
			ptr += n;
			kmers_ptr += n;

			sort(big_bins.begin(), big_bins.end(), [](std::tuple<KMER_T*, KMER_T*, uint64> x, std::tuple<KMER_T*, KMER_T*, uint64> y) { return get<2>(x) > get<2>(y); });

			//		uint32 n_threads_for_big_bins = 2 * n_threads / 3;
			uint32 n_threads_for_big_bins = uint32(ceil(n_threads * n_rec_in_big_bins * 5.0 / (4 * n_total_recs)));
			if (n_threads_for_big_bins > n_threads)
				n_threads_for_big_bins = n_threads;

			uint32 n_threads_for_small_bins_running;
			uint32 n_threads_for_small_bins = n_threads - n_threads_for_big_bins;
			//		cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;

			std::vector<CRadixSorterMSD<KMER_T, COUNTER_TYPE>*> sorters;
			for (n_threads_for_small_bins_running = 0; n_threads_for_small_bins_running < n_threads_for_small_bins; ++n_threads_for_small_bins_running)
			{
				sorters.push_back(new CRadixSorterMSD<KMER_T, COUNTER_TYPE>(tasks_queue, pmm_radix_buf, n_recs / 4096));
				threads.push_back(std::thread(std::ref(*sorters.back())));
			}
			//		cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;

			//process big bins (only in first radix pass, for later big_bins.size() equals 0)
			for (auto& big_bin : big_bins)
			{
				RadixSortMSD_impl<KMER_T, COUNTER_TYPE>(get<0>(big_bin), get<1>(big_bin), get<2>(big_bin), byte - 1, n_threads_for_big_bins, pmm_radix_buf, false,
					is_big_threshold, n_total_recs);
			}
			//		cout << n_threads_for_small_bins_running << " " << n_threads_for_small_bins << " " << n_threads_for_big_bins << " " << n_recs << " " << n_rec_in_big_bins << endl;
			//now i can use threads left after processing big bins to process small ones
			for (; n_threads_for_small_bins_running < n_threads; ++n_threads_for_small_bins_running)
			{
				sorters.push_back(new CRadixSorterMSD<KMER_T, COUNTER_TYPE>(tasks_queue, pmm_radix_buf, n_recs / 4096));
				threads.push_back(std::thread(std::ref(*sorters.back())));
			}

			for (auto& t : threads)
				t.join();
			for (auto s : sorters)
				delete s;
			//---------------------------

			//----------------------------
		}
	}

#if defined(__AVX2__)
#define RADULS_RADIX_SORT_FUNNAME RadixSortMSD_AVX2
#elif defined (__AVX__)
#define RADULS_RADIX_SORT_FUNNAME RadixSortMSD_AVX
#elif defined(__SSE4_1__)
#define RADULS_RADIX_SORT_FUNNAME RadixSortMSD_SSE41
#elif defined(__SSE2__)
#define RADULS_RADIX_SORT_FUNNAME RadixSortMSD_SSE2
#endif


	template<typename KMER_T>
	void RADULS_RADIX_SORT_FUNNAME(KMER_T* kmers, KMER_T* tmp, uint64 n_recs, uint32 byte, uint32 n_threads, CMemoryPool* pmm_radix_buf)
	{
		if (n_recs >= (1ull << 31))
			RadixSortMSD_impl<KMER_T, int64>(kmers, tmp, n_recs, byte, n_threads, pmm_radix_buf, true, 2 * n_recs / (3 * n_threads), n_recs);
		else
			RadixSortMSD_impl<KMER_T, int32>(kmers, tmp, n_recs, byte, n_threads, pmm_radix_buf, true, 2 * n_recs / (3 * n_threads), n_recs);
	}

	template<unsigned SIZE>
	class InstantiateTempl
	{
		friend class InstantiateTempl<SIZE + 1>;
		void inst()
		{
			volatile auto ptr = RADULS_RADIX_SORT_FUNNAME<CKmer<SIZE>>;
			(void)ptr; //suppress `unused` warning
			InstantiateTempl<SIZE - 1>().inst();
		}
	};

	template<>
	class InstantiateTempl<0>
	{
		friend class InstantiateTempl<0 + 1>;
		void inst()
		{
		}
	};
	template class InstantiateTempl<KMER_WORDS>;
}

#endif

// ***** EOF