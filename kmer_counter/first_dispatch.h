/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
*/
#ifndef _FIRST_DISPATCH_H
#define _FIRST_DISPATCH_H

#include "queues.h"
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
#include <array>
#include "intr_copy.h"



class CRangeQueue
{
	std::vector<std::tuple<uint64, uint64, uint32>> range_queue;
	std::mutex m;
	uint32 cur_idx;
	bool done;

public:
	CRangeQueue(uint32 parts, uint64 num_rec)
	{
		uint64 delta = num_rec / parts;
		uint64 N1 = 0, N2 = 0;
		/*
		uint64 reminder = num_rec % parts;
		for (uint32 i = 0; i < parts; ++i) {
		N2 = N1 + delta;
		if (i == 0)
		N2 += reminder;
		range_queue.push_back(std::make_tuple(N1, N2, i));
		N1 = N2;
		}*/

		uint64 smallest_fraction = 8;
		uint64 start_size = delta / smallest_fraction;
		uint64 step = (2 * smallest_fraction - 2) * delta / smallest_fraction / parts;
		uint64 cur_delta = start_size;
		for (uint32 i = 0; i < parts; ++i)
		{
			N2 = N1 + cur_delta;
			cur_delta += step;
			if (i == parts - 1)
				N2 = num_rec;
			range_queue.push_back(std::make_tuple(N1, N2, parts - i - 1));
			N1 = N2;
		}
		reverse(range_queue.begin(), range_queue.end());

		cur_idx = 0;
		if (parts)
			done = false;
	}


	bool get(uint64 &n1, uint64 &n2, uint32 &part_id)
	{
		std::lock_guard<std::mutex> lg(m);

		if (!done)
		{
			//n1 = std::get<0>(range_queue.back());
			//n2 = std::get<1>(range_queue.back());
			//range_queue.pop_back();

			n1 = std::get<0>(range_queue[cur_idx]);
			n2 = std::get<1>(range_queue[cur_idx]);;
			part_id = std::get<2>(range_queue[cur_idx]);
			cur_idx++;
			if (cur_idx == range_queue.size())
				done = true;
			return true;
		}
		return false;
		//TODO: dodaæ obs³ugê badania poprawnych danych, tzn przynajmniej jednego rekordu
	}

	void reset_indices()
	{
		cur_idx = 0;
		done = false;  //TODO: dodaæ obs³ugê badania poprawnych danych
	}
};

template <typename KMER_T, typename COUNTER_TYPE>
void pierwsze_kolko_etap1(uint32_t /*th_id*/, KMER_T *kmers, uint64 /*n_recs*/, uint32_t /*n_threads*/,
	//	uint64_t per_thread, std::vector<ALIGN_ARRAY COUNTER_TYPE[256]> &histos,
	uint64_t /*per_thread*/, std::vector<std::array<COUNTER_TYPE, 256>> &histos,
	uint32 byte, CRangeQueue& rq)
	//(std::thread([th_id, kmers, n_recs, n_threads, per_thread, &histos, byte]
{
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*	CStopWatch radix_mid_timer;
	if (byte == 15)
	radix_mid_timer.startTimer();*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef MEASURE_TIMES
	CThreadWatch tw;
	tw.startTimer();
#endif

	//	ALIGN_ARRAY COUNTER_TYPE myHisto[256] = { 0 };
	COUNTER_TYPE myHisto[256] = { 0 };

	uint64 idx1, idx2;
	uint32 part_id;
	//	uint8_t* ptr;
	//	uint64_t n;

	while (rq.get(idx1, idx2, part_id))
	{
		memset(myHisto, 0, 256 * sizeof(COUNTER_TYPE));

		//--->>>>uint8_t* ptr = (uint8_t*)(kmers + th_id*per_thread) + byte;
		uint8_t* ptr = (uint8_t*)(kmers + idx1) + byte;
		//--->>>>uint64_t n = per_thread;
		uint64_t n = idx2 - idx1;
		//if (th_id == n_threads - 1)
		//	n = n_recs - th_id*per_thread;

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
			//		for (uint64 i = 0; i < n / 2; ++i)
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
			//histos[th_id][i] = myHisto[i];
			histos[part_id][i] = myHisto[i];
		}

	}
#ifdef MEASURE_TIMES
	tw.stopTimer();

	times_byte_total[byte] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
	times_satish_stages[byte][0] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
#endif
}


//----------------------------------------------------------------------
template <typename KMER_T, typename COUNTER_TYPE>
void pierwsze_kolko_etap2(uint32_t /*th_id*/, KMER_T *kmers, KMER_T* tmp,
	uint64 /*n_recs*/, uint32_t /*n_threads*/, uint64_t /*per_thread*/, uint32 byte,
	//	std::vector<ALIGN_ARRAY COUNTER_TYPE[256]> &histos,
	std::vector<std::array<COUNTER_TYPE, 256>> &histos,
	std::vector<uchar*> &_raw_buffers,
	//	std::vector<ALIGN_ARRAY COUNTER_TYPE[256]> &threads_histos,
	std::vector<std::array<COUNTER_TYPE, 256>> &threads_histos,
	CMemoryPool* pmm_radix_buf,
	CRangeQueue& rq)
	//std::thread([th_id, kmers, tmp, n_recs, n_threads, per_thread, byte, 
	//&histos, &_raw_buffers, &threads_histos, pmm_radix_buf]
{
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*	CStopWatch radix_mid_timer;
	if (byte == 15)
	radix_mid_timer.startTimer();*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef MEASURE_TIMES
	CThreadWatch tw;
	tw.startTimer();
#endif

	ALIGN_ARRAY COUNTER_TYPE myHisto[256];

	uint64 idx1, idx2;
	uint32 part_id;

	uint8_t* ptr;
	uint64_t n;
	KMER_T* src;
	constexpr uint32_t BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(KMER_T) / 8];
	constexpr uint32_t BUFFER_WIDTH_IN_128BIT_WORDS = BUFFER_WIDTH * sizeof(KMER_T) / 16;
	constexpr uint32_t BUFFER_16B_ALIGNED = sizeof(KMER_T) % 16 == 0;
	//constexpr uint32_t BUFFER_16B_ALIGNED = 1;		// At the 1st level the pointers are always aligned to 16B

	uchar* raw_buffer;
	uchar* buffer;
	KMER_T *Buffer;
	uint8_t byteValue = 0;
	int index_x = 0;

	while (rq.get(idx1, idx2, part_id))
	{
		for (int i = 0; i < 256; ++i)
			myHisto[i] = histos[part_id][i]; //myHisto[i] = histos[th_id][i];

											 //--->>>>uint8_t* ptr = (uint8_t*)(kmers + th_id*per_thread) + byte;
		ptr = (uint8_t*)(kmers + idx1) + byte;
		//--->>>>uint64_t n = per_thread;
		n = idx2 - idx1;

		pmm_radix_buf->reserve(raw_buffer);
		//raw_buffer = new uchar[69632];
		_raw_buffers[part_id] = raw_buffer;
		buffer = raw_buffer;
		while ((uint64_t)buffer % ALIGNMENT)
			++buffer;
		Buffer = (KMER_T*)buffer;



		//if (th_id == n_threads - 1)
		//n = n_recs - th_id*per_thread;

		src = kmers + idx1;   //kmers + th_id*per_thread;

		byteValue = 0;
		index_x = 0;
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
				//				memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
				IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);

			ptr += sizeof(KMER_T);
		case 2:
			byteValue = *ptr;
			index_x = myHisto[byteValue] % BUFFER_WIDTH;
			Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 2];
			myHisto[byteValue]++;
			if (index_x == (BUFFER_WIDTH - 1))
				//				memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
				IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);

			ptr += sizeof(KMER_T);
		case 1:
			byteValue = *ptr;
			index_x = myHisto[byteValue] % BUFFER_WIDTH;
			Buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 1];
			myHisto[byteValue]++;
			if (index_x == (BUFFER_WIDTH - 1))
				//				memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
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
				//				memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
				IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);

			ptr += sizeof(KMER_T);

			byteValue = *ptr;
			index_x = myHisto[byteValue] % BUFFER_WIDTH;
			Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 1];
			myHisto[byteValue]++;
			if (index_x == (BUFFER_WIDTH - 1))
				//				memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
				IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);

			ptr += sizeof(KMER_T);

			byteValue = *ptr;
			index_x = myHisto[byteValue] % BUFFER_WIDTH;
			Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 2];
			myHisto[byteValue]++;
			if (index_x == (BUFFER_WIDTH - 1))
				//				memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
				IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);

			ptr += sizeof(KMER_T);

			byteValue = *ptr;
			index_x = myHisto[byteValue] % BUFFER_WIDTH;
			Buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 3];
			myHisto[byteValue]++;
			if (index_x == (BUFFER_WIDTH - 1))
				//				memcpy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(KMER_T));
				IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH]);

			ptr += sizeof(KMER_T);
		}

		for (uint32 i = 0; i < 256; ++i)
			threads_histos[part_id][i] = myHisto[i];  //threads_histos[th_id][i] = myHisto[i];
	}

#ifdef MEASURE_TIMES
	tw.stopTimer();

	times_byte_total[byte] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
	times_satish_stages[byte][1] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
#endif
}
//-----------------------------------------------------------
template <typename KMER_T, typename COUNTER_TYPE>
void pierwsze_kolko_etap3(uint32_t /*th_id*/, KMER_T */*kmers*/, KMER_T* tmp,
	uint64 /*n_recs*/, uint32_t /*n_threads*/, uint64_t /*per_thread*/, uint32 /*byte*/,
	//	std::vector<ALIGN_ARRAY COUNTER_TYPE[256]> &histos,
	std::vector<std::array<COUNTER_TYPE, 256>> &histos,
	std::vector<uchar*> &_raw_buffers,
	//	std::vector<ALIGN_ARRAY COUNTER_TYPE[256]> &threads_histos,
	std::vector<std::array<COUNTER_TYPE, 256>> &threads_histos,
	CMemoryPool* pmm_radix_buf,
	CRangeQueue& rq)

	//(std::thread([th_id, kmers, tmp, n_recs, n_threads, per_thread, byte, &histos, &_raw_buffers, &threads_histos,
	//	pmm_radix_buf]
{
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*	CStopWatch radix_mid_timer;
	if (byte == 15)
	radix_mid_timer.startTimer();*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef MEASURE_TIMES
	CThreadWatch tw;
	tw.startTimer();
#endif


	uint64 idx1, idx2;
	uint32 part_id;
	ALIGN_ARRAY COUNTER_TYPE myHisto[256];
	uchar* raw_buffer;
	uchar* buffer;
	KMER_T *Buffer;

	int64_t elemInBuffer;
	int64_t index_stop;
	int64_t index_start;
	int64_t elemWrittenIntoBuffer;
	const uint32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(KMER_T) / 8];

	while (rq.get(idx1, idx2, part_id))
	{
		raw_buffer = _raw_buffers[part_id];
		buffer = raw_buffer;
		while ((uint64_t)buffer % ALIGNMENT)
			++buffer;
		Buffer = (KMER_T*)buffer;


		for (int i = 0; i < 256; ++i)
			myHisto[i] = threads_histos[part_id][i];
		for (uint32_t private_i = 0; private_i < 256; private_i++)
		{
			index_stop = myHisto[private_i] % BUFFER_WIDTH;
			index_start = histos[part_id][private_i] % BUFFER_WIDTH;
			elemWrittenIntoBuffer = myHisto[private_i] - histos[part_id][private_i];

			if ((index_stop - elemWrittenIntoBuffer) <= 0)
				elemInBuffer = index_stop;
			else
				elemInBuffer = index_stop - index_start;

			if (elemInBuffer != 0)
				//				memcpy(&tmp[myHisto[private_i] - elemInBuffer], &Buffer[private_i * BUFFER_WIDTH + (myHisto[private_i] - elemInBuffer) % BUFFER_WIDTH], (elemInBuffer)*sizeof(KMER_T));
				IntrCopy64fun(&tmp[myHisto[private_i] - elemInBuffer],
					&Buffer[private_i * BUFFER_WIDTH + (myHisto[private_i] - elemInBuffer) % BUFFER_WIDTH], elemInBuffer * sizeof(KMER_T) / 8);
		}
		pmm_radix_buf->free(raw_buffer);
	}
#ifdef MEASURE_TIMES
	tw.stopTimer();

	times_byte_total[byte] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
	times_satish_stages[byte][2] += (uint64_t)(tw.getElapsedTime() * 1000000000.0);
#endif
}
#endif

// ***** EOF