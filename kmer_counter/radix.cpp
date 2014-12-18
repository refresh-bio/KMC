#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.1
  Date   : 2014-12-18
*/

#include <stdio.h>
#include "radix.h"

//----------------------------------------------------------------------------------
/*Parallel radix sort. The input data to be sorted are divided evenly among threads. 
  Each thread is responsible for building a local histogram to enable sorting keys 
  according to a given digit. Then a global histogram is created as a combination
  of local ones and the write offset (location) to which each digit should be written
 is computed. Finally, threads scatter the data to the appropriate locations.*/
template<typename COUNTER_TYPE>
void RadixOMP_uint8(uint32 *SourcePtr, uint32 *DestPtr, const int64 SourceSize, unsigned rec_size, unsigned data_offset, unsigned data_size, const unsigned n_phases, const unsigned n_threads)
{
/* SourceSize - number of records */
/* rec_size - in bytes */
/* data_offset - in bytes*/
/* data_size - in bytes - not used now */

    
#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) COUNTER_TYPE ByteCounter[MAX_NUM_THREADS][256];
#else
	COUNTER_TYPE ByteCounter[MAX_NUM_THREADS][256] __attribute__((aligned(ALIGNMENT)));
#endif

#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) COUNTER_TYPE globalHisto[256];
#else
	COUNTER_TYPE globalHisto[256] __attribute__((aligned(ALIGNMENT)));
#endif

#pragma omp parallel num_threads(n_threads)
	{ 
		int myID = omp_get_thread_num();
		uint8_t ByteIndex = 0;
		long long i;					
		COUNTER_TYPE prevSum;
		COUNTER_TYPE temp;
		uint32 n;
		
		int private_i;
		int byteValue;
	
		int64 SourceSize_in_bytes = SourceSize * rec_size;

		uint8_t *char_ptr_tempSource = (uint8_t*)(SourcePtr);
		uint8_t *char_ptr_tempDest = (uint8_t*)(DestPtr);
		uint8_t *char_tempPtr;

#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) COUNTER_TYPE privateByteCounter[256] = {0};
#else
	__attribute__((aligned(ALIGNMENT)))  COUNTER_TYPE privateByteCounter[256] = {0};
#endif

		for(uint32 privatePhaseCounter = 0; privatePhaseCounter < n_phases; privatePhaseCounter++)
		{
			#pragma omp for private(i) schedule(static) 
			for(i = data_offset; i < SourceSize_in_bytes; i = i + rec_size)
			{
				byteValue = *(&char_ptr_tempSource[i] + ByteIndex);
				
				++privateByteCounter[byteValue];
			}	
			A_memcpy(&ByteCounter[myID][0], privateByteCounter, sizeof(privateByteCounter));

			#pragma omp barrier

			#pragma omp for schedule(static)
			for(i = 0; i < 256; ++i)
			{
				prevSum = 0; 
				for(n = 0; n < n_threads; n++)
				{
					temp = ByteCounter[n][i];
					ByteCounter[n][i] = prevSum;
					prevSum += temp; 
				}
				globalHisto[i] = prevSum;
			}	

			#pragma omp single
			{
				prevSum = 0; 
				for(i = 0; i < 256; ++i)
				{
					temp = globalHisto[i];
					globalHisto[i] = prevSum;
					prevSum += temp; 
				}	
			}


			for (private_i = 0; private_i < 256; private_i++)
				ByteCounter[myID][private_i] += globalHisto[private_i];

			A_memcpy(privateByteCounter, &ByteCounter[myID][0], sizeof(privateByteCounter));
		
			#pragma omp for schedule(static)
			for(i = data_offset; i < SourceSize_in_bytes; i = i + rec_size)
			{
				byteValue = *(&char_ptr_tempSource[i] + ByteIndex);

				memcpy(&char_ptr_tempDest[privateByteCounter[byteValue] * rec_size], &char_ptr_tempSource[i - data_offset], rec_size);

				(privateByteCounter[byteValue])++;
			}

		
			#pragma omp barrier

			char_tempPtr = char_ptr_tempDest;
			char_ptr_tempDest = char_ptr_tempSource;
			char_ptr_tempSource = char_tempPtr;
			ByteIndex++;
			memset(privateByteCounter, 0, sizeof(privateByteCounter));
		}
	}
}

//----------------------------------------------------------------------------------
void RadixSort_uint8(uint32 *&data_ptr, uint32 *&tmp_ptr, uint64 size, unsigned rec_size, unsigned data_offset, unsigned data_size, const unsigned n_phases, const unsigned n_threads)
{
	if(size * rec_size >= (1ull << 32))
		RadixOMP_uint8<uint64>(data_ptr, tmp_ptr, size, rec_size, data_offset, data_size, n_phases, n_threads);
	else
		RadixOMP_uint8<uint32>(data_ptr, tmp_ptr, size, rec_size, data_offset, data_size, n_phases, n_threads);
}


//----------------------------------------------------------------------------------
/*Parallel radix sort. Parallelization scheme taken from 
  Satish, N., Kim, C., Chhugani, J., Nguyen, A.D., Lee, V.W., Kim, D., Dubey, P. (2010). 
  Fast Sort on CPUs and GPUs. A Case for Bandwidth Oblivious SIMD Sort. 
  Proc. of the 2010 Int. Conf. on Management of data, pp. 351–362. 
  The usage of software-managed buffers in the writting phase results in diminishing
  the influence of irregular memory accesses. As the number of cache conflict misses
  is reduced better efficiency is reached.*/
template<typename COUNTER_TYPE, typename INT_TYPE>
	void RadixOMP_buffer(CMemoryPool *pmm_radix_buf, uint64 *Source, uint64 *Dest, const int64 SourceSize, const unsigned n_phases, const unsigned n_threads)
{
#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) COUNTER_TYPE ByteCounter[MAX_NUM_THREADS][256];
#else
	COUNTER_TYPE ByteCounter[MAX_NUM_THREADS][256] __attribute__((aligned(ALIGNMENT)));
#endif

#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) COUNTER_TYPE globalHisto[256];
#else
	COUNTER_TYPE globalHisto[256] __attribute__((aligned(ALIGNMENT)));
#endif
	
#pragma omp parallel num_threads(n_threads)
	{ 
		int myID = omp_get_thread_num();
		uint8_t ByteIndex = 0;
		long long i;
		COUNTER_TYPE prevSum;
		COUNTER_TYPE temp;

		uint32 n;
		
		int index_x;
		int private_i;
		int byteValue;
		uint64 *tempSource = Source;
		uint64 *tempDest = Dest;
		uint64 *tempPtr;
					
		uint64 *raw_Buffer;
		pmm_radix_buf->reserve(raw_Buffer);
		uint64 *Buffer = raw_Buffer;     
		
		while(((unsigned long long) Buffer) % ALIGNMENT)
			Buffer++; 

#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) COUNTER_TYPE privateByteCounter[256] = {0};
#else
	__attribute__((aligned(ALIGNMENT)))  COUNTER_TYPE privateByteCounter[256] = {0};
#endif
		
   		for(uint32 privatePhaseCounter = 0; privatePhaseCounter < n_phases; privatePhaseCounter++)
		{
			#pragma omp for private(i) schedule(static) 
			for(i = 0; i < SourceSize; ++i)
			{
				byteValue = *(reinterpret_cast<const uint8_t*>(&tempSource[i]) + ByteIndex);
				++privateByteCounter[byteValue];
			}	
			A_memcpy(&ByteCounter[myID][0], privateByteCounter, sizeof(privateByteCounter));

			#pragma omp barrier

			#pragma omp for schedule(static)
			for(i = 0; i < 256; ++i)
			{
				prevSum = 0; 
				for(n = 0; n < n_threads; n++)
				{
					temp = ByteCounter[n][i];
					ByteCounter[n][i] = prevSum;
					prevSum += temp; 
				}
				globalHisto[i] = prevSum;
			}	

			#pragma omp single
			{
				prevSum = 0; 
				for(i = 0; i < 256; ++i)
				{
					temp = globalHisto[i];
					globalHisto[i] = prevSum;
					prevSum += temp; 
				}	
			}

			for (private_i = 0; private_i < 256; private_i++)
				ByteCounter[myID][private_i] += globalHisto[private_i];

			A_memcpy(privateByteCounter, &ByteCounter[myID][0], sizeof(privateByteCounter));
		

			#pragma omp for schedule(static)
			for(i = 0; i < SourceSize; ++i)
			{
				byteValue = *(reinterpret_cast<const uint8_t*>(&tempSource[i]) + ByteIndex);

				index_x = privateByteCounter[byteValue] % BUFFER_WIDTH;

				Buffer[byteValue * BUFFER_WIDTH + index_x] = tempSource[i];
			
				privateByteCounter[byteValue]++;

				if(index_x == (BUFFER_WIDTH -1))
					A_memcpy ( &tempDest[privateByteCounter[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(uint64) );
			} //end_for

			INT_TYPE elemInBuffer;
			INT_TYPE index_stop;
			INT_TYPE index_start;
			INT_TYPE elemWrittenIntoBuffer;
		
			for(private_i = 0; private_i < 256; private_i++)
			{
				index_stop = privateByteCounter[private_i] % BUFFER_WIDTH;
				index_start = ByteCounter[myID][private_i] % BUFFER_WIDTH;
				elemWrittenIntoBuffer = privateByteCounter[private_i] - ByteCounter[myID][private_i];

				if((index_stop - elemWrittenIntoBuffer) <= 0)
					elemInBuffer = index_stop;
				else
					elemInBuffer = index_stop - index_start;

				if(elemInBuffer != 0)
					A_memcpy ( &tempDest[privateByteCounter[private_i] - elemInBuffer], &Buffer[private_i * BUFFER_WIDTH + (privateByteCounter[private_i] - elemInBuffer)%BUFFER_WIDTH], (elemInBuffer)*sizeof(uint64) );
			
			}
			#pragma omp barrier

			tempPtr = tempDest;
			tempDest = tempSource;
			tempSource = tempPtr;
			ByteIndex++;
			memset(privateByteCounter, 0, sizeof(privateByteCounter));
		}
		pmm_radix_buf->free(raw_Buffer);
	}
}

//----------------------------------------------------------------------------------
void RadixSort_buffer(CMemoryPool *pmm_radix_buf, uint64 *&data, uint64 *&tmp, uint64 size, const unsigned n_phases, const unsigned n_threads)
{
	if(size >= (1ull << 31))
		RadixOMP_buffer<uint64, int64>(pmm_radix_buf, data, tmp, size, n_phases, n_threads);
	else
		RadixOMP_buffer<uint32, int32>(pmm_radix_buf, data, tmp, size, n_phases, n_threads);
}

// ***** EOF
