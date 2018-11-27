/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _SPLITTER_H
#define _SPLITTER_H

#include "defs.h"
#include "kmer.h"
#include "kb_collector.h"
#include "queues.h"
#include "s_mapper.h"
#include "mmer.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "small_k_buf.h"
#include "bam_utils.h"

using namespace std;

//************************************************************************************************************
// CSplitter class - splits kmers into bins according to their signatures
//************************************************************************************************************
class CSplitter {
	CMemoryMonitor *mm;
	uint64 total_kmers = 0;	
	uchar *part;
	uint64_t part_size, part_pos;
	CKmerBinCollector **bins;
	CBinPartQueue *bin_part_queue;
	CBinDesc *bd;
	CMemoryPool *pmm_reads;
	int64 mem_part_pmm_bins;
	int64 mem_part_pmm_reads;

	char codes[256];
	input_type file_type;
	bool both_strands;

	uint32 kmer_len;
	//uint32 prefix_len;
	uint32 signature_len;
	uint32 n_bins;	
	uint64 n_reads;//for multifasta its a sequences counter	

	CSignatureMapper* s_mapper;

	bool GetSeq(char *seq, uint32 &seq_size);

public:
	static uint32 MAX_LINE_SIZE;

	
	CSplitter(CKMCParams &Params, CKMCQueues &Queues); 
	void InitBins(CKMCParams &Params, CKMCQueues &Queues);	
	void CalcStats(uchar* _part, uint64 _part_size, uint32* _stats);
	bool ProcessReads(uchar *_part, uint64 _part_size);
	template<typename COUNTER_TYPE> bool ProcessReadsSmallK(uchar *_part, uint64 _part_size, CSmallKBuf<COUNTER_TYPE>& small_k_buf);	
	void Complete();
	inline void GetTotal(uint64 &_n_reads);
	inline uint64 GetTotalKmers();
	~CSplitter();
};

//----------------------------------------------------------------------------------
// Return the number of reads processed by splitter
void CSplitter::GetTotal(uint64 &_n_reads)
{
	_n_reads = n_reads;
}

//----------------------------------------------------------------------------------
// Return the number of kmers processed by splitter (!!! only for small k optimization)
uint64 CSplitter::GetTotalKmers()
{
	return total_kmers;
}

//************************************************************************************************************
// CWSplitter class - wrapper for multithreading purposes
//************************************************************************************************************

//----------------------------------------------------------------------------------
class CWSplitter {
	CPartQueue *pq;
	CBinPartQueue *bpq;
	CMemoryPool *pmm_fastq;

	CSplitter *spl;
	uint64 n_reads;

public:
	CWSplitter(CKMCParams &Params, CKMCQueues &Queues);	
	void operator()();
	void GetTotal(uint64 &_n_reads);
	~CWSplitter();
};

//************************************************************************************************************
// CWStatsSplitter class - wrapper for multithreading purposes
//************************************************************************************************************

//----------------------------------------------------------------------------------
class CWStatsSplitter {
	CStatsPartQueue *spq;
	CMemoryPool *pmm_fastq, *pmm_stats;
	uint32 *stats;
	CSplitter *spl;
	uint32 signature_len;

public:
	CWStatsSplitter(CKMCParams &Params, CKMCQueues &Queues);
	~CWStatsSplitter();

	void operator()();
	void GetStats(uint32* _stats);
};


//************************************************************************************************************
// CWSmallKSplitter class - wrapper for multithreading purposes
//************************************************************************************************************
//----------------------------------------------------------------------------------
template <typename COUNTER_TYPE> class CWSmallKSplitter {
	CPartQueue *pq;	
	CMemoryPool *pmm_fastq, *pmm_small_k;	
	CSmallKBuf<COUNTER_TYPE> small_k_buf;

	CSplitter *spl;
	uint64 n_reads;
	uint64 total_kmers;
	uint32 kmer_len;

public:
	CWSmallKSplitter(CKMCParams &Params, CKMCQueues &Queues);
	~CWSmallKSplitter();

	void operator()();
	void GetTotal(uint64 &_n_reads);

	CSmallKBuf<COUNTER_TYPE> GetResult()
	{
		return small_k_buf;
	}

	uint64 GetTotalKmers()
	{
		if (spl)
			return spl->GetTotalKmers();
		return total_kmers;
	}

	void Release()
	{
		pmm_small_k->free(small_k_buf.buf);
	}
};


#endif

// ***** EOF
