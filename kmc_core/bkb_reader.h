
/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.1
Date   : 2019-05-19
*/

#ifndef _BKB_READER_H_
#define  _BKB_READER_H_

#include "params.h"

//************************************************************************************************************
// CBigKmerBinReader - reader of bins from distribution phase. Only in strict memory mode
//************************************************************************************************************

class CBigKmerBinReader
{
	CTooLargeBinsQueue * tlbq;
	CDiskLogger* disk_logger;
	CBinDesc* bd;
	CBigBinPartQueue* bbpq;
	CMemoryPool* sm_pmm_input_file;

	uint64 sm_mem_part_input_file;
public:
	CBigKmerBinReader(CKMCParams& Params, CKMCQueues& Queues);
	~CBigKmerBinReader();
	void ProcessBigBin();
};


//************************************************************************************************************
// CWBigKmerBinReader - wrapper for multithreading purposes
//************************************************************************************************************
class CWBigKmerBinReader
{
	CBigKmerBinReader* bkb_reader;
public:
	CWBigKmerBinReader(CKMCParams& Params, CKMCQueues& Queues);
	~CWBigKmerBinReader();
	void operator()();
};

#endif 

// ***** EOF