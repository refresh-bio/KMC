/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _BKB_WRITER_H
#define _BKB_WRITER_H

#include "params.h"

//************************************************************************************************************
// CBigKmerBinWriter - Write sub bins to  HDD
//************************************************************************************************************
class CBigKmerBinWriter
{
	int32 bin_id, sub_bin_id;
	CBigBinSortedPartQueue* bbspq;
	CCompletedBinsCollector* sm_cbc;
	CDiskLogger* disk_logger;
	CMemoryPool * sm_pmm_sorter_suffixes;
	CMemoryPool * sm_pmm_sorter_lut;
	
	string working_directory;
	CBigBinDesc* bbd;
	string GetName();
public:
	CBigKmerBinWriter(CKMCParams& Params, CKMCQueues& Queues);
	void Process();	
};

//************************************************************************************************************
// CWBigKmerBinWriter - wrapper for multithreading purposes
//************************************************************************************************************
class CWBigKmerBinWriter
{
	CBigKmerBinWriter* bkb_writer;
public:
	CWBigKmerBinWriter(CKMCParams& Params, CKMCQueues& Queues);
	void operator()();
	~CWBigKmerBinWriter();
	
};

#endif


// ***** EOF 