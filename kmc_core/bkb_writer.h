/*
    This file is a part of KMC software distributed under GNU GPL 3 licence.
    The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

    Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

    Version: 2.2.0
    Date   : 2015-04-15
*/

#ifndef _BKB_WRITER_H
#define _BKB_WRITER_H
#include "../kmc/definitions.h"
class CBigBinDesc;
class CBigBinSortedPartQueue;
class CCompletedBinsCollector;
class CDiskLogger;
class CMemoryPool;
struct CKMCParams;
struct CKMCQueues;

//************************************************************************************************************
// CBigKmerBinWriter - Write sub bins to  HDD
//************************************************************************************************************
class CBigKmerBinWriter {
	int32 bin_id, sub_bin_id;
	CBigBinSortedPartQueue* bbspq;
	CCompletedBinsCollector* sm_cbc;
	CDiskLogger* disk_logger;
	CMemoryPool * sm_pmm_sorter_suffixes;
	CMemoryPool * sm_pmm_sorter_lut;

	std::string working_directory;
	CBigBinDesc* bbd;
	std::string GetName();
  public:
	CBigKmerBinWriter(CKMCParams& Params, CKMCQueues& Queues);
	void Process();
};

//************************************************************************************************************
// CWBigKmerBinWriter - wrapper for multithreading purposes
//************************************************************************************************************
class CWBigKmerBinWriter {
	CBigKmerBinWriter* bkb_writer;
  public:
	CWBigKmerBinWriter(CKMCParams& Params, CKMCQueues& Queues);
	void operator()();
	~CWBigKmerBinWriter();

};

#endif


// ***** EOF
