/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/

#ifndef _KB_STORER_H
#define _KB_STORER_H

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include "radix.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <tuple>
#include <stdio.h>

using namespace std;

//************************************************************************************************************
// CKmerBinStorer - storer of bins of k-mers
//************************************************************************************************************
class CKmerBinStorer
{
	uint64 total_size; 
	CMemoryPool *pmm_bins;
	string working_directory;
	int n_bins;
	CBinPartQueue *q_part;
	CBinDesc *bd;
	CExpanderPackDesc *epd;
	uint64 buffer_size_bytes;
	uint64 max_mem_buffer;
	uint64 max_mem_single_package;

	CDiskLogger *disk_logger;
	std::unique_ptr<uchar[]> tmp_buff;

	std::vector<uint64> buf_sizes;
	uint64 max_buf_size;
	uint32 max_buf_size_id;
	CTmpFilesOwner* tmp_files_owner;

	typedef list<tuple<uchar *, uint32, uint32>> elem_t; 
	std::vector<std::unique_ptr<elem_t>> buffer;

	void Release();
	string GetName(int n);
	void CheckBuffer();
	void ReleaseBuffer();
	void PutBinToTmpFile(uint32 n);
	
public:
	void GetTotal(uint64& _total)
	{
		_total = total_size;
	}
	CKmerBinStorer(CKMCParams &Params, CKMCQueues &Queues);
	~CKmerBinStorer();

	bool OpenFiles();
	void ProcessQueue();
};

//************************************************************************************************************
// CWKmerBinStorer - wrapper for multithreading purposes
//************************************************************************************************************
class CWKmerBinStorer {
	std::unique_ptr<CKmerBinStorer> kbs;

public:
	void GetTotal(uint64& _total)
	{
		kbs->GetTotal(_total);
	}
	CWKmerBinStorer(CKMCParams &Params, CKMCQueues &Queues);

	void operator()();
};

#endif

// ***** EOF
