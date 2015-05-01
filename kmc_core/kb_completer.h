/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.2.0
  Date   : 2015-04-15
*/
#ifndef _KB_COMPLETER_H
#define _KB_COMPLETER_H

#include "../kmc/definitions.h"
#include "params.h"
#include "kmer.h"
#include "radix.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <cstdio>


//************************************************************************************************************
// CKmerBinCompleter - complete the sorted bins and store in a file
//************************************************************************************************************
class CKmerBinCompleter {
	CMemoryMonitor *mm;
	std::string file_name, kmer_file_name, lut_file_name;
	CKmerQueue *kq;
	CBinDesc *bd;
	CSignatureMapper *s_mapper;
	uint32 *sig_map;
	uint64 _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total;
	uint64 n_recs;

	FILE *out_kmer, *out_lut;
	uint32 lut_pos;
	uint32 sig_map_size;
	uint64 counter_size;

	CMemoryBins *memory_bins;

	bool use_strict_mem;
	CBigBinKmerPartQueue* bbkpq;
	CMemoryPool *sm_pmm_merger_lut, *sm_pmm_merger_suff;

	uint32 lut_prefix_len;
	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total;
	uint32 kmer_t_size;
	int32 cutoff_min, cutoff_max;
	int32 counter_max;
	int32 kmer_len;
	int32 signature_len;
	bool use_quake;

	bool store_uint(FILE *out, uint64 x, uint32 size);	

public:
	CKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues);
	~CKmerBinCompleter();

	void ProcessBinsFirstStage();
	void ProcessBinsSecondStage();
	void GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total);
	void InitStage2(CKMCParams& Params, CKMCQueues& Queues);
};


//************************************************************************************************************
// CWKmerBinCompleter - wrapper for multithreading purposes
//************************************************************************************************************
class CWKmerBinCompleter {
	CKmerBinCompleter *kbc;

public:
	CWKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues);
	~CWKmerBinCompleter();

	void operator()(bool first_stage);

	void GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total);
	void InitStage2(CKMCParams& Params, CKMCQueues& Queues);
};

#endif

// ***** EOF
