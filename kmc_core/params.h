/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _PARAMS_H
#define _PARAMS_H

#include "defs.h"
#include "kmc_runner.h"
#include "queues.h"
#include "s_mapper.h"
#include <vector>
#include <string>

//typedef enum {fasta, fastq, multiline_fasta, bam} input_type;

using InputType = KMC::InputFileType;
using OutputType = KMC::OutputFileType;

using namespace std;

// Structure for passing KMC parameters
struct CKMCParams {
//TODO: reconsider !!!!
#ifdef DEVELOP_MODE
	bool p_verbose_log = false;         // verbose log
#endif
	// File names
	vector<string> input_file_names;
	string output_file_name;
	string working_directory;
	InputType file_type;
	OutputType output_type;

	string json_summary_file_name = "";
	bool without_output = false;

	uint32 lut_prefix_len;

	uint32 KMER_T_size;

	// Memory sizes
	int64 max_mem_size;				// maximum amount of memory to be used in GBs;	default: 30GB
	int64 max_mem_storer;			// maximum amount of memory for internal buffers of KmerStorer
	int64 max_mem_stage2;			// maximum amount of memory in stage 2
	int64 max_mem_storer_pkg;		// maximum amount of memory for single package

	int64 mem_tot_pmm_bins;			// maximal amount of memory per pool memory manager (PMM) of bin parts
	int64 mem_part_pmm_bins;		// maximal amount of memory per single part of memory maintained by PMM of bin parts
	int64 mem_tot_pmm_fastq;
	int64 mem_part_pmm_fastq;
	int64 mem_part_pmm_reads;
	int64 mem_tot_pmm_reads;
	int64 mem_part_pmm_radix_buf;
	int64 mem_tot_pmm_radix_buf;	
	int64 mem_part_pmm_cnts_sort;	
	int64 mem_tot_pmm_stats;
	int64 mem_part_pmm_stats;
	
	int64 mem_part_pmm_binary_file_reader;
	int64 mem_tot_pmm_binary_file_reader;
	int64 mem_part_small_k_buf;
	int64 mem_tot_small_k_buf;
	int64 mem_part_small_k_completer;
	int64 mem_tot_small_k_completer;

	KMC::ILogger* verboseLogger;
	KMC::IPercentProgressObserver* percentProgressObserver;
	KMC::ILogger* warningsLogger;
#ifdef DEVELOP_MODE
	bool verbose_log;
#endif

	int kmer_len;			// kmer length
	int signature_len;
	int cutoff_min;			// exclude k-mers occurring less than times
	int64 cutoff_max;			// exclude k-mers occurring more than times
	int64 counter_max;		// maximal counter value	
	bool use_strict_mem;	// use strict memory limit mode
	bool homopolymer_compressed; //count homopolymer compressed k-mers
	bool both_strands;		// find canonical representation of each k-mer
	bool mem_mode;			// use RAM instead of disk

	int n_bins;				// number of bins;
	int bin_part_size;		// size of a bin part; fixed: 2^15
	int fastq_buffer_size;	// size of FASTQ file buffer; fixed: 2^23

	int n_threads;			// number of cores
	int n_readers;			// number of FASTQ readers; default: 1
	int n_splitters;		// number of splitters; default: 1
	int n_sorters;			// number of sorters; default: 1
	//vector<int> n_sorting_threads;// number of OMP threads per sorters
	uint32 max_x;					//k+x-mers will be counted

	//params for strict memory mode
	int sm_n_uncompactors;
	int sm_n_sorting_threads;	
	int sm_n_mergers;

	int64 sm_mem_part_input_file;
	int64 sm_mem_tot_input_file;	
	int64 sm_mem_part_expand;	
	int64 sm_mem_tot_expand;
	int64 sm_mem_part_sort;
	int64 sm_mem_tot_sort;
	int64 sm_mem_part_suffixes;
	int64 sm_mem_tot_suffixes;
	int64 sm_mem_part_lut;
	int64 sm_mem_tot_lut;

	int64 sm_mem_part_sub_bin_lut;
	int64 sm_mem_tot_sub_bin_lut;
	int64 sm_mem_part_sub_bin_suff;
	int64 sm_mem_tot_sub_bin_suff;
	int64 sm_mem_part_merger_lut;
	int64 sm_mem_tot_merger_lut;
	int64 sm_mem_part_merger_suff;
	int64 sm_mem_tot_merger_suff;
};

// Structure for passing KMC queues and monitors to threads
struct CKMCQueues 
{
	//Signature mapper
	CSignatureMapper* s_mapper;
	// Memory monitors
	CMemoryMonitor *mm;

	vector<CBinaryPackQueue*> binary_pack_queues;

	CBamTaskManager* bam_task_manager = nullptr;

	// Queues
	CInputFilesQueue *input_files_queue;
	CPartQueue *part_queue;
	CStatsPartQueue* stats_part_queue;

	CBinPartQueue *bpq;
	CBinDesc *bd;
	CExpanderPackDesc* epd;
	CBinQueue *bq;
	CKmerQueue *kq;
	CMemoryPool *pmm_bins, *pmm_reads, *pmm_radix_buf, *pmm_stats, *pmm_binary_file_reader;
	CMemoryPoolWithBamSupport *pmm_fastq;
	CMissingEOL_at_EOF_counter* missingEOL_at_EOF_counter{};
	CMemoryBins *memory_bins;
	CMemoryPool* pmm_small_k_buf, *pmm_small_k_completer;


	CDiskLogger* disk_logger;

	//for strict memory mode
	CTooLargeBinsQueue* tlbq;
	CBigBinPartQueue* bbpq;
	CBigBinKXmersQueue* bbkq;
	CBigBinDesc* bbd;
	CBigBinKmerPartQueue* bbkpq;
	CBigBinSortedPartQueue* bbspq;
	CKMCQueues() {}
	CMemoryPool* sm_pmm_input_file, *sm_pmm_expand, *sm_pmm_sort, *sm_pmm_sorter_suffixes, *sm_pmm_sorter_lut, *sm_pmm_sub_bin_lut, *sm_pmm_sub_bin_suff, *sm_pmm_merger_lut, *sm_pmm_merger_suff;
	
	CCompletedBinsCollector* sm_cbc;
	CSortersManager* sorters_manager = nullptr;
};

#endif

// ***** EOF
