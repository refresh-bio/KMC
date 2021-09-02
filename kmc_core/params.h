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
#include <memory>
#include "libs/ntHash/ntHashWrapper.h"
#include "tmp_files_owner.h"

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
	uint32 max_x;			//k+x-mers will be counted

	KMC::EstimateHistogramCfg estimateHistogramCfg = KMC::EstimateHistogramCfg::DONT_ESTIMATE;
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
	std::unique_ptr<CSignatureMapper> s_mapper;

	std::vector<std::unique_ptr<CBinaryPackQueue>> binary_pack_queues;

	std::unique_ptr<CBamTaskManager> bam_task_manager;

	// Queues
	std::unique_ptr<CInputFilesQueue> input_files_queue;
	std::unique_ptr<CPartQueue> part_queue;
	std::unique_ptr<CStatsPartQueue> stats_part_queue;

	std::unique_ptr<CBinPartQueue> bpq;
	std::unique_ptr<CBinDesc> bd;
	std::unique_ptr<CExpanderPackDesc> epd;
	std::unique_ptr<CBinQueue> bq;
	std::unique_ptr<CKmerQueue> kq;
	std::unique_ptr<CMemoryPool> pmm_bins;
	std::unique_ptr<CMemoryPool> pmm_reads;
	std::unique_ptr<CMemoryPool> pmm_radix_buf;
	std::unique_ptr<CMemoryPool> pmm_stats;
	std::unique_ptr<CMemoryPool> pmm_binary_file_reader;

	std::unique_ptr<CMemoryPoolWithBamSupport> pmm_fastq;
	std::unique_ptr<CMissingEOL_at_EOF_counter> missingEOL_at_EOF_counter{};
	std::unique_ptr<CMemoryBins> memory_bins;
	std::unique_ptr<CMemoryPool> pmm_small_k_buf;
	std::unique_ptr<CMemoryPool> pmm_small_k_completer;


	std::unique_ptr<CDiskLogger> disk_logger;

	//for strict memory mode
	std::unique_ptr<CTooLargeBinsQueue> tlbq;
	std::unique_ptr<CBigBinPartQueue> bbpq;
	std::unique_ptr<CBigBinKXmersQueue> bbkq;
	std::unique_ptr<CBigBinDesc> bbd;
	std::unique_ptr<CBigBinKmerPartQueue> bbkpq;
	std::unique_ptr<CBigBinSortedPartQueue> bbspq;
	CKMCQueues() {}

	std::unique_ptr<CMemoryPool> sm_pmm_input_file;
	std::unique_ptr<CMemoryPool> sm_pmm_expand;
	std::unique_ptr<CMemoryPool> sm_pmm_sort;
	std::unique_ptr<CMemoryPool> sm_pmm_sorter_suffixes;
	std::unique_ptr<CMemoryPool> sm_pmm_sorter_lut;
	std::unique_ptr<CMemoryPool> sm_pmm_sub_bin_lut;
	std::unique_ptr<CMemoryPool> sm_pmm_sub_bin_suff;
	std::unique_ptr<CMemoryPool> sm_pmm_merger_lut;
	std::unique_ptr<CMemoryPool> sm_pmm_merger_suff;
	
	std::unique_ptr<CCompletedBinsCollector> sm_cbc;
	std::unique_ptr<CSortersManager> sorters_manager;

	std::unique_ptr<CntHashEstimator> ntHashEstimator;

	std::unique_ptr<CTmpFilesOwner> tmp_files_owner;
};

#endif

// ***** EOF
