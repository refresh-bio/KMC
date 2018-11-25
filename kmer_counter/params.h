/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _PARAMS_H
#define _PARAMS_H

#include "defs.h"
#include "queues.h"
#include "s_mapper.h"
#include <vector>
#include <string>

typedef enum {fasta, fastq, multiline_fasta, bam} input_type;


using namespace std;

// Structure for passing KMC parameters
struct CKMCParams {
	
	// Input parameters
	int p_m;							// max. total RAM usage
	int p_k;							// k-mer length
	int p_t;							// no. of threads
	int p_sf;							// no. of reading threads
	int p_sp;							// no. of splitting threads
	int p_sr;							// no. of threads for 2nd stage
	int p_ci;							// do not count k-mers occurring less than
	int64 p_cx;							// do not count k-mers occurring more than
	int64 p_cs;							// maximal counter value	
	bool p_strict_mem;					// use strict memory limit mode
	bool p_mem_mode;					// use RAM instead of disk	
	input_type p_file_type;				// input in FASTA format
	bool p_verbose;						// verbose mode
	bool p_without_output = false;		// do not create output files 
#ifdef DEVELOP_MODE
	bool p_verbose_log = false;         // verbose log
#endif
	bool p_both_strands;				// compute canonical k-mer representation
	int p_p1;							// signature length	
	int p_n_bins;						// no. of bins
	int p_smso;							// no. of OpenMP threads for sorting in strict memory mode
	int p_smun;							// no. of uncompacting threads in strict memory mode
	int p_smme;							// no. of merging threads in strict memory mode

	// File names
	vector<string> input_file_names;
	string output_file_name;
	string working_directory;
	input_type file_type;

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

	bool verbose;	
#ifdef DEVELOP_MODE
	bool verbose_log;
#endif

	int kmer_len;			// kmer length
	int signature_len;
	int cutoff_min;			// exclude k-mers occurring less than times
	int64 cutoff_max;			// exclude k-mers occurring more than times
	int64 counter_max;		// maximal counter value	
	bool use_strict_mem;	// use strict memory limit mode
	bool both_strands;		// find canonical representation of each k-mer
	bool mem_mode;			// use RAM instead of disk

	int n_bins;				// number of bins;
	int bin_part_size;		// size of a bin part; fixed: 2^15
	int fastq_buffer_size;	// size of FASTQ file buffer; fixed: 2^23

	int n_threads;			// number of cores
	int n_readers;			// number of FASTQ readers; default: 1
	int n_splitters;		// number of splitters; default: 1
	int n_sorters;			// number of sorters; default: 1
	vector<int> n_sorting_threads;// number of OMP threads per sorters
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

	CKMCParams()
	{
		p_m = 12;
		p_k = 25;
		p_t = 0;
		p_sf = 0;
		p_sp = 0;
		p_sr = 0;
		p_smme = p_smso = p_smun = 0;
		p_ci = 2;
		p_cx = 1000000000;
		p_cs = 255;		
		p_strict_mem = false;
		p_mem_mode = false;		
		p_file_type = fastq;
		p_verbose = false;
		p_both_strands = true;
		p_p1 = 9;	
		p_n_bins = 512;		
	}
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
