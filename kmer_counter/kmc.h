/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _KMC_H
#define _KMC_H

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>
#include "queues.h"
#include "timer.h"
#include "fastq_reader.h"
#include "kb_collector.h"
#include "kb_completer.h"
#include "kb_reader.h"
#include "kb_sorter.h"
#include "kb_storer.h"
#include "s_mapper.h"
#include "splitter.h"
#include "cpu_info.h"
#include "small_sort.h"

#ifdef DEVELOP_MODE
#include "develop.h"
#endif
#include "bkb_reader.h"
#include "bkb_uncompactor.h"
#include "bkb_sorter.h"
#include "bkb_merger.h"
#include "bkb_writer.h"
#include "binary_reader.h"
#include "kmc_runner.h"

using namespace std;

template <unsigned SIZE> class CKMC {
	bool initialized;

	CStopWatch w0, heuristic_time , w1, w2, w3;//w3 - strict memory time

	// Parameters (input and internal)
	CKMCParams Params;

	// Memory monitor and queues
	CKMCQueues Queues;

	// Thread groups
	vector<thread> gr0_1, gr0_2;
	vector<thread> gr1_1, gr1_2, gr1_3, gr1_4, gr1_5, gr1_6;		// thread groups for 1st stage
	vector<thread> gr2_1, gr2_2, gr2_3;						// thread groups for 2nd stage

	uint64 n_reads;
	bool was_small_k_opt;

	// Threads
	vector<CWStatsFastqReader*> w_stats_fastqs;
	vector<CWStatsSplitter*> w_stats_splitters;
	vector<CWFastqReader*> w_fastqs;
	vector<CWSplitter*> w_splitters;

	CWBinaryFilesReader* w_bin_file_reader;

	CWKmerBinStorer *w_storer;

	CWKmerBinReader<SIZE>* w_reader;
	vector<CWKmerBinSorter<SIZE>*> w_sorters;
	CWKmerBinCompleter *w_completer;

	void SetThreads1Stage(const KMC::Stage1Params& stage1Params);
	void SetThreads2Stage(const KMC::Stage2Params& stage2Params);
	void SetThreadsStrictMemoryMode(const KMC::Stage2Params& stage2Params);

	void AdjustMemoryLimitsStrictMemoryMode();
	bool AdjustMemoryLimits();
	void AdjustMemoryLimitsStage2();

	void ShowSettingsStage1();
	void ShowSettingsStage2();
	void ShowSettingsSmallKOpt();


	//small k
	bool AdjustMemoryLimitsSmallK();
	vector<CWSmallKSplitter<uint64_t>*> w_small_k_splitters;
	KMC::Stage1Results ProcessSmallKOptimization_Stage1();
	KMC::Stage2Results ProcessSmallKOptimization_Stage2();


	void CheckAndReportMissingEOLs();
	
public:
	CKMC();
	~CKMC();
	
	void SetParamsStage1(const KMC::Stage1Params& stage1Params);
	void SetParamsStage2(const KMC::Stage2Params& stage2Params);
	KMC::Stage1Results ProcessStage1();
	KMC::Stage2Results ProcessStage2();
};

//----------------------------------------------------------------------------------
template <unsigned SIZE> CKMC<SIZE>::CKMC()
{
	initialized   = false;
	Params.kmer_len      = 0;
	Params.n_readers     = 1;
	Params.n_splitters   = 1;
	Params.n_sorters     = 1;	
	Queues.s_mapper = nullptr;
}

//----------------------------------------------------------------------------------
template <unsigned SIZE> CKMC<SIZE>::~CKMC()
{
}

//----------------------------------------------------------------------------------
// Set params of the first stage of k-mer counter
template <unsigned SIZE> void CKMC<SIZE>::SetParamsStage1(const KMC::Stage1Params& stage1Params)
{
	Params.input_file_names = stage1Params.GetInputFiles();
	Params.working_directory = stage1Params.GetTmpPath();
	Params.kmer_len = stage1Params.GetKmerLen();
	Params.file_type = stage1Params.GetInputFileType();
	Params.n_bins = stage1Params.GetNBins();

	if (Params.kmer_len % 32 == 0)
		Params.max_x = 0;
	else
		Params.max_x = MIN(31 - (Params.kmer_len % 32), KMER_X);

	Params.verboseLogger = stage1Params.GetVerboseLogger();
	Params.percentProgressObserver = stage1Params.GetPercentProgressObserver();
	Params.warningsLogger = stage1Params.GetWarningsLogger();

	// Technical parameters related to temporary files
	Params.signature_len = stage1Params.GetSignatureLen();
	Params.bin_part_size = 1 << 16;

	

//TODO: reconsider !!!!
//#ifdef DEVELOP_MODE
//	Params.verbose_log = Params.p_verbose_log;
//#endif

	Params.both_strands = stage1Params.GetCanonicalKmers();
	Params.homopolymer_compressed = stage1Params.GetHomopolymerCompressed();
	Params.mem_mode = stage1Params.GetRamOnlyMode();

	if (stage1Params.GetNReaders() && stage1Params.GetNSplitters())
	{
		Params.n_readers = NORM(stage1Params.GetNReaders(), 1, 32);
		Params.n_splitters = NORM(stage1Params.GetNSplitters(), 1, 32);
	}
	else
	{
		Params.n_threads = stage1Params.GetNThreads();
		if (stage1Params.GetNThreads() > stage1Params.GetMaxRamGB() * 64)
		{
			Params.n_threads = stage1Params.GetMaxRamGB() * 64;
			std::ostringstream ostr;
			ostr << "number of threads is reduced to " << Params.n_threads << " (maximum numer of threads equals 64 * MaxRamGB)";
			Params.warningsLogger->Log(ostr.str());
		}
		SetThreads1Stage(stage1Params);
	}

	Params.max_mem_size = NORM(((uint64)stage1Params.GetMaxRamGB()) * 1000000000ull, (uint64)MIN_MEM * 1000000000ull, 1024ull * 1000000000ull);
	Params.KMER_T_size = sizeof(CKmer<SIZE>);

	initialized = true;
}

//----------------------------------------------------------------------------------
// Set params of the second stage of k-mer counter
template <unsigned SIZE> void CKMC<SIZE>::SetParamsStage2(const KMC::Stage2Params& stage2Params)
{
	Params.output_type = stage2Params.GetOutputFileType();
	Params.output_file_name = stage2Params.GetOutputFileName();

	// Thresholds for counters
	Params.cutoff_min = stage2Params.GetCutoffMin();
	Params.cutoff_max = stage2Params.GetCutoffMax();
	Params.counter_max = stage2Params.GetCounterMax();
	
	if (Params.kmer_len > 9)
	{
		if ((uint64)Params.cutoff_max > ((1ull << 32) - 1))
		{
			std::ostringstream ostr;
			ostr << "for k > 9 maximum value of cutoff_max is 4294967295";
			Params.warningsLogger->Log(ostr.str());
			Params.cutoff_max = 4294967295;
		}
		if ((uint64)Params.counter_max > ((1ull << 32) - 1))
		{
			std::ostringstream ostr;
			ostr << "for k > 9 maximum value of counter_max is 4294967295";
			Params.warningsLogger->Log(ostr.str());
			Params.counter_max = 4294967295;
		}
	}

	if (Params.counter_max == 1)
	{
		std::ostringstream ostr;
		ostr << "using counter_max == 1 will cause not storying counters in KMC output file, all counters will be assumed to be 1. This is experimental and is not currently supported in kmc_tools. Will be implemented soon.";
		Params.warningsLogger->Log(ostr.str());
	}

	Params.without_output = stage2Params.GetWithoutOutput();
	Params.use_strict_mem = stage2Params.GetStrictMemoryMode();

	Params.max_mem_size = NORM(((uint64)stage2Params.GetMaxRamGB()) * 1000000000ull, (uint64)MIN_MEM * 1000000000ull, 1024ull * 1000000000ull);
	SetThreads2Stage(stage2Params);
	if(Params.use_strict_mem)
		SetThreadsStrictMemoryMode(stage2Params);
}

//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKMC<SIZE>::SetThreads1Stage(const KMC::Stage1Params& stage1Params)
{
	if (!stage1Params.GetNReaders() || !stage1Params.GetNSplitters())
	{
		int cores = Params.n_threads;
		bool gz_bz2 = false;
		vector<uint64> file_sizes;

		for (auto& p : Params.input_file_names)
		{
			if (p.size() > 3 && string(p.end() - 3, p.end()) == ".gz")
				gz_bz2 = true;
			else if (p.size() > 4 && string(p.end() - 4, p.end()) == ".bz2")
				gz_bz2 = true;

			FILE* tmp = my_fopen(p.c_str(), "rb");
			if (!tmp)
			{
				cerr << "Error: cannot open file: " << p.c_str();
				exit(1);
			}
			my_fseek(tmp, 0, SEEK_END);
			file_sizes.push_back(my_ftell(tmp));
			fclose(tmp);
		}
		if (gz_bz2)
		{
			sort(file_sizes.begin(), file_sizes.end(), greater<uint64>());
			uint64 file_size_threshold = (uint64)(file_sizes.front() * 0.05);
			int32 n_allowed_files = 0;
			for (auto& p : file_sizes)
				if (p > file_size_threshold)
					++n_allowed_files;
			Params.n_readers = MIN(n_allowed_files, MAX(1, cores / 2));
		}
		else if (Params.file_type == InputType::BAM)
		{
			Params.n_readers = MAX(1, Params.n_threads / 2); //TODO: for now, split distribute threads equally for bam input, it seems to work quite well but I didnt perform detailed tests
		}
		else
			Params.n_readers = 1;
		Params.n_splitters = MAX(1, cores - Params.n_readers);
	}
}

//----------------------------------------------------------------------------------
template<unsigned SIZE> void CKMC<SIZE>::SetThreads2Stage(const KMC::Stage2Params& stage2Params)
{
	Params.n_sorters = stage2Params.GetNThreads();
}
//----------------------------------------------------------------------------------
template<unsigned SIZE> void CKMC<SIZE>::SetThreadsStrictMemoryMode(const KMC::Stage2Params& stage2Params)
{
	Params.sm_n_mergers = stage2Params.GetStrictMemoryNMergers();
	Params.sm_n_uncompactors = stage2Params.GetStrictMemoryNUncompactors();
	Params.sm_n_sorting_threads = stage2Params.GetStrictMemoryNSortingThreadsPerSorters();

	if (!Params.sm_n_sorting_threads)
		Params.sm_n_sorting_threads = Params.n_threads;
	if (!Params.sm_n_uncompactors)
		Params.sm_n_uncompactors = 1;
	if (!Params.sm_n_mergers)
		Params.sm_n_mergers = 1;
}

//----------------------------------------------------------------------------------

template <unsigned SIZE> void CKMC<SIZE>::AdjustMemoryLimitsStrictMemoryMode()
{
	int64 m_rest = Params.max_mem_size;

	Params.sm_mem_part_input_file = 1ull << 26;
	Params.sm_mem_tot_input_file = Params.sm_mem_part_input_file * (Params.sm_n_uncompactors + 1);

	m_rest -= Params.sm_mem_tot_input_file;

	Params.sm_mem_part_expand = Params.sm_mem_part_input_file;
	Params.sm_mem_tot_expand = Params.sm_mem_part_expand * (Params.sm_n_uncompactors + 1); 

	m_rest -= Params.sm_mem_tot_expand;

	Params.sm_mem_part_suffixes = 1 << 25;
	Params.sm_mem_tot_suffixes = Params.sm_mem_part_suffixes * 2;

	m_rest -= Params.sm_mem_tot_suffixes;

	Params.sm_mem_part_lut = (1 << 2 * 12) * sizeof(uint64); //12 is max lut prefix len for strict memory sub bins
	Params.sm_mem_tot_lut = Params.sm_mem_part_lut * 2;

	m_rest -= Params.sm_mem_tot_lut;

	Params.sm_mem_part_merger_suff = 1 << 24;
	Params.sm_mem_tot_merger_suff = (Params.sm_n_mergers + 1) * Params.sm_mem_part_merger_suff;

	m_rest -= Params.sm_mem_tot_merger_suff;

	Params.sm_mem_part_merger_lut = 1 << 24;
	Params.sm_mem_tot_merger_lut = (Params.sm_n_mergers + 1) * Params.sm_mem_part_merger_lut; 
	m_rest -= Params.sm_mem_tot_merger_lut;

	const uint32 PREDICTET_NO_OF_SUBBINS = 3;

	Params.sm_mem_part_sub_bin_suff = (1ull << 24) * PREDICTET_NO_OF_SUBBINS;
	Params.sm_mem_tot_sub_bin_suff = Params.sm_mem_part_sub_bin_suff * Params.sm_n_mergers;
	
	m_rest -= Params.sm_mem_tot_sub_bin_suff;

	Params.sm_mem_part_sub_bin_lut = (1ull << 24) * PREDICTET_NO_OF_SUBBINS;
	Params.sm_mem_tot_sub_bin_lut = Params.sm_mem_part_sub_bin_lut * Params.sm_n_mergers;
	
	m_rest -= Params.sm_mem_tot_sub_bin_lut;

	//cout << "Memory left for sorter: " << (m_rest >> 20) << "MB"<< endl;

	Params.sm_mem_part_sort = m_rest;	
	Params.sm_mem_tot_sort = Params.sm_mem_part_sort;	
}

//----------------------------------------------------------------------------------
template <unsigned SIZE> void CKMC<SIZE>::AdjustMemoryLimitsStage2()
{
	// Memory for 2nd stage
	// Settings for memory manager of radix internal buffers
	Params.mem_part_pmm_radix_buf = (256 * BUFFER_WIDTHS[sizeof(CKmer<SIZE>)/8] + ALIGNMENT) * sizeof(CKmer<SIZE>);

	Params.mem_tot_pmm_radix_buf = Params.mem_part_pmm_radix_buf * Params.n_sorters * MAGIC_NUMBER;

	Params.max_mem_stage2 = Params.max_mem_size - Params.mem_tot_pmm_radix_buf;
}

//----------------------------------------------------------------------------------
// Adjust the memory limits for queues and other large data structures
template <unsigned SIZE> bool CKMC<SIZE>::AdjustMemoryLimits()
{
	// Memory for splitter internal buffers
	int64 m_rest = Params.max_mem_size;  

	Params.mem_part_pmm_stats = ((1 << Params.signature_len * 2) + 1) * sizeof(uint32);
	Params.mem_tot_pmm_stats = (Params.n_splitters + 1 + 1) * Params.mem_part_pmm_stats; //1 merged in main thread, 1 for sorting indices

	
	// Settings for memory manager of FASTQ buffers
	Params.fastq_buffer_size = 32 << 20;

	do {
		if(Params.fastq_buffer_size & (Params.fastq_buffer_size-1))
			Params.fastq_buffer_size &= Params.fastq_buffer_size - 1;
		else
			Params.fastq_buffer_size = Params.fastq_buffer_size / 2 + Params.fastq_buffer_size / 4;
		Params.mem_part_pmm_fastq = Params.fastq_buffer_size + CFastqReader::OVERHEAD_SIZE;
		Params.mem_tot_pmm_fastq  = Params.mem_part_pmm_fastq * (Params.n_readers + Params.n_splitters + 96);
	} while(Params.mem_tot_pmm_fastq > m_rest * 0.17);
	
	Params.mem_part_pmm_binary_file_reader = 1ll << 28; 
	Params.mem_tot_pmm_binary_file_reader = Params.mem_part_pmm_binary_file_reader * Params.n_readers * 3;
	while (Params.mem_tot_pmm_binary_file_reader > m_rest * 0.10)
	{
		Params.mem_part_pmm_binary_file_reader = Params.mem_part_pmm_binary_file_reader / 2 + Params.mem_part_pmm_binary_file_reader / 4;
		Params.mem_tot_pmm_binary_file_reader = Params.mem_part_pmm_binary_file_reader * Params.n_readers * 3;
	}

	if (Params.mem_part_pmm_binary_file_reader < (1ll << 23)) //8 MiB is a required minimum
	{
		Params.mem_part_pmm_binary_file_reader = 1ll << 23;
		Params.mem_tot_pmm_binary_file_reader = Params.mem_part_pmm_binary_file_reader * Params.n_readers * 3;
	}

	m_rest -= Params.mem_tot_pmm_fastq;
	m_rest -= Params.mem_tot_pmm_binary_file_reader;

	// Subtract memory for bin collectors internal buffers
	m_rest -= Params.n_splitters * Params.bin_part_size * sizeof(CKmer<SIZE>);

	// Settings for memory manager of reads
	Params.mem_part_pmm_reads = (CSplitter::MAX_LINE_SIZE + 1) * sizeof(double);
	Params.mem_tot_pmm_reads  = Params.mem_part_pmm_reads * 2 * Params.n_splitters;
	m_rest -= Params.mem_tot_pmm_reads;

	

	Params.mem_part_pmm_bins = Params.bin_part_size;

	Params.mem_tot_pmm_bins = m_rest;

	
	// memory for storer internal buffer
	if(Params.max_mem_size >= 16ll << 30)
		Params.max_mem_storer = (int64) (Params.mem_tot_pmm_bins * 0.75);
	else
		Params.max_mem_storer = (int64) (Params.mem_tot_pmm_bins * 0.65);

	if(Params.max_mem_storer < (1ll << 28))
		return false;

	//each splitter has n_bins collectors, 
	uint64 pmm_bins_req_parts = Params.n_splitters * Params.n_bins;
	uint64 available_pmm_bins_parts = (Params.mem_tot_pmm_bins - Params.max_mem_storer) / Params.mem_part_pmm_bins;
	while (available_pmm_bins_parts < pmm_bins_req_parts)
	{
		Params.bin_part_size = Params.bin_part_size / 2 + Params.bin_part_size / 4;
		Params.mem_part_pmm_bins = Params.bin_part_size;
		available_pmm_bins_parts = (Params.mem_tot_pmm_bins - Params.max_mem_storer) / Params.mem_part_pmm_bins;
	}

	// Max. memory for single package
	Params.max_mem_storer_pkg = 2 * Params.max_mem_storer / Params.n_bins;

	return true;
}

//----------------------------------------------------------------------------------
// Show the settings of the KMC (in verbose mode only)
template <unsigned SIZE> void CKMC<SIZE>::ShowSettingsStage1()
{
	std::ostringstream ostr;

	ostr << "\n********** Used parameters for Stage 1 : **********\n";

	ostr << "No. of input files           : " << Params.input_file_names.size() << "\n";
	ostr << "Output file name             : " << Params.output_file_name << "\n";
	ostr << "No. of working directories   : " << 1 << "\n";
	ostr << "Input format                 : ";
	switch (Params.file_type)
	{
	case InputType::FASTA:
		ostr << "FASTA\n";
		break;
	case InputType::FASTQ:
		ostr << "FASTQ\n";
		break;
	case InputType::MULTILINE_FASTA:
		ostr << "MULTI LINE FASTA\n";
		break;
	case InputType::BAM:
		ostr << "BAM\n";
		break;
	}
	ostr << "Output format                : ";
	switch (Params.output_type)		
	{
	case OutputType::KFF:
		ostr << "KFF\n";
		break;
	case OutputType::KMC:
		ostr << "KMC\n";
		break;	
	}
	ostr << "\n";
	ostr << "k-mer length                 : " << Params.kmer_len << "\n";
	ostr << "Max. k-mer length            : " << MAX_K << "\n";
	ostr << "Signature length             : " << Params.signature_len << "\n";
	ostr << "Both strands                 : " << (Params.both_strands ? "true\n" : "false\n");
	ostr << "RAM only mode                : " << (Params.mem_mode ? "true\n" : "false\n");

	ostr << "\n******* Stage 1 configuration: *******\n";
	ostr << "\n";
	ostr << "No. of bins                  : " << Params.n_bins << "\n";
	ostr << "Bin part size                : " << Params.bin_part_size << "\n";
	ostr << "Input buffer size            : " << Params.fastq_buffer_size << "\n";
	
	ostr << "\n";

	ostr << "No. of readers               : " << Params.n_readers << "\n";
	ostr << "No. of splitters             : " << Params.n_splitters << "\n";
	ostr << "\n";

	ostr << "Max. mem. size               : " << setw(5) << (Params.max_mem_size / 1000000) << "MB\n";
	ostr << "Max. mem. per storer         : " << setw(5) << (Params.max_mem_storer / 1000000) << "MB\n";
	ostr << "Max. mem. for single package : " << setw(5) << (Params.max_mem_storer_pkg / 1000000) << "MB\n";
	ostr << "\n";

	ostr << "Max. mem. for PMM (bin parts): " << setw(5) << (Params.mem_tot_pmm_bins / 1000000) << "MB\n";
	ostr << "Max. mem. for PMM (FASTQ)    : " << setw(5) << (Params.mem_tot_pmm_fastq / 1000000) << "MB\n";
	ostr << "Max. mem. for PMM (reads)    : " << setw(5) << (Params.mem_tot_pmm_reads / 1000000) << "MB\n";
	ostr << "Max. mem. for PMM (b. reader): " << setw(5) << (Params.mem_tot_pmm_binary_file_reader / 1000000) << "MB\n";

	ostr << "\n";

	Params.verboseLogger->Log(ostr.str());
}

//----------------------------------------------------------------------------------
// Show the settings of the KMC (in verbose mode only)
template <unsigned SIZE> void CKMC<SIZE>::ShowSettingsStage2()
{
	std::ostringstream ostr;

	ostr << "\n********** Used parameters for Stage 2 : **********\n";
	ostr << "Min. count threshold         : " << Params.cutoff_min << "\n";
	ostr << "Max. count threshold         : " << Params.cutoff_max << "\n";
	ostr << "Max. counter value           : " << Params.counter_max << "\n";

	ostr << "\n******* Stage 2 configuration: *******\n";

	ostr << "No. of threads               : " << Params.n_sorters << "\n";
	
	ostr << "\n";

	ostr << "Max. mem. for 2nd stage      : " << setw(5) << (Params.max_mem_stage2 / 1000000) << "MB\n";
	ostr << "\n";

	Params.verboseLogger->Log(ostr.str());
}

//----------------------------------------------------------------------------------
// Show the settings of the KMC (in verbose mode only)
template <unsigned SIZE> void CKMC<SIZE>::ShowSettingsSmallKOpt()
{
	std::ostringstream ostr;

	ostr << "\n******* configuration for small k mode: *******\n";

	ostr << "No. of input files           : " << Params.input_file_names.size() << "\n";
	ostr << "Output file name             : " << Params.output_file_name << "\n";
	ostr << "Input format                 : ";
	switch (Params.file_type)
	{
	case InputType::FASTA:
		ostr << "FASTA\n";
		break;
	case InputType::FASTQ:
		ostr << "FASTQ\n";
		break;
	case InputType::MULTILINE_FASTA:
		ostr << "MULTI LINE FASTA\n";
		break;
	case InputType::BAM:
		ostr << "BAM\n";
		break;
	}
	ostr << "Output format                 : ";
	switch (Params.output_type)
	{
	case OutputType::KFF:
		ostr << "KFF\n";
		break;
	case OutputType::KMC:
		ostr << "KMC\n";
		break;
	}

	ostr << "\n";
	ostr << "k-mer length                 : " << Params.kmer_len << "\n";
	ostr << "Max. k-mer length            : " << MAX_K << "\n";
	ostr << "Min. count threshold         : " << Params.cutoff_min << "\n";
	ostr << "Max. count threshold         : " << Params.cutoff_max << "\n";
	ostr << "Max. counter value           : " << Params.counter_max << "\n";
	
	ostr << "Both strands                 : " << (Params.both_strands ? "true\n" : "false\n");
	ostr << "Input buffer size            : " << Params.fastq_buffer_size << "\n";

	ostr << "\n";

	ostr << "No. of readers               : " << Params.n_readers << "\n";
	ostr << "No. of splitters             : " << Params.n_splitters << "\n";
	ostr << "\n";

	ostr << "Max. mem. size               : " << setw(5) << (Params.max_mem_size / 1000000) << "MB\n";
	ostr << "\n";

	ostr << "Max. mem. for PMM (FASTQ)    : " << setw(5) << (Params.mem_tot_pmm_fastq / 1000000) << "MB\n";
	ostr << "Part. mem. for PMM (FASTQ)   : " << setw(5) << (Params.mem_part_pmm_fastq / 1000000) << "MB\n";
	ostr << "Max. mem. for PMM (reads)    : " << setw(5) << (Params.mem_tot_pmm_reads / 1000000) << "MB\n";
	ostr << "Part. mem. for PMM (reads)   : " << setw(5) << (Params.mem_part_pmm_reads / 1000000) << "MB\n";
	ostr << "Max. mem. for PMM (b. reader): " << setw(5) << (Params.mem_tot_pmm_binary_file_reader / 1000000) << "MB\n";
	ostr << "Part. mem. for PMM (b. reader): " << setw(5) << (Params.mem_part_pmm_binary_file_reader / 1000000) << "MB\n";

	ostr << "\n";

	Params.verboseLogger->Log(ostr.str());
}
//----------------------------------------------------------------------------------
template <unsigned SIZE> bool CKMC<SIZE>::AdjustMemoryLimitsSmallK() 
{
	if (Params.kmer_len > 13) 
		return false;

	bool small_k_opt_required = Params.kmer_len < Params.signature_len;	

	int tmp_n_splitters = Params.n_splitters;
	int tmp_n_readers = Params.n_readers;
	int tmp_fastq_buffer_size = 0;
	int64 tmp_mem_part_pmm_fastq = 0;
	int64 tmp_mem_tot_pmm_fastq = 0;
	int64 tmp_mem_part_pmm_reads = (CSplitter::MAX_LINE_SIZE + 1) * sizeof(double);
	int64 tmp_mem_tot_pmm_reads = 0;
	int64 tmp_mem_part_pmm_binary_file_reader = 0;
	int64 tmp_mem_tot_pmm_binary_file_reader = 0;

	int64 tmp_mem_part_small_k_buf = (1ll << 2 * Params.kmer_len) * sizeof(uint64_t);//no of possible k-mers * counter size
	int64 tmp_mem_tot_small_k_buf = 0;

	int64 additional_buffers = 96;

	while (true)
	{
		int64 mim_mem_for_readers = tmp_n_readers * (16 << 20);
		tmp_fastq_buffer_size = 32 << 20; //important for bam reading performance, previously 1 << 24
		tmp_mem_part_pmm_fastq = tmp_fastq_buffer_size + CFastqReader::OVERHEAD_SIZE;
		tmp_mem_tot_pmm_fastq = tmp_mem_part_pmm_fastq * (tmp_n_readers + tmp_n_splitters + additional_buffers);

		tmp_mem_part_pmm_binary_file_reader = 1ll << 27;
		tmp_mem_tot_pmm_binary_file_reader = tmp_mem_part_pmm_binary_file_reader * tmp_n_readers * 3;

		int64 mem_rest = Params.max_mem_size - tmp_mem_tot_pmm_fastq - tmp_mem_tot_pmm_binary_file_reader;

		tmp_mem_tot_pmm_reads = tmp_mem_part_pmm_reads * 3 * tmp_n_splitters;
		tmp_mem_tot_small_k_buf = tmp_mem_part_small_k_buf * tmp_n_splitters;

		if (tmp_mem_tot_pmm_reads + tmp_mem_tot_small_k_buf + mim_mem_for_readers < mem_rest)
			break;

		if (additional_buffers)
			additional_buffers = additional_buffers / 2 + additional_buffers / 4;
		else if (tmp_n_readers < tmp_n_splitters)
			--tmp_n_splitters;
		else
			--tmp_n_readers;

		if (!tmp_n_readers || !tmp_n_splitters)
		{
			if (small_k_opt_required)
			{
				//Should never be here
				cerr << "Error: Internal error occurred during small k adjustment, please report this via https://github.com/refresh-bio/KMC/issues";
				exit(1);
			}
			return false;
		}
	}

	Params.n_splitters = tmp_n_splitters;
	Params.n_readers = tmp_n_readers;
	Params.fastq_buffer_size = tmp_fastq_buffer_size;
	Params.mem_part_pmm_fastq = tmp_mem_part_pmm_fastq;
	Params.mem_part_small_k_completer = Params.mem_tot_small_k_completer = Params.mem_tot_pmm_fastq = tmp_mem_tot_pmm_fastq;
	Params.mem_part_pmm_binary_file_reader = tmp_mem_part_pmm_binary_file_reader;
	Params.mem_tot_pmm_binary_file_reader = tmp_mem_tot_pmm_binary_file_reader;
	Params.mem_part_pmm_reads = tmp_mem_part_pmm_reads;
	Params.mem_tot_pmm_reads = tmp_mem_tot_pmm_reads;
	Params.mem_part_small_k_buf = tmp_mem_part_small_k_buf;
	Params.mem_tot_small_k_buf = tmp_mem_tot_small_k_buf;
	
	return true;
}

//----------------------------------------------------------------------------------
template <unsigned SIZE>
void CKMC<SIZE>::CheckAndReportMissingEOLs()
{
	auto n_missing = Queues.missingEOL_at_EOF_counter->Get();
	if (n_missing)
	{
		std::ostringstream ostr;
		ostr << "in " << n_missing << " input file(s) there was not end of line character at EOF.";
		Params.warningsLogger->Log(ostr.str());
	}
}

//----------------------------------------------------------------------------------
template <unsigned SIZE>
KMC::Stage1Results CKMC<SIZE>::ProcessSmallKOptimization_Stage1()
{
	KMC::Stage1Results results;
	results.wasSmallKOptUsed = true;
	ShowSettingsSmallKOpt();

	w1.startTimer();
	Queues.input_files_queue = new CInputFilesQueue(Params.input_file_names);
	Queues.part_queue = new CPartQueue(Params.n_readers);

	Queues.pmm_fastq = new CMemoryPoolWithBamSupport(Params.mem_tot_pmm_fastq, Params.mem_part_pmm_fastq);
	Queues.pmm_binary_file_reader = new CMemoryPool(Params.mem_tot_pmm_binary_file_reader, Params.mem_part_pmm_binary_file_reader);
	Queues.pmm_reads = new CMemoryPool(Params.mem_tot_pmm_reads, Params.mem_part_pmm_reads);
	Queues.pmm_small_k_buf = new CMemoryPool(Params.mem_tot_small_k_buf, Params.mem_part_small_k_buf);

	Queues.missingEOL_at_EOF_counter = new CMissingEOL_at_EOF_counter;

	w_small_k_splitters.resize(Params.n_splitters);

	for (int i = 0; i < Params.n_splitters; ++i)
	{
		w_small_k_splitters[i] = new CWSmallKSplitter<uint64_t>(Params, Queues);
		gr1_2.push_back(thread(std::ref(*w_small_k_splitters[i])));
	}

	w_fastqs.resize(Params.n_readers);
	if (Params.file_type != InputType::BAM)
	{
		Queues.binary_pack_queues.resize(Params.n_readers);
		for (int i = 0; i < Params.n_readers; ++i)
		{
			Queues.binary_pack_queues[i] = new CBinaryPackQueue;
			w_fastqs[i] = new CWFastqReader(Params, Queues, Queues.binary_pack_queues[i]);
			gr1_1.push_back(thread(std::ref(*w_fastqs[i])));
		}
	}
	else
	{
		Queues.bam_task_manager = new CBamTaskManager;
		for (int i = 0; i < Params.n_readers; ++i)
		{
			w_fastqs[i] = new CWFastqReader(Params, Queues, nullptr);
			gr1_1.push_back(thread(std::ref(*w_fastqs[i])));
		}
	}

	w_bin_file_reader = new CWBinaryFilesReader(Params, Queues);
	thread bin_file_reader_th(std::ref(*w_bin_file_reader));

	for (auto& t : gr1_1)
		t.join();

	for (auto& t : gr1_2)
		t.join();

	for (auto r : w_fastqs)
		delete r;

	bin_file_reader_th.join();
	delete w_bin_file_reader;

	for (auto& ptr : Queues.binary_pack_queues)
		delete ptr;

	if (Params.file_type == InputType::BAM)
		delete Queues.bam_task_manager;

	uint64 tmp_n_reads;
	n_reads = 0;
	results.nTotalSuperKmers = 0;
	for (auto s : w_small_k_splitters)
	{
		s->GetTotal(tmp_n_reads);
		n_reads += tmp_n_reads;
	}
	results.nSeqences = n_reads;

	w1.stopTimer();

	results.time = w1.getElapsedTime();
	results.tmpSize = 0;

	return results;
}
//----------------------------------------------------------------------------------
template <unsigned SIZE>
KMC::Stage2Results CKMC<SIZE>::ProcessSmallKOptimization_Stage2()
{
	KMC::Stage2Results results;
	w2.startTimer();

	vector<CSmallKBuf<uint64_t>> count_results(Params.n_splitters);

	for (int i = 0; i < Params.n_splitters; ++i)
	{
		count_results[i] = w_small_k_splitters[i]->GetResult();
	}

	uint64 n_kmers = 0;

	for (int j = 1; j < Params.n_splitters; ++j)
	{
		for (int i = 0; i < (1 << 2 * Params.kmer_len); ++i)
			count_results[0].buf[i] += count_results[j].buf[i];
	}

	results.nTotalKmers = 0;


	for (int j = 0; j < (1 << 2 * Params.kmer_len); ++j)
		if (count_results[0].buf[j]) ++n_kmers;

	for (auto s : w_small_k_splitters)
	{
		results.nTotalKmers += s->GetTotalKmers();
		s->Release();
		delete s;
	}


	Queues.pmm_fastq->release();
	delete Queues.pmm_fastq;
	delete Queues.pmm_binary_file_reader;


	uint32 best_lut_prefix_len = 0;
	uint64 best_mem_amount = 1ull << 62;

	if (Params.output_type == OutputType::KMC)
	{
		uint32 counter_size = 0;

		counter_size = calc_counter_size(Params.cutoff_max, Params.counter_max);

		for (Params.lut_prefix_len = 1; Params.lut_prefix_len < 16; ++Params.lut_prefix_len)
		{
			uint32 suffix_len;
			if (Params.lut_prefix_len > (uint32)Params.kmer_len)
				suffix_len = 0;
			else
				suffix_len = Params.kmer_len - Params.lut_prefix_len;

			if (suffix_len % 4)
				continue;

			uint64 suf_mem = n_kmers * (suffix_len / 4 + counter_size);
			uint64 lut_mem = (1ull << (2 * Params.lut_prefix_len)) * sizeof(uint64);

			if (suf_mem + lut_mem < best_mem_amount)
			{
				best_lut_prefix_len = Params.lut_prefix_len;
				best_mem_amount = suf_mem + lut_mem;
			}
		}

		Params.lut_prefix_len = best_lut_prefix_len;
	}
	else if (Params.output_type == OutputType::KFF)
		Params.lut_prefix_len = 0;
	else
	{
		std::cerr << "Error: not implemented, plase contact authors showind this message" << __FILE__ << "\t" << __LINE__ << "\n";
		exit(1);
	}

	Queues.pmm_small_k_completer = new CMemoryPool(Params.mem_tot_small_k_completer, Params.mem_part_small_k_completer);

	CSmallKCompleter small_k_completer(Params, Queues);
	small_k_completer.Complete(count_results[0]);
	small_k_completer.GetTotal(results.nUniqueKmers, results.nBelowCutoffMin, results.nAboveCutoffMax);

	Queues.pmm_reads->release();
	Queues.pmm_small_k_buf->release();
	Queues.pmm_small_k_completer->release();
	delete Queues.pmm_small_k_completer;
	delete Queues.pmm_reads;
	delete Queues.pmm_small_k_buf;
	w2.stopTimer();
	//cerr << "\n";


	results.maxDiskUsage = 0;

	CheckAndReportMissingEOLs();
	delete Queues.missingEOL_at_EOF_counter;

	return results;
}

//----------------------------------------------------------------------------------
// Run the counter stage 1
template <unsigned SIZE> KMC::Stage1Results CKMC<SIZE>::ProcessStage1()
{
	w1.startTimer();
	KMC::Stage1Results results{};
	if (!initialized)
		throw std::runtime_error("Error: kmc was not initialized!\n");

	was_small_k_opt = false;
	if (AdjustMemoryLimitsSmallK())
	{
		was_small_k_opt = true;

		Params.verboseLogger->Log("\nInfo: Small k optimization on!\n");

		return ProcessSmallKOptimization_Stage1();
	}

	if (!AdjustMemoryLimits())
		throw std::runtime_error("Error: cannot adjust memory, please contact authors");

	// Create monitors
	Queues.mm = new CMemoryMonitor(Params.max_mem_stage2);


	// Create queues
	Queues.input_files_queue = new CInputFilesQueue(Params.input_file_names);
	Queues.part_queue = new CPartQueue(Params.n_readers);
	Queues.bpq = new CBinPartQueue(Params.n_splitters);
	Queues.bd = new CBinDesc;
	Queues.epd = new CExpanderPackDesc(Params.n_bins);
	Queues.bq = new CBinQueue(1);


	// Create memory manager
	Queues.pmm_binary_file_reader = new CMemoryPool(Params.mem_tot_pmm_binary_file_reader, Params.mem_part_pmm_binary_file_reader);
	Queues.pmm_bins = new CMemoryPool(Params.mem_tot_pmm_bins, Params.mem_part_pmm_bins);
	Queues.pmm_fastq = new CMemoryPoolWithBamSupport(Params.mem_tot_pmm_fastq, Params.mem_part_pmm_fastq);
	Queues.pmm_reads = new CMemoryPool(Params.mem_tot_pmm_reads, Params.mem_part_pmm_reads);
	Queues.pmm_stats = new CMemoryPool(Params.mem_tot_pmm_stats, Params.mem_part_pmm_stats);

	Queues.missingEOL_at_EOF_counter = new CMissingEOL_at_EOF_counter;

	if (Params.file_type != InputType::BAM)
	{
		Queues.binary_pack_queues.resize(Params.n_readers);
		for (int i = 0; i < Params.n_readers; ++i)
		{
			Queues.binary_pack_queues[i] = new CBinaryPackQueue;
		}
	}
	else //for bam
	{
		Queues.bam_task_manager = new CBamTaskManager;
	}
	w_bin_file_reader = new CWBinaryFilesReader(Params, Queues, false);
	Queues.stats_part_queue = new CStatsPartQueue(Params.n_readers, MAX(STATS_FASTQ_SIZE, w_bin_file_reader->GetPredictedSize() / 100));

	Queues.s_mapper = new CSignatureMapper(Queues.pmm_stats, Params.signature_len, Params.n_bins
#ifdef DEVELOP_MODE
		, Params.verbose_log
#endif
	);
	Queues.disk_logger = new CDiskLogger;
	// ***** Stage 0 *****
	w0.startTimer();
	w_stats_splitters.resize(Params.n_splitters);

	for (int i = 0; i < Params.n_splitters; ++i)
	{
		w_stats_splitters[i] = new CWStatsSplitter(Params, Queues);
		gr0_2.push_back(thread(std::ref(*w_stats_splitters[i])));
	}

	w_stats_fastqs.resize(Params.n_readers);
	for (int i = 0; i < Params.n_readers; ++i)
	{
		w_stats_fastqs[i] = new CWStatsFastqReader(Params, Queues, Params.file_type == InputType::BAM ? nullptr : Queues.binary_pack_queues[i]);
		gr0_1.push_back(thread(std::ref(*w_stats_fastqs[i])));
	}
	thread bin_file_reader_th(std::ref(*w_bin_file_reader));

	for (auto p = gr0_1.begin(); p != gr0_1.end(); ++p)
		p->join();

	for (auto p = gr0_2.begin(); p != gr0_2.end(); ++p)
		p->join();

	bin_file_reader_th.join();
	delete w_bin_file_reader;

	for (auto& ptr : Queues.binary_pack_queues)
		delete ptr;

	if (Params.file_type == InputType::BAM)
		delete Queues.bam_task_manager;

	uint32* stats;
	Queues.pmm_stats->reserve(stats);
	fill_n(stats, (1 << Params.signature_len * 2) + 1, 0);


	for (int i = 0; i < Params.n_readers; ++i)
		delete w_stats_fastqs[i];

	for (int i = 0; i < Params.n_splitters; ++i)
	{
		w_stats_splitters[i]->GetStats(stats);
		delete w_stats_splitters[i];
	}

	delete Queues.stats_part_queue;
	Queues.stats_part_queue = nullptr;
	delete Queues.input_files_queue;
	Queues.input_files_queue = new CInputFilesQueue(Params.input_file_names);

	heuristic_time.startTimer();
	Queues.s_mapper->Init(stats);
	heuristic_time.stopTimer();

	cerr << "\n";
	w0.stopTimer();


	Queues.pmm_stats->free(stats);
	Queues.pmm_stats->release();
	delete Queues.pmm_stats;
	Queues.pmm_stats = nullptr;

	// ***** Stage 1 *****
	ShowSettingsStage1();
	Queues.missingEOL_at_EOF_counter->Reset();

	w_splitters.resize(Params.n_splitters);

	for (int i = 0; i < Params.n_splitters; ++i)
	{
		w_splitters[i] = new CWSplitter(Params, Queues);
		gr1_2.push_back(thread(std::ref(*w_splitters[i])));
	}

	w_storer = new CWKmerBinStorer(Params, Queues);
	gr1_3.push_back(thread(std::ref(*w_storer)));

	w_fastqs.resize(Params.n_readers);
	if (Params.file_type != InputType::BAM)
	{
		Queues.binary_pack_queues.resize(Params.n_readers);
		for (int i = 0; i < Params.n_readers; ++i)
		{
			Queues.binary_pack_queues[i] = new CBinaryPackQueue;
			w_fastqs[i] = new CWFastqReader(Params, Queues, Queues.binary_pack_queues[i]);
			gr1_1.push_back(thread(std::ref(*w_fastqs[i])));
		}
	}
	else //bam
	{
		Queues.bam_task_manager = new CBamTaskManager;
		for (int i = 0; i < Params.n_readers; ++i)
		{
			w_fastqs[i] = new CWFastqReader(Params, Queues, nullptr);
			gr1_1.push_back(thread(std::ref(*w_fastqs[i])));
		}
	}

	w_bin_file_reader = new CWBinaryFilesReader(Params, Queues);
	bin_file_reader_th = thread(std::ref(*w_bin_file_reader));

	for (auto p = gr1_1.begin(); p != gr1_1.end(); ++p)
		p->join();
	for (auto p = gr1_2.begin(); p != gr1_2.end(); ++p)
		p->join();

	bin_file_reader_th.join();
	delete w_bin_file_reader;

	for (auto& ptr : Queues.binary_pack_queues)
		delete ptr;

	if (Params.file_type == InputType::BAM)
		delete Queues.bam_task_manager;

	Queues.pmm_fastq->release();
	Queues.pmm_reads->release();

	delete Queues.pmm_fastq;
	delete Queues.pmm_reads;
	delete Queues.pmm_binary_file_reader;

	for (auto p = gr1_3.begin(); p != gr1_3.end(); ++p)
		p->join();

	n_reads = 0;

	thread* release_thr_st1_1 = new thread([&] {
		for (int i = 0; i < Params.n_readers; ++i)
			delete w_fastqs[i];

		for (int i = 0; i < Params.n_splitters; ++i)
		{
			uint64 _n_reads;
			w_splitters[i]->GetTotal(_n_reads);
			n_reads += _n_reads;
			delete w_splitters[i];
		}

		delete w_storer;
	});

	thread* release_thr_st1_2 = new thread([&] {
		Queues.pmm_bins->release();
		delete Queues.pmm_bins;
	});


	release_thr_st1_1->join();
	results.nSeqences = n_reads;
	release_thr_st1_2->join();

	delete release_thr_st1_1;
	delete release_thr_st1_2;

	// ***** Getting disk usage statistics *****

	results.tmpSize = 0;
	results.nTotalSuperKmers = 0;
	int32 bin_id;
	CMemDiskFile* file;
	string name;
	uint64 size;
	uint64 n_rec;
	uint64 n_plus_x_recs;
	uint64 n_super_kmers;
	while ((bin_id = Queues.bd->get_next_bin()) >= 0)
	{
		Queues.bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, n_super_kmers);
		results.tmpSize += size;
		results.nTotalSuperKmers += n_super_kmers;
	}

	w1.stopTimer();

	CheckAndReportMissingEOLs();
	delete Queues.missingEOL_at_EOF_counter;

	// ***** End of Stage 1 *****
	results.time = w1.getElapsedTime();
	return results;
}
//----------------------------------------------------------------------------------
// Run the counter stage 2
template <unsigned SIZE> KMC::Stage2Results CKMC<SIZE>::ProcessStage2()
{
	KMC::Stage2Results results;
	if (was_small_k_opt)
	{
		return ProcessSmallKOptimization_Stage2();
	}
	w2.startTimer();

	if (!initialized)
		throw std::runtime_error("Error: kmc was not initialized!\n");

	int32 bin_id;
	CMemDiskFile* file;
	string name;
	uint64 size;
	uint64 n_rec;
	uint64 n_plus_x_recs;
	uint64 n_super_kmers;


	// Adjust RAM for 2nd stage
	// Calculate LUT size

	if (Params.output_type == OutputType::KMC)
	{
		uint32 best_lut_prefix_len = 0;
		uint64 best_mem_amount = 1ull << 62;

		for (Params.lut_prefix_len = 2; Params.lut_prefix_len < 16; ++Params.lut_prefix_len)
		{
			uint32 suffix_len = Params.kmer_len - Params.lut_prefix_len;
			if (suffix_len % 4)
				continue;

			uint64 est_suf_mem = n_reads * suffix_len;
			uint64 lut_mem = Params.n_bins * (1ull << (2 * Params.lut_prefix_len)) * sizeof(uint64);

			if (est_suf_mem + lut_mem < best_mem_amount)
			{
				best_lut_prefix_len = Params.lut_prefix_len;
				best_mem_amount = est_suf_mem + lut_mem;
			}
		}

		Params.lut_prefix_len = best_lut_prefix_len;
	}
	else if (Params.output_type == OutputType::KFF)
		Params.lut_prefix_len = 0;
	else
	{
		std::cerr << "Error: not implemented, plase contact authors showind this message" << __FILE__ << "\t" << __LINE__ << "\n";
		exit(1);
	}

	Queues.bd->reset_reading();
	vector<int64> bin_sizes;
	while ((bin_id = Queues.bd->get_next_bin()) >= 0)
	{
		Queues.bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, n_super_kmers);
		if (Params.max_x)
			bin_sizes.push_back(n_plus_x_recs * 2 * sizeof(CKmer<SIZE>));			// estimation of RAM for sorting bins
		else
			bin_sizes.push_back(n_rec * 2 * sizeof(CKmer<SIZE>));
	}

	sort(bin_sizes.begin(), bin_sizes.end(), greater<int64>());
	
	AdjustMemoryLimitsStage2();

	if (Params.use_strict_mem)
	{
		Queues.tlbq = new CTooLargeBinsQueue;
		Queues.bbkpq = new CBigBinKmerPartQueue(Params.sm_n_mergers);
	}
	else
	{
		Queues.tlbq = nullptr;
		Queues.bbkpq = nullptr;
	}


	int64 stage2_size = 0;
	for (auto bin_size : bin_sizes)
		stage2_size += bin_size;
	stage2_size = MAX(stage2_size, 16 << 20);
	Params.max_mem_stage2 = MIN(Params.max_mem_stage2, stage2_size);

	ShowSettingsStage2();

	// ***** Stage 2 *****
	Queues.bd->reset_reading();
	Queues.pmm_radix_buf = new CMemoryPool(Params.mem_tot_pmm_radix_buf, Params.mem_part_pmm_radix_buf);
	Queues.memory_bins = new CMemoryBins(Params.max_mem_stage2, Params.n_bins, Params.use_strict_mem, Params.n_threads);

	auto sorted_bins = Queues.bd->get_sorted_req_sizes(Params.max_x, sizeof(CKmer<SIZE>), Params.cutoff_min, Params.cutoff_max, Params.counter_max, Params.lut_prefix_len);
	Queues.bd->init_sort(sorted_bins);

#ifdef DEVELOP_MODE
	if (Params.verbose_log)
		save_bins_stats(Queues, Params, sizeof(CKmer<SIZE>), n_reads, Params.signature_len, Queues.s_mapper->GetMapSize(), Queues.s_mapper->GetMap());
#endif

	SortFunction<CKmer<SIZE>> sort_func;
#ifdef __APPLE__
	sort_func = RadixSort::RadixSortMSD<CKmer<SIZE>, SIZE>;
	CSmallSort<SIZE>::Adjust(384);
#else	
	auto proc_name = CCpuInfo::GetBrand();
	bool is_intel = CCpuInfo::GetVendor() == "GenuineIntel";
	bool at_least_avx = CCpuInfo::AVX_Enabled();
	std::transform(proc_name.begin(), proc_name.end(), proc_name.begin(), ::tolower);
	bool is_xeon = proc_name.find("xeon") != string::npos;
	if (is_xeon || (is_intel && at_least_avx))
	{
		if (CCpuInfo::AVX2_Enabled())
			sort_func = RadulsSort::RadixSortMSD_AVX2<CKmer<SIZE>>;
		else if (CCpuInfo::AVX_Enabled())
			sort_func = RadulsSort::RadixSortMSD_AVX<CKmer<SIZE>>;
		else if (CCpuInfo::SSE41_Enabled())
			sort_func = RadulsSort::RadixSortMSD_SSE41<CKmer<SIZE>>;
		else if (CCpuInfo::SSE2_Enabled())
			sort_func = RadulsSort::RadixSortMSD_SSE2<CKmer<SIZE>>;
		else
		{
			//probably never because x64 always supports sse2 as far as I know
			std::cerr << "Error: At least SSE2 must be supported\n";
		}
	}
	else
	{
		sort_func = RadixSort::RadixSortMSD<CKmer<SIZE>, SIZE>;
		CSmallSort<SIZE>::Adjust(384);
	}
#endif

	{
		auto _2nd_stage_threads = Params.n_sorters;
		Queues.kq = new CKmerQueue(Params.n_bins, _2nd_stage_threads);
		//max memory is specified by user
		int64 max_mem_size = Queues.memory_bins->GetTotalSize();

		//but in some cases we will need more memory to process some big bins
		//so max memory size will be higher
		//but it is not true in strict memory mode, where such big bins are processed after stage 2
		if (max_mem_size < sorted_bins.front().second && !Params.use_strict_mem)
			max_mem_size = sorted_bins.front().second;

		Queues.sorters_manager = new CSortersManager(Params.n_bins, _2nd_stage_threads, Queues.bq, max_mem_size, sorted_bins);

		w_sorters.resize(_2nd_stage_threads);

		for (int i = 0; i < _2nd_stage_threads; ++i)
		{
			w_sorters[i] = new CWKmerBinSorter<SIZE>(Params, Queues, sort_func);
			gr2_2.push_back(thread(std::ref(*w_sorters[i])));
		}
	}

	w_reader = new CWKmerBinReader<SIZE>(Params, Queues);
	gr2_1.push_back(thread(std::ref(*w_reader)));

	w_completer = new CWKmerBinCompleter(Params, Queues);
	gr2_3.push_back(thread(std::ref(*w_completer), true));



	for (auto p = gr2_1.begin(); p != gr2_1.end(); ++p)
		p->join();
	for (auto p = gr2_2.begin(); p != gr2_2.end(); ++p)
		p->join();

	//Finishing first stage of completer
	for (auto p = gr2_3.begin(); p != gr2_3.end(); ++p)
		p->join();

	gr2_3.clear();

	thread* release_thr_st2_1 = new thread([&] {
		delete Queues.mm;
		//Queues.pmm_radix_buf->release();
		Queues.memory_bins->release();
		//delete Queues.pmm_radix_buf;
		delete Queues.memory_bins;
	});

	//process big bins if necessary (only in strict memory limit mode)
	thread* release_thr_sm = nullptr;

	if (Params.use_strict_mem)
	{
		w2.stopTimer();
		results.time = w2.getElapsedTime();
		w3.startTimer();
		release_thr_st2_1->join(); //need to be sure that memory_bins is released
		AdjustMemoryLimitsStrictMemoryMode();

		cerr << "\n";

		Queues.sm_pmm_input_file = new CMemoryPool(Params.sm_mem_tot_input_file, Params.sm_mem_part_input_file);
		Queues.sm_pmm_expand = new CMemoryPool(Params.sm_mem_tot_expand, Params.sm_mem_part_expand);
		Queues.sm_pmm_sort = new CMemoryPool(Params.sm_mem_tot_sort, Params.sm_mem_part_sort);
		Queues.sm_pmm_sorter_suffixes = new CMemoryPool(Params.sm_mem_tot_suffixes, Params.sm_mem_part_suffixes);
		Queues.sm_pmm_sorter_lut = new CMemoryPool(Params.sm_mem_tot_lut, Params.sm_mem_part_lut);

		Queues.sm_pmm_merger_lut = new CMemoryPool(Params.sm_mem_tot_merger_lut, Params.sm_mem_part_merger_lut);
		Queues.sm_pmm_merger_suff = new CMemoryPool(Params.sm_mem_tot_merger_suff, Params.sm_mem_part_merger_suff);
		Queues.sm_pmm_sub_bin_lut = new CMemoryPool(Params.sm_mem_tot_sub_bin_lut, Params.sm_mem_part_sub_bin_lut);
		Queues.sm_pmm_sub_bin_suff = new CMemoryPool(Params.sm_mem_tot_sub_bin_suff, Params.sm_mem_part_sub_bin_suff);

		Queues.bbpq = new CBigBinPartQueue();
		Queues.bbkq = new CBigBinKXmersQueue(Params.sm_n_uncompactors);
		Queues.bbd = new CBigBinDesc();
		Queues.bbspq = new CBigBinSortedPartQueue(1);
		Queues.sm_cbc = new CCompletedBinsCollector(1);

		CWBigKmerBinReader* w_bkb_reader = new CWBigKmerBinReader(Params, Queues);
		thread bkb_reader(std::ref(*w_bkb_reader));

		vector<CWBigKmerBinUncompactor<SIZE>*> w_bkb_uncompactors(Params.sm_n_uncompactors);
		vector<thread> bkb_uncompactors;
		for (int32 i = 0; i < Params.sm_n_uncompactors; ++i)
		{
			w_bkb_uncompactors[i] = new CWBigKmerBinUncompactor<SIZE>(Params, Queues);
			bkb_uncompactors.push_back(thread(std::ref(*w_bkb_uncompactors[i])));
		}

		CWBigKmerBinSorter<SIZE>* w_bkb_sorter = new CWBigKmerBinSorter<SIZE>(Params, Queues, sort_func);
		thread bkb_sorter(std::ref(*w_bkb_sorter));

		CWBigKmerBinWriter* w_bkb_writer = new CWBigKmerBinWriter(Params, Queues);
		thread bkb_writer(std::ref(*w_bkb_writer));

		vector<CWBigKmerBinMerger<SIZE>*> w_bkb_mergers(Params.sm_n_mergers);
		vector<thread> bkb_mergers;
		for (int32 i = 0; i < Params.sm_n_mergers; ++i)
		{
			w_bkb_mergers[i] = new CWBigKmerBinMerger<SIZE>(Params, Queues);
			bkb_mergers.push_back(thread(std::ref(*w_bkb_mergers[i])));
		}

		w_completer->InitStage2(Params, Queues);
		gr2_3.push_back(thread(std::ref(*w_completer), false));
		for (auto& m : bkb_mergers)
			m.join();

		bkb_sorter.join();
		bkb_writer.join();
		for (auto& u : bkb_uncompactors)
			u.join();
		bkb_reader.join();
		delete w_bkb_reader;
		for (auto& u : w_bkb_uncompactors)
			delete u;

		delete w_bkb_sorter;
		delete w_bkb_writer;


		for (auto& m : w_bkb_mergers)
			delete m;

		release_thr_sm = new thread([&] {
			delete Queues.bbpq;
			delete Queues.bbkq;
			delete Queues.sm_cbc;
			delete Queues.bbspq;
		});
	}
	else
	{
		gr2_3.push_back(thread(std::ref(*w_completer), false));
	}

	Queues.pmm_radix_buf->release();
	delete Queues.pmm_radix_buf;

	for (auto p = gr2_3.begin(); p != gr2_3.end(); ++p)
		p->join();


	if (Params.use_strict_mem)
	{
		Queues.sm_pmm_input_file->release();
		Queues.sm_pmm_expand->release();
		Queues.sm_pmm_sort->release();
		Queues.sm_pmm_sorter_suffixes->release();
		Queues.sm_pmm_sorter_lut->release();

		delete Queues.sm_pmm_input_file;
		delete Queues.sm_pmm_expand;
		delete Queues.sm_pmm_sort;
		delete Queues.sm_pmm_sorter_suffixes;
		delete Queues.sm_pmm_sorter_lut;

		Queues.sm_pmm_merger_lut->release();
		Queues.sm_pmm_merger_suff->release();
		Queues.sm_pmm_sub_bin_lut->release();
		Queues.sm_pmm_sub_bin_suff->release();

		delete Queues.sm_pmm_merger_lut;
		delete Queues.sm_pmm_merger_suff;
		delete Queues.sm_pmm_sub_bin_lut;
		delete Queues.sm_pmm_sub_bin_suff;
	}

	// ***** End of Stage 2 *****
	w_completer->GetTotal(results.nUniqueKmers, results.nBelowCutoffMin, results.nAboveCutoffMax, results.nTotalKmers);
	
	uint64 stat_n_plus_x_recs, stat_n_recs, stat_n_recs_tmp, stat_n_plus_x_recs_tmp;
	stat_n_plus_x_recs = stat_n_recs = stat_n_recs_tmp = stat_n_plus_x_recs_tmp = 0;

	thread* release_thr_st2_2 = new thread([&] {

		delete w_reader;
		for (int i = 0; i < Params.n_sorters; ++i)
		{
			w_sorters[i]->GetDebugStats(stat_n_recs_tmp, stat_n_plus_x_recs_tmp);
			stat_n_plus_x_recs += stat_n_plus_x_recs_tmp;
			stat_n_recs += stat_n_recs_tmp;
			delete w_sorters[i];
		}
		delete w_completer;

		delete Queues.sorters_manager;

		delete Queues.input_files_queue;
		delete Queues.bq;
		delete Queues.part_queue;
		delete Queues.bpq;
		delete Queues.kq;
		delete Queues.tlbq;
	});

	Queues.bd->reset_reading();
	while ((bin_id = Queues.bd->get_next_bin()) >= 0)
	{
		Queues.bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, n_super_kmers);
		delete file;
	}

	delete Queues.bd;
	delete Queues.epd;

	results.tmpSizeStrictMemory = 0;
	if (!Params.use_strict_mem)
	{
		release_thr_st2_1->join();
	}
	else
	{
		release_thr_sm->join();
		Queues.bbd->reset_reading();
		int32 sub_bin_id = 0;
		uint32 lut_prefix_len = 0;
		uint64 n_kmers = 0;
		uint64 file_size = 0;
		uint32 size = 0;
		FILE* file = nullptr;
		while (Queues.bbd->next_bin(bin_id, size))
		{
			while (Queues.bbd->next_sub_bin(bin_id, sub_bin_id, lut_prefix_len, n_kmers, file, name, file_size))
			{
				results.tmpSizeStrictMemory += file_size;
			}
		}
		delete Queues.bbd;
		delete release_thr_sm;
	}
	release_thr_st2_2->join();

	delete release_thr_st2_1;
	delete release_thr_st2_2;
	delete Queues.s_mapper;
	results.maxDiskUsage = Queues.disk_logger->get_max();

	delete Queues.disk_logger;
	if (!Params.use_strict_mem)
	{
		w2.stopTimer();
		results.time = w2.getElapsedTime();
	}
	else
	{
		w3.stopTimer();
		results.timeStrictMem = w3.getElapsedTime();
	}

	return results;
}


#endif

// ***** EOF
