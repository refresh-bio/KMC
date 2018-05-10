/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
*/

#ifndef _FASTQ_FILTER_H
#define _FASTQ_FILTER_H

#include "config.h"
#include "queues.h"

#include "../kmc_api/kmc_file.h"

//************************************************************************************************************
// CFastqFilter - filter of reads
//************************************************************************************************************
class CFastqFilter
{
private:
	CFilteringParams::FilterMode mode;
	CPartQueue *input_part_queue, *filtered_part_queue;
	CMemoryPool *pmm_fastq_reader;
	CMemoryPool *pmm_fastq_filter;
	CFilteringParams::file_type input_file_type, output_file_type;
	CKMCFile& kmc_api;
	uint64 output_part_size;

	uchar* input_part;
	uint64 input_part_size;
	uint64 input_part_pos;
	uchar* output_part;
	uint64 output_part_pos;

	std::vector<uint32> counters;
	std::string read;
	struct {
		uint64 read_header_start;
		uint64 read_header_end;
		uint64 read_start;
		uint64 read_end;
		uint64 quality_header_start;
		uint64 quality_header_end;
		uint64 quality_start;
		uint64 quality_end;
		uint64 end;
	}seq_desc;

	uint32 trim_len; //for trim mode

	bool use_float_value;
	float f_max_kmers;
	float f_min_kmers;
	uint32 n_max_kmers;
	uint32 n_min_kmers;
	uint32 kmer_len;

	template<class Helper> void ProcessImpl();



	bool NextSeqFastq();
	bool NextSeqFasta();
	bool FilterRead();
	bool FilterReadTrim();
	void HardMask();
public:
	CFastqFilter(CFilteringParams& Params, CFilteringQueues& Queues, CKMCFile& kmc_api);
	void Process();


private: //Helpers classes for ProcessImpl

	class FastqToFastqHelper;
	class FastqToFastaHelper;
	class FastaToFastaHelper;

	class TrimFastqToFastqHelper;
	class TrimFastqToFastaHelper;
	class TrimFastaToFastaHelper;

	class HardMaskFastqToFastqHelper;
	class HardMaskFastqToFastaHelper;
	class HardMaskFastaToFastaHelper;

};

//************************************************************************************************************
// CWFastqFilter - wrapper for CFastqFilter class - for multithreading purposes
//************************************************************************************************************
class CWFastqFilter
{
	std::unique_ptr<CFastqFilter> ff;
public:
	CWFastqFilter(CFilteringParams& Params, CFilteringQueues& Queues, CKMCFile& kmc_api);
	void operator()();
};


#endif

// ***** EOF
