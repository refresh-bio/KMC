/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#include "stdafx.h"
#include "fastq_filter.h"
#include "asmlib_wrapper.h"
#include <numeric>

using namespace std;

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
CFastqFilter::CFastqFilter(CFilteringParams& Params, CFilteringQueues& Queues, CKMCFile& kmc_api) :
	kmc_api(kmc_api)
{
	input_part_queue      =   Queues.input_part_queue;
	filtered_part_queue   =   Queues.filtered_part_queue;
	pmm_fastq_reader      =   Queues.pmm_fastq_reader;	
	pmm_fastq_filter      =   Queues.pmm_fastq_filter;
	input_file_type       =   Params.input_file_type;
	output_file_type	  =   Params.output_file_type;
	use_float_value		  =   Params.use_float_value;
	f_max_kmers			  =   Params.f_max_kmers;
	f_min_kmers			  =   Params.f_min_kmers;
	n_max_kmers			  =   Params.n_max_kmers;
	n_min_kmers			  =   Params.n_min_kmers;
	kmer_len			  =   Params.kmer_len;
	output_part_size      =   Params.mem_part_pmm_fastq_reader;
}

/*****************************************************************************************************************************/
CWFastqFilter::CWFastqFilter(CFilteringParams& Params, CFilteringQueues& Queues, CKMCFile& kmc_api)
{
	ff = make_unique<CFastqFilter>(Params, Queues, kmc_api);
}


/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
void CFastqFilter::Process()
{
	if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fastq)
		ProcessFastqToFastq();
	else if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fasta)
		ProcessFastqToFasta();
	else if (input_file_type == CFilteringParams::file_type::fasta && output_file_type == CFilteringParams::file_type::fasta)
		ProcessFastaToFasta();
	else
	{
		cout << "Error: this file type is not supported by filter operation\n";
		exit(1);
	}
}

/*****************************************************************************************************************************/
void CWFastqFilter::operator()()
{
	ff->Process();
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
bool CFastqFilter::FilterRead() 
{
	uint32 read_len = static_cast<uint32>(seq_desc.read_end - seq_desc.read_start);	
	read.assign((char*)input_part + seq_desc.read_start, read_len);

	kmc_api.GetCountersForRead(read, counters);
	uint32 valid_kmers = 0;
	for(auto counter : counters)
		if (counter)
			++valid_kmers;

	if (use_float_value)
	{
		uint32 min = static_cast<uint32>(f_min_kmers * (read_len - kmer_len + 1));
		uint32 max = static_cast<uint32>(f_max_kmers * (read_len - kmer_len + 1));
		if (valid_kmers >= min && valid_kmers <= max)
			return true;
		return false;
	}
	else
	{
		if (valid_kmers >= n_min_kmers && valid_kmers <= n_max_kmers)
			return true;
		return false;
	}
}

/*****************************************************************************************************************************/
bool CFastqFilter::NextSeqFasta()
{
	// Title
	char c;
	if (input_part_pos >= input_part_size)
		return false;
	c = input_part[input_part_pos++];
	if (c != '>')
		return false;

	seq_desc.read_header_start = input_part_pos - 1;

	for (; input_part_pos < input_part_size;)
	{
		c = input_part[input_part_pos++];
		if (c < 32)					// newliners
			break;
	}
	seq_desc.read_header_end = input_part_pos - 1;

	if (input_part_pos >= input_part_size)
		return false;

	c = input_part[input_part_pos++];
	if (c >= 32)
		input_part_pos--;
	else if (input_part_pos >= input_part_size)
		return false;

	seq_desc.read_start = input_part_pos;
	// Sequence
	for (; input_part_pos < input_part_size;)
	{
		c = input_part[input_part_pos++];
		if (c < 32)					// newliners
			break;		
	}
	seq_desc.read_end = input_part_pos - 1;

	seq_desc.end = input_part_pos;

	if (input_part_pos >= input_part_size)
		return true;

	seq_desc.end++;
	if (input_part[input_part_pos++] >= 32)
	{
		input_part_pos--;
		seq_desc.end--;
	}

	else if (input_part_pos >= input_part_size)
		return true;

	return (c == '\n' || c == '\r');
}

/*****************************************************************************************************************************/
bool CFastqFilter::NextSeqFastq()
{
	char c;
	// Title
	if (input_part_pos >= input_part_size)
		return false;

	
	c = input_part[input_part_pos++];
	if (c != '@')
		return false;

	seq_desc.read_header_start = input_part_pos - 1;

	for (; input_part_pos < input_part_size;)
	{
		c = input_part[input_part_pos++];
		if (c < 32)					// newliners
			break;
	}
	seq_desc.read_header_end = input_part_pos - 1;

	if (input_part_pos >= input_part_size)
		return false;

	c = input_part[input_part_pos++];
	if (c >= 32)
		input_part_pos--;
	else if (input_part_pos >= input_part_size)
		return false;

	seq_desc.read_start = input_part_pos;
	// Sequence
	for (; input_part_pos < input_part_size;)
	{
		c = input_part[input_part_pos++];
		if (c < 32)					// newliners
			break;		
	}
	seq_desc.read_end = input_part_pos - 1;

	if (input_part_pos >= input_part_size)
		return false;

	c = input_part[input_part_pos++];
	if (c >= 32)
		input_part_pos--;
	else if (input_part_pos >= input_part_size)
		return false;

	// Plus	
	c = input_part[input_part_pos++];
	if (input_part_pos >= input_part_size)
		return false;
	if (c != '+')
		return false;

	seq_desc.quality_header_start = input_part_pos - 1;

	for (; input_part_pos < input_part_size;)
	{
		c = input_part[input_part_pos++];
		if (c < 32)					// newliners
			break;
	}
	seq_desc.quality_header_end = input_part_pos - 1;

	if (input_part_pos >= input_part_size)
		return false;

	c = input_part[input_part_pos++];
	if (c >= 32)
		input_part_pos--;
	else if (input_part_pos >= input_part_size)
		return false;

	// Quality
	seq_desc.quality_start = input_part_pos;

	input_part_pos += seq_desc.read_end - seq_desc.read_start;
	if (input_part_pos >= input_part_size)
		return false;
	c = input_part[input_part_pos++];

	seq_desc.quality_end = input_part_pos - 1;

	seq_desc.end = input_part_pos;

	if (input_part_pos >= input_part_size)
		return true;

	seq_desc.end++;
	if (input_part[input_part_pos++] >= 32)
	{
		input_part_pos--;
		seq_desc.end--;

	}
	else if (input_part_pos >= input_part_size)
		return true;

	return c == '\n' || c == '\r';
}

/*****************************************************************************************************************************/
void CFastqFilter::ProcessFastaToFasta()
{
	pmm_fastq_filter->reserve(output_part);
	output_part_pos = 0;
	uint64 required_size;
	while (input_part_queue->pop(input_part, input_part_size))
	{		
		input_part_pos = 0;
		while (NextSeqFasta())
		{
			if (FilterRead())
			{
				required_size = seq_desc.end - seq_desc.read_header_start;
				if (output_part_pos + required_size > output_part_size)
				{
					filtered_part_queue->push(output_part, output_part_pos);
					pmm_fastq_filter->reserve(output_part);
					output_part_pos = 0;
				}
				A_memcpy(output_part + output_part_pos, input_part + seq_desc.read_header_start, required_size);
				output_part_pos += required_size;
			}
		}
		pmm_fastq_reader->free(input_part);
	}
	filtered_part_queue->push(output_part, output_part_pos);
	filtered_part_queue->mark_completed();
}

/*****************************************************************************************************************************/
void CFastqFilter::ProcessFastqToFastq()
{
	pmm_fastq_filter->reserve(output_part);
	output_part_pos = 0;
	uint64 required_size;
	while (input_part_queue->pop(input_part, input_part_size))
	{
		input_part_pos = 0;
		while (NextSeqFastq())
		{
			if (FilterRead())
			{
				required_size = seq_desc.quality_header_start - seq_desc.read_header_start + 1 + seq_desc.end - seq_desc.quality_header_end;
				if (output_part_pos + required_size > output_part_size)
				{
					filtered_part_queue->push(output_part, output_part_pos);
					pmm_fastq_filter->reserve(output_part);
					output_part_pos = 0;
				}
				A_memcpy(output_part + output_part_pos, input_part + seq_desc.read_header_start, seq_desc.quality_header_start - seq_desc.read_header_start + 1); 
				output_part_pos += seq_desc.quality_header_start - seq_desc.read_header_start + 1;
				A_memcpy(output_part + output_part_pos, input_part + seq_desc.quality_header_end, seq_desc.end - seq_desc.quality_header_end); 
				output_part_pos += seq_desc.end - seq_desc.quality_header_end;
			}
		}
		pmm_fastq_reader->free(input_part);
	}
	filtered_part_queue->push(output_part, output_part_pos);
	filtered_part_queue->mark_completed();
}

/*****************************************************************************************************************************/
void CFastqFilter::ProcessFastqToFasta()
{
	pmm_fastq_filter->reserve(output_part);
	output_part_pos = 0;
	uint64 required_size;
	while (input_part_queue->pop(input_part, input_part_size))
	{
		input_part_pos = 0;
		while (NextSeqFastq())
		{
			if (FilterRead())
			{
				required_size = seq_desc.quality_header_start - seq_desc.read_header_start;
				if (output_part_pos + required_size > output_part_size)
				{
					filtered_part_queue->push(output_part, output_part_pos);
					pmm_fastq_filter->reserve(output_part);
					output_part_pos = 0;
				}
				input_part[seq_desc.read_header_start] = '>';
				A_memcpy(output_part + output_part_pos, input_part + seq_desc.read_header_start, seq_desc.quality_header_start - seq_desc.read_header_start);
				output_part_pos += seq_desc.quality_header_start - seq_desc.read_header_start;
			}
		}
		pmm_fastq_reader->free(input_part);
	}
	filtered_part_queue->push(output_part, output_part_pos);
	filtered_part_queue->mark_completed();
}



// ***** EOF