/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
*/

#include "stdafx.h"
#include "fastq_filter.h"
#include <numeric>

using namespace std;

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
CFastqFilter::CFastqFilter(CFilteringParams& Params, CFilteringQueues& Queues, CKMCFile& kmc_api) :
kmc_api(kmc_api)
{
	input_part_queue = Queues.input_part_queue;
	filtered_part_queue = Queues.filtered_part_queue;
	pmm_fastq_reader = Queues.pmm_fastq_reader;
	pmm_fastq_filter = Queues.pmm_fastq_filter;
	input_file_type = Params.input_file_type;
	output_file_type = Params.output_file_type;
	use_float_value = Params.use_float_value;
	f_max_kmers = Params.f_max_kmers;
	f_min_kmers = Params.f_min_kmers;
	n_max_kmers = Params.n_max_kmers;
	n_min_kmers = Params.n_min_kmers;
	kmer_len = Params.kmer_len;
	output_part_size = Params.mem_part_pmm_fastq_reader;
	mode = Params.filter_mode;
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
	if (mode == CFilteringParams::FilterMode::trim)
	{
		if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fastq)
			ProcessImpl<TrimFastqToFastqHelper>();
		else if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fasta)
			ProcessImpl<TrimFastqToFastaHelper>();
		else if (input_file_type == CFilteringParams::file_type::fasta && output_file_type == CFilteringParams::file_type::fasta)
			ProcessImpl<TrimFastaToFastaHelper>();
		else
		{
			cerr << "Error: this file type is not supported by filter operation\n";
			exit(1);
		}
	}
	else if (mode == CFilteringParams::FilterMode::hard_mask)
	{
		if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fastq)
			ProcessImpl<HardMaskFastqToFastqHelper>();
		else if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fasta)
			ProcessImpl<HardMaskFastqToFastaHelper>();
		else if (input_file_type == CFilteringParams::file_type::fasta && output_file_type == CFilteringParams::file_type::fasta)
			ProcessImpl<HardMaskFastaToFastaHelper>();
		else
			cerr << "Error: this file type is not supported by filter operation\n";
	}
	else
	{
		if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fastq)
			ProcessImpl<FastqToFastqHelper>();
		else if (input_file_type == CFilteringParams::file_type::fastq && output_file_type == CFilteringParams::file_type::fasta)
			ProcessImpl<FastqToFastaHelper>();
		else if (input_file_type == CFilteringParams::file_type::fasta && output_file_type == CFilteringParams::file_type::fasta)
			ProcessImpl<FastaToFastaHelper>();
		else
		{
			cerr << "Error: this file type is not supported by filter operation\n";
			exit(1);
		}
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
	for (auto counter : counters)
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
bool CFastqFilter::FilterReadTrim()
{
	uint32 read_len = static_cast<uint32>(seq_desc.read_end - seq_desc.read_start);
	read.assign((char*)input_part + seq_desc.read_start, read_len);

	kmc_api.GetCountersForRead(read, counters);
	if (counters[0] < n_min_kmers)
		return false;

	trim_len = kmer_len;
	for (uint32 i = 1; i < counters.size(); ++i, ++trim_len)
	{
		if (counters[i] < n_min_kmers)
			break;
	}


	return true;
}
void CFastqFilter::HardMask()
{
	uint32 read_len = static_cast<uint32>(seq_desc.read_end - seq_desc.read_start);
	read.assign((char*)input_part + seq_desc.read_start, read_len);
	kmc_api.GetCountersForRead(read, counters);
	uchar* read_in = input_part + seq_desc.read_start;
	uint32 read_in_pos = 0;
	for (uint32 counter_pos = 0; counter_pos < counters.size(); ++counter_pos)
	{
		if (counters[counter_pos] < n_min_kmers)
		{
			while (read_in_pos < counter_pos + kmer_len)
			{
				output_part[output_part_pos++] = 'N';
				read_in_pos++;
			}
		}
		else if (read_in_pos <= counter_pos)
			output_part[output_part_pos++] = read_in[read_in_pos++];
	}
	while (read_in_pos < read_len)
		output_part[output_part_pos++] = read_in[read_in_pos++];
	output_part[output_part_pos++] = '\n';
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
template<class Helper>
void CFastqFilter::ProcessImpl()
{
	Helper helper(*this);
	pmm_fastq_filter->reserve(output_part);
	output_part_pos = 0;
	uint64 required_size;
	while (input_part_queue->pop(input_part, input_part_size))
	{
		input_part_pos = 0;
		while (helper.NextSeq())
		{
			if (helper.FilterRead())
			{
				required_size = helper.GetReqSize();
				if (output_part_pos + required_size > output_part_size)
				{
					filtered_part_queue->push(output_part, output_part_pos);
					pmm_fastq_filter->reserve(output_part);
					output_part_pos = 0;
				}
				helper.SendToOutBuf();

			}
		}
		pmm_fastq_reader->free(input_part);
	}
	filtered_part_queue->push(output_part, output_part_pos);
	filtered_part_queue->mark_completed();
}


/*****************************************************************************************************************************/
/************************************************************ HELPERS CLASSES ************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
class CFastqFilter::FastqToFastqHelper
{
	CFastqFilter& owner;
public:
	FastqToFastqHelper(CFastqFilter& owner) :owner(owner){}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.quality_header_start - owner.seq_desc.read_header_start + 1 + owner.seq_desc.end - owner.seq_desc.quality_header_end;
	}
	void SendToOutBuf() const
	{
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.quality_header_start - owner.seq_desc.read_header_start + 1);
		owner.output_part_pos += owner.seq_desc.quality_header_start - owner.seq_desc.read_header_start + 1;
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.quality_header_end, owner.seq_desc.end - owner.seq_desc.quality_header_end);
		owner.output_part_pos += owner.seq_desc.end - owner.seq_desc.quality_header_end;
	}

	bool NextSeq() const
	{
		return owner.NextSeqFastq();
	}

	bool FilterRead() const
	{
		return owner.FilterRead();
	}
};

/*****************************************************************************************************************************/
class CFastqFilter::FastqToFastaHelper
{
	CFastqFilter& owner;
public:
	FastqToFastaHelper(CFastqFilter& owner) :owner(owner){}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.quality_header_start - owner.seq_desc.read_header_start;
	}
	void SendToOutBuf() const
	{
		owner.input_part[owner.seq_desc.read_header_start] = '>';
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.quality_header_start - owner.seq_desc.read_header_start);
		owner.output_part_pos += owner.seq_desc.quality_header_start - owner.seq_desc.read_header_start;
	}
	bool NextSeq() const
	{
		return owner.NextSeqFastq();
	}

	bool FilterRead() const
	{
		return owner.FilterRead();
	}
};

/*****************************************************************************************************************************/
class CFastqFilter::FastaToFastaHelper
{
	CFastqFilter& owner;
public:
	FastaToFastaHelper(CFastqFilter& owner) :owner(owner){}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.end - owner.seq_desc.read_header_start;
	}

	void SendToOutBuf() const
	{
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.end - owner.seq_desc.read_header_start);
		owner.output_part_pos += owner.seq_desc.end - owner.seq_desc.read_header_start;
	}
	bool NextSeq() const
	{
		return owner.NextSeqFasta();
	}

	bool FilterRead() const
	{
		return owner.FilterRead();
	}
};

/*****************************************************************************************************************************/
class CFastqFilter::TrimFastqToFastqHelper
{
	CFastqFilter& owner;
public:
	TrimFastqToFastqHelper(CFastqFilter& owner) :owner(owner){}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.read_header_end - owner.seq_desc.read_header_start + 1 + owner.trim_len + 3 + owner.trim_len + 1;
	}
	void SendToOutBuf() const
	{
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.read_header_end - owner.seq_desc.read_header_start);
		owner.output_part_pos += owner.seq_desc.read_header_end - owner.seq_desc.read_header_start;
		owner.output_part[owner.output_part_pos++] = '\n';
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_start, owner.trim_len);
		owner.output_part_pos += owner.trim_len;
		owner.output_part[owner.output_part_pos++] = '\n';
		owner.output_part[owner.output_part_pos++] = '+';
		owner.output_part[owner.output_part_pos++] = '\n';
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.quality_start, owner.trim_len);
		owner.output_part_pos += owner.trim_len;
		owner.output_part[owner.output_part_pos++] = '\n';
	}

	bool NextSeq() const
	{
		return owner.NextSeqFastq();
	}

	bool FilterRead() const
	{
		return owner.FilterReadTrim();
	}
};


/*****************************************************************************************************************************/
class CFastqFilter::TrimFastqToFastaHelper
{
	CFastqFilter& owner;
public:
	TrimFastqToFastaHelper(CFastqFilter& owner) :owner(owner){}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.read_header_end - owner.seq_desc.read_header_start + 1 + owner.trim_len + 1;
	}
	void SendToOutBuf() const
	{
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.read_header_end - owner.seq_desc.read_header_start);
		owner.output_part[owner.output_part_pos] = '>';
		owner.output_part_pos += owner.seq_desc.read_header_end - owner.seq_desc.read_header_start;
		owner.output_part[owner.output_part_pos++] = '\n';
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_start, owner.trim_len);
		owner.output_part_pos += owner.trim_len;
		owner.output_part[owner.output_part_pos++] = '\n';
	}

	bool NextSeq() const
	{
		return owner.NextSeqFastq();
	}

	bool FilterRead() const
	{
		return owner.FilterReadTrim();
	}
};



/*****************************************************************************************************************************/
class CFastqFilter::TrimFastaToFastaHelper
{
	CFastqFilter& owner;
public:
	TrimFastaToFastaHelper(CFastqFilter& owner) :owner(owner){}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.read_header_end - owner.seq_desc.read_header_start + 1 + owner.trim_len + 1;
	}
	void SendToOutBuf() const
	{
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.read_header_end - owner.seq_desc.read_header_start);
		owner.output_part_pos += owner.seq_desc.read_header_end - owner.seq_desc.read_header_start;
		owner.output_part[owner.output_part_pos++] = '\n';
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_start, owner.trim_len);
		owner.output_part_pos += owner.trim_len;
		owner.output_part[owner.output_part_pos++] = '\n';
	}

	bool NextSeq() const
	{
		return owner.NextSeqFasta();
	}

	bool FilterRead() const
	{
		return owner.FilterReadTrim();
	}
};

class CFastqFilter::HardMaskFastqToFastqHelper
{
	CFastqFilter& owner;
public:
	HardMaskFastqToFastqHelper(CFastqFilter& owner) : owner(owner) {}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.read_header_end - owner.seq_desc.read_header_start + 1 +
			owner.seq_desc.read_end - owner.seq_desc.read_start + 1 +
			2 +
			owner.seq_desc.quality_end - owner.seq_desc.quality_start + 1;
	}
	void SendToOutBuf() const
	{
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.read_header_end - owner.seq_desc.read_header_start);
		owner.output_part_pos += owner.seq_desc.read_header_end - owner.seq_desc.read_header_start;
		owner.output_part[owner.output_part_pos++] = '\n';
		owner.HardMask();
		owner.output_part[owner.output_part_pos++] = '+';
		owner.output_part[owner.output_part_pos++] = '\n';
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.quality_start, owner.seq_desc.quality_end - owner.seq_desc.quality_start);
		owner.output_part_pos += owner.seq_desc.quality_end - owner.seq_desc.quality_start;
		owner.output_part[owner.output_part_pos++] = '\n';		
	}
	bool NextSeq() const
	{
		return owner.NextSeqFastq();
	}
	bool FilterRead() const
	{
		return true;
	}
};
class CFastqFilter::HardMaskFastqToFastaHelper
{
	CFastqFilter& owner;
public:
	HardMaskFastqToFastaHelper(CFastqFilter& owner) : owner(owner) {}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.read_header_end - owner.seq_desc.read_header_start + 1 +
			owner.seq_desc.read_end - owner.seq_desc.read_start + 1;
	}
	void SendToOutBuf() const
	{
		owner.input_part[owner.seq_desc.read_header_start] = '>';
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.read_header_end - owner.seq_desc.read_header_start);
		owner.output_part_pos += owner.seq_desc.read_header_end - owner.seq_desc.read_header_start;
		owner.output_part[owner.output_part_pos++] = '\n';
		owner.HardMask();		
	}
	bool NextSeq() const
	{
		return owner.NextSeqFastq();
	}
	bool FilterRead() const
	{
		return true;
	}
};
class CFastqFilter::HardMaskFastaToFastaHelper
{
	CFastqFilter& owner;
public:
	HardMaskFastaToFastaHelper(CFastqFilter& owner) : owner(owner) {}
	uint64 GetReqSize() const
	{
		return owner.seq_desc.read_header_end - owner.seq_desc.read_header_start + 1 +
			owner.seq_desc.read_end - owner.seq_desc.read_start + 1;
	}
	void SendToOutBuf() const
	{
		memcpy(owner.output_part + owner.output_part_pos, owner.input_part + owner.seq_desc.read_header_start, owner.seq_desc.read_header_end - owner.seq_desc.read_header_start);
		owner.output_part_pos += owner.seq_desc.read_header_end - owner.seq_desc.read_header_start;
		owner.output_part[owner.output_part_pos++] = '\n';
		owner.HardMask();	
	}
	bool NextSeq() const
	{
		return owner.NextSeqFasta();
	}
	bool FilterRead() const
	{
		return true;
	}
};
// ***** EOF
