/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#ifndef _SPLITTER_H
#define _SPLITTER_H

#include "defs.h"
#include "kmer.h"
#include "kb_storer.h"
#include "kb_collector.h"
#include "prob_qual.h"
#include "kb_reader.h"
#include "kb_sorter.h"
#include "kb_completer.h"
#include "queues.h"
#include "s_mapper.h"
#include "mmer.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "small_k_buf.h"

using namespace std;

//************************************************************************************************************
//************************************************************************************************************
template <bool QUAKE_MODE> class CSplitter_Impl;

//************************************************************************************************************
// CSplitter class - splits kmers into bins according to their prefix
//************************************************************************************************************
template <bool QUAKE_MODE> class CSplitter {
	CMemoryMonitor *mm;
	uint64 total_kmers = 0;
	//CExKmer ex_kmer;
	uchar *part;
	uint64 part_size, part_pos;
	CKmerBinCollector **bins;
	CBinPartQueue *bin_part_queue;
	CBinDesc *bd;
	CMemoryPool *pmm_reads;
	int64 mem_part_pmm_bins;
	int64 mem_part_pmm_reads;

	char codes[256];
	bool use_quake;
	input_type file_type;
	int lowest_quality;
	bool both_strands;

	uint32 kmer_len;
	//uint32 prefix_len;
	uint32 signature_len;
	uint32 n_bins;	
	uint64 n_reads;//for multifasta its a sequences counter	

	CSignatureMapper* s_mapper;

	inline bool GetSeq(char *seq, uint32 &seq_size);
	inline bool GetSeq(char *seq, char *quals, uint32 &seq_size);

	friend class CSplitter_Impl<QUAKE_MODE>;

public:
	inline void CalcStats(uchar* _part, uint64 _part_size, uint32* _stats);

	static uint32 MAX_LINE_SIZE;

	CSplitter(CKMCParams &Params, CKMCQueues &Queues); 
	void InitBins(CKMCParams &Params, CKMCQueues &Queues);
	
	bool ProcessReads(uchar *_part, uint64 _part_size);
	
	template<typename COUNTER_TYPE>
	bool ProcessReadsSmallK(uchar *_part, uint64 _part_size, CSmallKBuf<COUNTER_TYPE>& small_k_buf);

	void Complete();

	void GetTotal(uint64 &_n_reads);

	uint64 GetTotalKmers();

	~CSplitter();
};

template <bool QUAKE_MODE> uint32 CSplitter<QUAKE_MODE>::MAX_LINE_SIZE = 1 << 14;


//************************************************************************************************************
// Implementation of ProcessReads and Complete methods for various types and sizes of kmer class
//************************************************************************************************************
template <bool QUAKE_MODE> class CSplitter_Impl {
public: 
	static bool ProcessReads(CSplitter<QUAKE_MODE> &ptr, uchar *_part, uint64 _part_size);
	template<typename COUNTER_TYPE>
	static bool ProcessReadsSmallK(CSplitter<QUAKE_MODE> &ptr, uchar *_part, uint64 _part_size, CSmallKBuf<COUNTER_TYPE>& small_k_buf);
};

template <> class CSplitter_Impl<false> {
public: 	
	static bool ProcessReads(CSplitter<false> &ptr, uchar *_part, uint64 _part_size);
	template<typename COUNTER_TYPE>
	static bool ProcessReadsSmallK(CSplitter<false> &ptr, uchar *_part, uint64 _part_size, CSmallKBuf<COUNTER_TYPE>& small_k_buf);
};

template <> class CSplitter_Impl<true> {
public: 
	static bool ProcessReads(CSplitter<true> &ptr, uchar *_part, uint64 _part_size);		
	static bool ProcessReadsSmallK(CSplitter<true> &ptr, uchar *_part, uint64 _part_size, CSmallKBuf<float>& small_k_buf);
};

//----------------------------------------------------------------------------------
// Return a single record from FASTA/FASTQ data
template <bool QUAKE_MODE> bool CSplitter<QUAKE_MODE>::GetSeq(char *seq, uint32 &seq_size)
{
	uchar c = 0;
	uint32 pos = 0;
	
	if(file_type == fasta)
	{
		// Title
		if(part_pos >= part_size)
			return false;
		c = part[part_pos++];
		if(c != '>')
			return false;
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Sequence
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
			seq[pos++] = codes[c];
		}
		seq_size = pos;

		if(part_pos >= part_size)
			return true;

		if(part[part_pos++] >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return true;
	}
	else if(file_type == fastq)
	{
		// Title
		if(part_pos >= part_size)
			return false;
		c = part[part_pos++];
		if(c != '@')
			return false;
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Sequence
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
			seq[pos++] = codes[c];
		}
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Plus
		c = part[part_pos++];
		if(part_pos >= part_size)
			return false;
		if(c != '+')
			return false;
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Quality
		part_pos += pos;
		if(part_pos >= part_size)
			return false;
		c = part[part_pos++];
		seq_size = pos;

		if(part_pos >= part_size)
			return true;

		if(part[part_pos++] >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return true;
	}
	else if(file_type == multiline_fasta)
	{
		if(part_pos >= part_size)
			return false;
		if(part[part_pos] == '>')//need to ommit header
		{
			++n_reads;
			for(;part_pos < part_size && part[part_pos] != '\n' && part[part_pos] != '\r';++part_pos);//find EOF
			++part_pos;
			if(part[part_pos] == '\n' || part[part_pos] == '\r')
				++part_pos;
		}
		for(;part_pos < part_size && pos < mem_part_pmm_reads && part[part_pos] != '>';)
		{
			seq[pos++] = codes[part[part_pos++]];
		}
		seq_size = pos;
		if(part_pos < part_size && part[part_pos] != '>')//need to copy last k-1 kmers 
		{
			part_pos -= kmer_len - 1;
		}
		return true;
		
	}

	return (c == '\n' || c == '\r');
}

//----------------------------------------------------------------------------------
// Return a single record with quality codes from FASTA/FASTQ data
template <bool QUAKE_MODE> bool CSplitter<QUAKE_MODE>::GetSeq(char *seq, char *quals, uint32 &seq_size)
{
	uchar c;
	uint32 pos = 0;

	if(file_type == fasta || file_type == multiline_fasta)
	{
		return false;		// FASTA file does not store quality values
	}
	else
	{
		// Title
		if(part_pos >= part_size)
			return false;
		c = part[part_pos++];
		if(c != '@')
			return false;
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Sequence
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
			seq[pos++] = codes[c];
		}
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Plus
		c = part[part_pos++];
		if(part_pos >= part_size)
			return false;
		if(c != '+')
			return false;
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Quality
		copy(part+part_pos, part+part_pos+pos, quals);

		part_pos += pos;
		if(part_pos >= part_size)
			return false;
		c = part[part_pos++];
		seq_size = pos;

		if(part_pos >= part_size)
			return true;

		if(part[part_pos++] >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return true;
	}

	return (c == '\n' || c == '\r');
}


template <bool QUAKE_MODE> void CSplitter<QUAKE_MODE>::CalcStats(uchar* _part, uint64 _part_size, uint32* _stats)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char *seq;
	uint32 seq_size;
	pmm_reads->reserve(seq);

	uint32 signature_start_pos;
	CMmer current_signature(signature_len), end_mmer(signature_len);

	uint32 i;
	uint32 len;//length of extended kmer

	while (GetSeq(seq, seq_size))
	{
		i = 0;
		len = 0;
		while (i + kmer_len - 1 < seq_size)
		{
			bool contains_N = false;
			//building first signature after 'N' or at the read begining
			for (uint32 j = 0; j < signature_len; ++j, ++i)
			if (seq[i] < 0)//'N'
			{
				contains_N = true;
				break;
			}
			//signature must be shorter than k-mer so if signature contains 'N', k-mer will contains it also
			if (contains_N)
			{
				++i;
				continue;
			}
			len = signature_len;
			signature_start_pos = i - signature_len;
			current_signature.insert(seq + signature_start_pos);
			end_mmer.set(current_signature);
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)//'N'
				{
					if (len >= kmer_len)
						_stats[current_signature.get()] += 1 + len - kmer_len;
					len = 0;
					++i;
					break;
				}
				end_mmer.insert(seq[i]);
				if (end_mmer < current_signature)//signature at the end of current k-mer is lower than current
				{
					if (len >= kmer_len)
					{
						_stats[current_signature.get()] += 1 + len - kmer_len;
						len = kmer_len - 1;
					}
					current_signature.set(end_mmer);
					signature_start_pos = i - signature_len + 1;
				}
				else if (end_mmer == current_signature)
				{
					current_signature.set(end_mmer);
					signature_start_pos = i - signature_len + 1;
				}
				else if (signature_start_pos + kmer_len - 1 < i)//need to find new signature
				{
					_stats[current_signature.get()] += 1 + len - kmer_len;
					len = kmer_len - 1;
					//looking for new signature
					++signature_start_pos;
					//building first signature in current k-mer
					end_mmer.insert(seq + signature_start_pos);
					current_signature.set(end_mmer);
					for (uint32 j = signature_start_pos + signature_len; j <= i; ++j)
					{
						end_mmer.insert(seq[j]);
						if (end_mmer <= current_signature)
						{
							current_signature.set(end_mmer);
							signature_start_pos = j - signature_len + 1;
						}
					}
				}
				++len;
			}
		}
		if (len >= kmer_len)//last one in read
			_stats[current_signature.get()] += 1 + len - kmer_len;
	}

	cerr << '*';
	cerr.flush();

	pmm_reads->free(seq);
}
//----------------------------------------------------------------------------------
// Assigns queues and monitors
template <bool QUAKE_MODE> CSplitter<QUAKE_MODE>::CSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	mm             = Queues.mm;
	file_type       = Params.file_type;
	use_quake	   = Params.use_quake;
	lowest_quality = Params.lowest_quality;
	both_strands   = Params.both_strands;	

	bin_part_queue = Queues.bpq;
	bd			   = Queues.bd;
	pmm_reads	   = Queues.pmm_reads;
	kmer_len		  = Params.kmer_len;
	signature_len     = Params.signature_len;
	
	mem_part_pmm_bins = Params.mem_part_pmm_bins;

	mem_part_pmm_reads = Params.mem_part_pmm_reads;

	s_mapper = Queues.s_mapper;

	part = NULL;

	// Prepare encoding of symbols
	for(int i = 0; i < 256; ++i)
		codes[i] = -1;
	codes['A'] = codes['a'] = 0;
	codes['C'] = codes['c'] = 1;
	codes['G'] = codes['g'] = 2;
	codes['T'] = codes['t'] = 3;

	n_reads = 0;
	bins = NULL;
}


template <bool QUAKE_MODE> void CSplitter<QUAKE_MODE>::InitBins(CKMCParams &Params, CKMCQueues &Queues)
{
	n_bins = Params.n_bins;
	uint32 buffer_size = Params.bin_part_size;
	// Create objects for all bin
	bins = new CKmerBinCollector*[n_bins];
	for (uint32 i = 0; i < n_bins; ++i)
	{
		bins[i] = new CKmerBinCollector(Queues, Params, buffer_size, i);
		bd->insert(i, NULL, "", 0, 0, 0, 0, buffer_size, kmer_len);
	}
}
//----------------------------------------------------------------------------------
// Release memory
template <bool QUAKE_MODE> CSplitter<QUAKE_MODE>::~CSplitter()
{
	if (bins)
	{
		for (uint32 i = 0; i < n_bins; ++i)
		if (bins[i])
			delete bins[i];
		delete[] bins;
	}
}

//----------------------------------------------------------------------------------
// Finish the processing of input file
template <bool QUAKE_MODE> void CSplitter<QUAKE_MODE>::Complete()
{
	if (bins)
		for (uint32 i = 0; i < n_bins; ++i)
			if(bins[i])
				bins[i]->Flush();
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part in small k optimization mode
template<bool QUAKE_MODE>
template<typename COUNTER_TYPE>
bool CSplitter<QUAKE_MODE>::ProcessReadsSmallK(uchar *_part, uint64 _part_size, CSmallKBuf<COUNTER_TYPE>& small_k_buf)
{
	return CSplitter_Impl<QUAKE_MODE>::ProcessReadsSmallK(*this, _part, _part_size, small_k_buf);
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part
template <bool QUAKE_MODE> bool CSplitter<QUAKE_MODE>::ProcessReads(uchar *_part, uint64 _part_size)
{
	return CSplitter_Impl<QUAKE_MODE>::ProcessReads(*this, _part, _part_size);
}

//----------------------------------------------------------------------------------
// Return the number of reads processed by splitter
template <bool QUAKE_MODE> void CSplitter<QUAKE_MODE>::GetTotal(uint64 &_n_reads)
{
	_n_reads = n_reads;
}

//----------------------------------------------------------------------------------
// Return the number of kmers processed by splitter (!!! only for small k optimization)
template <bool QUAKE_MODE> uint64 CSplitter<QUAKE_MODE>::GetTotalKmers()
{
	return total_kmers;
}

//************************************************************************************************************
// Implementation of specific splitter methods for various types and sizes of kmers
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part in small k optimization mode
template<typename COUNTER_TYPE>
bool CSplitter_Impl<false>::ProcessReadsSmallK(CSplitter<false> &ptr, uchar *_part, uint64 _part_size, CSmallKBuf<COUNTER_TYPE>& small_k_buf)
{
	ptr.part = _part;
	ptr.part_size = _part_size;
	ptr.part_pos = 0;

	char *seq;
	uint32 seq_size;
	int omit_next_n_kmers;
	CKmer<1> kmer_str, kmer_rev, kmer_can;
	uint32 i;
	CKmer<1> kmer_mask;
	ptr.pmm_reads->reserve(seq);
	kmer_mask.set_n_1(2 * ptr.kmer_len);

	uint32 kmer_len_shift = (ptr.kmer_len - 1) * 2;

	if (ptr.both_strands)
	while (ptr.GetSeq(seq, seq_size))
	{
		if (ptr.file_type != multiline_fasta)
			ptr.n_reads++;

		// Init k-mer
		kmer_str.clear();
		kmer_rev.clear();

		// Process first k-1 symbols of a read
		uint32 str_pos = kmer_len_shift - 2;
		uint32 rev_pos = 2;

		omit_next_n_kmers = 0;

		for (i = 0; i < ptr.kmer_len - 1; ++i, str_pos -= 2, rev_pos += 2)
		{
			if (seq[i] < 0)
			{
				seq[i] = 0;
				omit_next_n_kmers = i + 1;
			}
			kmer_str.set_2bits(seq[i], str_pos);
			kmer_rev.set_2bits(3 - seq[i], rev_pos);
		}

		// Process next part of a read
		for (; i < seq_size; ++i)
		{
			if (seq[i] < 0)		// N in a read
			{
				seq[i] = 0;
				omit_next_n_kmers = ptr.kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
			}
			kmer_str.SHL_insert_2bits(seq[i]);
			kmer_str.mask(kmer_mask);
			kmer_rev.SHR_insert_2bits(3 - seq[i], kmer_len_shift);

			// If necessary ommit next symbols
			if (omit_next_n_kmers > 0)
			{
				omit_next_n_kmers--;
				continue;
			}

			// Find canonical kmer representation
			kmer_can = (kmer_str < kmer_rev) ? kmer_str : kmer_rev;

			++small_k_buf.buf[kmer_can.data];
			++ptr.total_kmers;
		}
	}
	else
	while (ptr.GetSeq(seq, seq_size))
	{
		if (ptr.file_type != multiline_fasta)
			ptr.n_reads++;

		// Init k-mer
		kmer_str.clear();

		// Process first k-1 symbols of a read
		uint32 str_pos = kmer_len_shift - 2;

		omit_next_n_kmers = 0;

		for (i = 0; i < ptr.kmer_len - 1; ++i, str_pos -= 2)
		{
			if (seq[i] < 0)
			{
				seq[i] = 0;
				omit_next_n_kmers = i + 1;
			}
			kmer_str.set_2bits(seq[i], str_pos);
		}

		// Process next part of a read
		for (; i < seq_size; ++i)
		{
			if (seq[i] < 0)		// N in a read
			{
				seq[i] = 0;
				omit_next_n_kmers = ptr.kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
			}
			kmer_str.SHL_insert_2bits(seq[i]);
			kmer_str.mask(kmer_mask);

			// If necessary ommit next symbols
			if (omit_next_n_kmers > 0)
			{
				omit_next_n_kmers--;
				continue;
			}

			++small_k_buf.buf[kmer_str.data];
			++ptr.total_kmers;
		}
	}
	//putchar('*');
	//fflush(stdout);

	ptr.pmm_reads->free(seq);
	return true;
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part
bool CSplitter_Impl<false>::ProcessReads(CSplitter<false> &ptr, uchar *_part, uint64 _part_size)
{
	ptr.part      = _part;
	ptr.part_size = _part_size;
	ptr.part_pos  = 0;

	char *seq;
	uint32 seq_size;
	ptr.pmm_reads->reserve(seq);

	uint32 signature_start_pos;
	CMmer current_signature(ptr.signature_len), end_mmer(ptr.signature_len);
	uint32 bin_no;
	
	uint32 i;
	uint32 len;//length of extended kmer
	
	while (ptr.GetSeq(seq, seq_size))
	{
		if (ptr.file_type != multiline_fasta)
			ptr.n_reads++;
		i = 0;
		len = 0;
		while (i + ptr.kmer_len - 1 < seq_size)
		{
			bool contains_N = false;
			//building first signature after 'N' or at the read begining
			for (uint32 j = 0; j < ptr.signature_len; ++j, ++i)
			if (seq[i] < 0)//'N'
			{
				contains_N = true;
				break;
			}
			//signature must be shorter than k-mer so if signature contains 'N', k-mer will contains it also
			if (contains_N)
			{
				++i;
				continue;
			}
			len = ptr.signature_len;
			signature_start_pos = i - ptr.signature_len;
			current_signature.insert(seq + signature_start_pos);
			end_mmer.set(current_signature);
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)//'N'
				{
					if (len >= ptr.kmer_len)
					{
						bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
						ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, len);
					}
					len = 0;
					++i;
					break;
				}
				end_mmer.insert(seq[i]);
				if (end_mmer < current_signature)//signature at the end of current k-mer is lower than current
				{
					if (len >= ptr.kmer_len)
					{
						bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
						ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, len);
						len = ptr.kmer_len - 1;
					}
					current_signature.set(end_mmer);
					signature_start_pos = i - ptr.signature_len + 1;
				}
				else if (end_mmer == current_signature)
				{
					current_signature.set(end_mmer);
					signature_start_pos = i - ptr.signature_len + 1;
				}
				else if (signature_start_pos + ptr.kmer_len - 1 < i)//need to find new signature
				{
					bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
					ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, len);
					len = ptr.kmer_len - 1;
					//looking for new signature
					++signature_start_pos;
					//building first signature in current k-mer
					end_mmer.insert(seq + signature_start_pos);
					current_signature.set(end_mmer);
					for (uint32 j = signature_start_pos + ptr.signature_len; j <= i; ++j)
					{
						end_mmer.insert(seq[j]);
						if (end_mmer <= current_signature)
						{
							current_signature.set(end_mmer);
							signature_start_pos = j - ptr.signature_len + 1;
						}
					}
				}
				++len;
				if (len == ptr.kmer_len + 255) //one byte is used to store counter of additional symbols in extended k-mer
				{
					bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
					ptr.bins[bin_no]->PutExtendedKmer(seq + i + 1 - len, len);
					i -= ptr.kmer_len - 2;
					len = 0;
					break;
				}

			}
		}
		if (len >= ptr.kmer_len)//last one in read
		{
			bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
			ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, len);
		}
	}
		
	//putchar('*');
	//fflush(stdout);

	ptr.pmm_reads->free(seq);


	return true;
}

bool CSplitter_Impl<true>::ProcessReadsSmallK(CSplitter<true> &ptr, uchar *_part, uint64 _part_size, CSmallKBuf<float>& small_k_buf)
{	
	ptr.part = _part;
	ptr.part_size = _part_size;
	ptr.part_pos = 0;

	char *seq;
	char *quals;
	double *raw_inv_probs;

	ptr.pmm_reads->reserve(seq);
	ptr.pmm_reads->reserve(quals);
	ptr.pmm_reads->reserve(raw_inv_probs);

	double *inv_probs = raw_inv_probs + 1;
	inv_probs[-1] = 1.0;			// !!! Correct

	uint32 seq_size;
	int omit_next_n_kmers;
	CKmer<1> kmer_str, kmer_rev, kmer_can;
	double kmer_prob;

	uint32 i;
	CKmer<1> kmer_mask;	

	kmer_mask.set_n_1(2 * ptr.kmer_len);

	uint32 kmer_len_shift = (ptr.kmer_len - 1) * 2;

	if (ptr.both_strands)
	while (ptr.GetSeq(seq, quals, seq_size))
	{
		ptr.n_reads++;

		// Init k-mer
		kmer_str.clear();
		kmer_rev.clear();

		// Process first k-1 symbols of a read
		uint32 str_pos = kmer_len_shift - 2;
		uint32 rev_pos = 2;

		omit_next_n_kmers = 0;
		kmer_prob = 1.0;

		for (i = 0; i < ptr.kmer_len - 1; ++i, str_pos -= 2, rev_pos += 2)
		{
			if (seq[i] < 0)
			{
				seq[i] = 0;
				omit_next_n_kmers = i + 1;

			}
			inv_probs[i] = CProbQual::inv_prob_qual[quals[i] - ptr.lowest_quality];

			kmer_str.set_2bits(seq[i], str_pos);
			kmer_rev.set_2bits(3 - seq[i], rev_pos);
			kmer_prob *= CProbQual::prob_qual[quals[i] - ptr.lowest_quality];
		}

		// Process next part of a read
		for (; i < seq_size; ++i)
		{
			if (seq[i] < 0)		// N in a read
			{
				seq[i] = 0;
				omit_next_n_kmers = ptr.kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
			}
			inv_probs[i] = CProbQual::inv_prob_qual[quals[i] - ptr.lowest_quality];

			kmer_str.SHL_insert_2bits(seq[i]);
			kmer_str.mask(kmer_mask);
			kmer_rev.SHR_insert_2bits(3 - seq[i], kmer_len_shift);
			kmer_prob *= CProbQual::prob_qual[quals[i] - ptr.lowest_quality] * inv_probs[(int)i - (int)ptr.kmer_len];

			// If necessary ommit next symbols
			if (omit_next_n_kmers > 0)
			{
				omit_next_n_kmers--;
				continue;
			}

			if (kmer_prob < CProbQual::MIN_PROB_QUAL_VALUE)
				continue;

			// Find canonical kmer representation
			kmer_can = (kmer_str < kmer_rev) ? kmer_str : kmer_rev;
			small_k_buf.buf[kmer_can.data] += static_cast<float>(kmer_prob);
			++ptr.total_kmers;
		}
	}
	else
	while (ptr.GetSeq(seq, quals, seq_size))
	{
		ptr.n_reads++;

		// Init k-mer
		kmer_str.clear();

		// Process first k-1 symbols of a read
		uint32 str_pos = kmer_len_shift - 2;

		omit_next_n_kmers = 0;
		kmer_prob = 1.0;

		for (i = 0; i < ptr.kmer_len - 1; ++i, str_pos -= 2)
		{
			if (seq[i] < 0)
			{
				seq[i] = 0;
				omit_next_n_kmers = i + 1;

			}
			inv_probs[i] = CProbQual::inv_prob_qual[quals[i] - ptr.lowest_quality];

			kmer_str.set_2bits(seq[i], str_pos);
			kmer_prob *= CProbQual::prob_qual[quals[i] - ptr.lowest_quality];
		}

		// Process next part of a read
		for (; i < seq_size; ++i)
		{
			if (seq[i] < 0)		// N in a read
			{
				seq[i] = 0;
				omit_next_n_kmers = ptr.kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
			}
			inv_probs[i] = CProbQual::inv_prob_qual[quals[i] - ptr.lowest_quality];

			kmer_str.SHL_insert_2bits(seq[i]);
			kmer_str.mask(kmer_mask);
			kmer_prob *= CProbQual::prob_qual[quals[i] - ptr.lowest_quality] * inv_probs[(int)i - (int)ptr.kmer_len];

			// If necessary ommit next symbols
			if (omit_next_n_kmers > 0)
			{
				omit_next_n_kmers--;
				continue;
			}

			if (kmer_prob < CProbQual::MIN_PROB_QUAL_VALUE)
				continue;

			small_k_buf.buf[kmer_str.data] += static_cast<float>(kmer_prob);
			++ptr.total_kmers;
		}
	}

	//putchar('*');
	//fflush(stdout);

	ptr.pmm_reads->free(seq);
	ptr.pmm_reads->free(quals);
	ptr.pmm_reads->free(raw_inv_probs);

	return true;
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part
bool CSplitter_Impl<true>::ProcessReads(CSplitter<true> &ptr, uchar *_part, uint64 _part_size)
{
	ptr.part      = _part;
	ptr.part_size = _part_size;
	ptr.part_pos  = 0;

	char *seq; 
	char *quals;

	ptr.pmm_reads->reserve(seq);
	ptr.pmm_reads->reserve(quals);


	uint32 seq_size;

	uint32 signature_start_pos;
	CMmer current_signature(ptr.signature_len), end_mmer(ptr.signature_len);
	uint32 bin_no;

	uint32 i;
	uint32 len;//length of extended kmer



	while (ptr.GetSeq(seq, quals, seq_size))
	{
		if (ptr.file_type != multiline_fasta)
			ptr.n_reads++;
		i = 0;
		len = 0;
		while (i + ptr.kmer_len - 1 < seq_size)
		{
			bool contains_N = false;
			//building first signature after 'N' or at the read begining
			for (uint32 j = 0; j < ptr.signature_len; ++j, ++i)
			if (seq[i] < 0)//'N'
			{
				contains_N = true;
				break;
			}
			//signature must be shorter than k-mer so if signature contains 'N', k-mer will contains it also
			if (contains_N)
			{
				++i;
				continue;
			}
			len = ptr.signature_len;
			signature_start_pos = i - ptr.signature_len;
			current_signature.insert(seq + signature_start_pos);
			end_mmer.set(current_signature);
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)//'N'
				{
					if (len >= ptr.kmer_len)
					{
						bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
						ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, quals + i - len, len);
					}
					len = 0;
					++i;
					break;
				}
				end_mmer.insert(seq[i]);
				if (end_mmer < current_signature)//signature at the end of current k-mer is lower than current
				{
					if (len >= ptr.kmer_len)
					{
						bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
						ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, quals + i - len, len);
						len = ptr.kmer_len - 1;
					}
					current_signature.set(end_mmer);
					signature_start_pos = i - ptr.signature_len + 1;
				}
				else if (end_mmer == current_signature)
				{
					current_signature.set(end_mmer);
					signature_start_pos = i - ptr.signature_len + 1;
				}
				else if (signature_start_pos + ptr.kmer_len - 1 < i)//need to find new signature
				{
					bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
					ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, quals + i - len, len);
					len = ptr.kmer_len - 1;
					//looking for new signature
					++signature_start_pos;
					//building first signature in current k-mer
					end_mmer.insert(seq + signature_start_pos);
					current_signature.set(end_mmer);
					for (uint32 j = signature_start_pos + ptr.signature_len; j <= i; ++j)
					{
						end_mmer.insert(seq[j]);
						if (end_mmer <= current_signature)
						{
							current_signature.set(end_mmer);
							signature_start_pos = j - ptr.signature_len + 1;
						}
					}
				}
				++len;
				if (len == ptr.kmer_len + 255) //one byte is used to store counter of additional symbols in extended k-mer
				{
					bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
					ptr.bins[bin_no]->PutExtendedKmer(seq + i + 1 - len, quals + i + 1 - len, len);
					i -= ptr.kmer_len - 2;
					len = 0;
					break;
				}
			}
		}
		if (len >= ptr.kmer_len)//last one in read
		{
			bin_no = ptr.s_mapper->get_bin_id(current_signature.get());
			ptr.bins[bin_no]->PutExtendedKmer(seq + i - len, quals + i - len, len);
		}
	}


	//putchar('*');
	//fflush(stdout);

	
	ptr.pmm_reads->free(seq);
	ptr.pmm_reads->free(quals);	

	return true;
}

//************************************************************************************************************
// CWSplitter class - wrapper for multithreading purposes
//************************************************************************************************************

//----------------------------------------------------------------------------------
template <bool QUAKE_MODE> class CWSplitter {
	CPartQueue *pq;
	CBinPartQueue *bpq;
	CMemoryPool *pmm_fastq;

	CSplitter<QUAKE_MODE> *spl;
	uint64 n_reads;

public:
	CWSplitter(CKMCParams &Params, CKMCQueues &Queues);
	~CWSplitter();

	void operator()();
	void GetTotal(uint64 &_n_reads);
};

//----------------------------------------------------------------------------------
// Constructor
template <bool QUAKE_MODE> CWSplitter<QUAKE_MODE>::CWSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	pq		  = Queues.part_queue;
	bpq		  = Queues.bpq;
	pmm_fastq = Queues.pmm_fastq;
	spl = new CSplitter<QUAKE_MODE>(Params, Queues);
	spl->InitBins(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
template <bool QUAKE_MODE> CWSplitter<QUAKE_MODE>::~CWSplitter()
{
}

//----------------------------------------------------------------------------------
// Execution
template <bool QUAKE_MODE> void CWSplitter<QUAKE_MODE>::operator()()
{
	// Splitting parts
	while(!pq->completed())
	{
		uchar *part;
		uint64 size;
		if(pq->pop(part, size))
		{
			spl->ProcessReads(part, size);
			pmm_fastq->free(part);
		}
	}
	spl->Complete();
	bpq->mark_completed();

	spl->GetTotal(n_reads);

	delete spl;
	spl = NULL;
}

//----------------------------------------------------------------------------------
// Return statistics
template <bool QUAKE_MODE> void CWSplitter<QUAKE_MODE>::GetTotal(uint64 &_n_reads)
{
	if(spl)
		spl->GetTotal(n_reads);

	_n_reads = n_reads;
}



//************************************************************************************************************
// CWStatsSplitter class - wrapper for multithreading purposes
//************************************************************************************************************

//----------------------------------------------------------------------------------
template <bool QUAKE_MODE> class CWStatsSplitter {
	CStatsPartQueue *spq;
	CMemoryPool *pmm_fastq, *pmm_stats;
	uint32 *stats;
	CSplitter<QUAKE_MODE> *spl;
	uint32 signature_len;

public:
	CWStatsSplitter(CKMCParams &Params, CKMCQueues &Queues);
	~CWStatsSplitter();

	void operator()();
	void GetStats(uint32* _stats);
};

//----------------------------------------------------------------------------------
// Constructor
template <bool QUAKE_MODE> CWStatsSplitter<QUAKE_MODE>::CWStatsSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	spq = Queues.stats_part_queue;
	pmm_fastq = Queues.pmm_fastq;
	pmm_stats = Queues.pmm_stats;
	spl = new CSplitter<QUAKE_MODE>(Params, Queues);
	
	signature_len = Params.signature_len;
	pmm_stats->reserve(stats);
	fill_n(stats, (1 << signature_len * 2) + 1, 0);
}

//----------------------------------------------------------------------------------
// Destructor
template <bool QUAKE_MODE> CWStatsSplitter<QUAKE_MODE>::~CWStatsSplitter()
{
	pmm_stats->free(stats);
}

//----------------------------------------------------------------------------------
// Execution
template <bool QUAKE_MODE> void CWStatsSplitter<QUAKE_MODE>::operator()()
{
	// Splitting parts
	while (!spq->completed())
	{
		uchar *part;
		uint64 size;
		if (spq->pop(part, size))
		{
			spl->CalcStats(part, size, stats);
			pmm_fastq->free(part);
		}
	}

	delete spl;
	spl = NULL;
}

//----------------------------------------------------------------------------------
template <bool QUAKE_MODE> void CWStatsSplitter<QUAKE_MODE>::GetStats(uint32* _stats)
{
	uint32 size = (1 << signature_len * 2) + 1;
	for (uint32 i = 0; i < size; ++i)
		_stats[i] += stats[i];
}


//************************************************************************************************************
// CWSmallKSplitter class - wrapper for multithreading purposes
//************************************************************************************************************
//----------------------------------------------------------------------------------
template <bool QUAKE_MODE, typename COUNTER_TYPE> class CWSmallKSplitter {
	CPartQueue *pq;	
	CMemoryPool *pmm_fastq, *pmm_small_k;	
	CSmallKBuf<COUNTER_TYPE> small_k_buf;

	CSplitter<QUAKE_MODE> *spl;
	uint64 n_reads;
	uint64 total_kmers;
	uint32 kmer_len;

public:
	CWSmallKSplitter(CKMCParams &Params, CKMCQueues &Queues);
	~CWSmallKSplitter();

	void operator()();
	void GetTotal(uint64 &_n_reads);

	CSmallKBuf<COUNTER_TYPE> GetResult()
	{
		return small_k_buf;
	}

	uint64 GetTotalKmers()
	{
		if (spl)
			return spl->GetTotalKmers();
		return total_kmers;
	}

	void Release()
	{
		pmm_small_k->free(small_k_buf.buf);
	}
};

//----------------------------------------------------------------------------------
// Constructor
template <bool QUAKE_MODE, typename COUNTER_TYPE> CWSmallKSplitter<QUAKE_MODE, COUNTER_TYPE>::CWSmallKSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	pq = Queues.part_queue;	
	pmm_fastq = Queues.pmm_fastq;
	pmm_small_k = Queues.pmm_small_k_buf;
	kmer_len = Params.kmer_len;
	spl = new CSplitter<QUAKE_MODE>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
template <bool QUAKE_MODE, typename COUNTER_TYPE> CWSmallKSplitter<QUAKE_MODE, COUNTER_TYPE>::~CWSmallKSplitter()
{
}

//----------------------------------------------------------------------------------
// Execution
template <bool QUAKE_MODE, typename COUNTER_TYPE> void CWSmallKSplitter<QUAKE_MODE, COUNTER_TYPE>::operator()()
{
	pmm_small_k->reserve(small_k_buf.buf);
	memset(small_k_buf.buf, 0, (1ull << 2 * kmer_len) * sizeof(*small_k_buf.buf)); 

	// Splitting parts
	while (!pq->completed())
	{
		uchar *part;
		uint64 size;
		if (pq->pop(part, size))
		{
			spl->ProcessReadsSmallK(part, size, small_k_buf);
			pmm_fastq->free(part);
		}
	}
	spl->Complete();

	spl->GetTotal(n_reads);
	total_kmers = spl->GetTotalKmers();
	delete spl;
	spl = NULL;
}

//----------------------------------------------------------------------------------
// Return statistics
template <bool QUAKE_MODE, typename COUNTER_TYPE> void CWSmallKSplitter<QUAKE_MODE, COUNTER_TYPE>::GetTotal(uint64 &_n_reads)
{
	if (spl)
		spl->GetTotal(n_reads);

	_n_reads = n_reads;
}


#endif

// ***** EOF
