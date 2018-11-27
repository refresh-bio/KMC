/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
*/

#include "stdafx.h"
#include "kb_collector.h"


CKmerBinCollector::CKmerBinCollector(CKMCQueues& Queues, CKMCParams& Params, uint32 _buffer_size, uint32 _bin_no)
{
	bin_part_queue = Queues.bpq;
	kmer_len = Params.kmer_len;
	bd = Queues.bd;
	buffer_size = _buffer_size;
	pmm_bins = Queues.pmm_bins;
	max_x = Params.max_x;

	bin_no = _bin_no;
	n_recs = 0;
	n_super_kmers = 0;
	n_plus_x_recs = 0;
	buffer_pos = 0;
	pmm_bins->reserve(buffer);

	both_strands = Params.both_strands;
	kmer_bytes = (kmer_len + 3) / 4;
}
//---------------------------------------------------------------------------------
void CKmerBinCollector::PutExtendedKmer(char* seq, uint32 n)
{
	if (super_kmer_no >= max_super_kmers_expander_pack)
	{
		expander_parts.push_back(make_pair(buffer_pos - prev_pos, n_plus_x_recs - prev_n_plus_x_recs));
		prev_pos = buffer_pos;
		prev_n_plus_x_recs = n_plus_x_recs;
		super_kmer_no = 0;
	}
	uint32 bytes = 1 + (n + 3) / 4;
	if (buffer_pos + bytes > buffer_size)
	{
		//send current buff
		Flush();

		pmm_bins->reserve(buffer);
		buffer_pos = 0;
		n_recs = 0;
		n_super_kmers = 0;
		n_plus_x_recs = 0;
	}


	buffer[buffer_pos++] = n - kmer_len;
	for (uint32 i = 0, j = 0; i < n / 4; ++i, j += 4)
		buffer[buffer_pos++] = (seq[j] << 6) + (seq[j + 1] << 4) + (seq[j + 2] << 2) + seq[j + 3];
	switch (n % 4)
	{
	case 1:
		buffer[buffer_pos++] = (seq[n - 1] << 6);
		break;
	case 2:
		buffer[buffer_pos++] = (seq[n - 2] << 6) + (seq[n - 1] << 4);
		break;
	case 3:
		buffer[buffer_pos++] = (seq[n - 3] << 6) + (seq[n - 2] << 4) + (seq[n - 1] << 2);
		break;
	}

	++n_super_kmers;
	n_recs += n - kmer_len + 1;
	if (max_x) ///for max_x = 0 k-mers (not k+x-mers) will be sorted
	{
		if (!both_strands)
			n_plus_x_recs += 1 + (n - kmer_len) / (max_x + 1);
		else
		{
			switch (max_x)
			{
			case 1: update_n_plus_x_recs<2>(seq, n); break;
			case 2: update_n_plus_x_recs<3>(seq, n); break;
			case 3: update_n_plus_x_recs<4>(seq, n); break;
			}

		}
	}
}

//---------------------------------------------------------------------------------
void CKmerBinCollector::Flush()
{
	if (prev_pos < buffer_pos)
	{
		expander_parts.push_back(make_pair(buffer_pos - prev_pos, n_plus_x_recs - prev_n_plus_x_recs));
	}
	prev_pos = 0;
	prev_n_plus_x_recs = 0;
	super_kmer_no = 0;

	bin_part_queue->push(bin_no, buffer, buffer_pos, buffer_size, expander_parts);
	expander_parts.clear();
	bd->insert(bin_no, nullptr, "", buffer_pos, n_recs, n_plus_x_recs, n_super_kmers);
}

// ***** EOF