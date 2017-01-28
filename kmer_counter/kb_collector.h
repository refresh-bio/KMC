/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#ifndef _KB_COLLECTOR_H
#define _KB_COLLECTOR_H

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include "queues.h"
#include "radix.h"
#include "rev_byte.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <vector>
#include <stdio.h>

using namespace std;



//----------------------------------------------------------------------------------
// Class collecting kmers belonging to a single bin
class CKmerBinCollector
{
	list<pair<uint64, uint64>> expander_parts; //range, n_plus_x_recs_in_range
	uint64 prev_n_plus_x_recs = 0;
	uint64 prev_pos = 0;

	enum comparision_state  { kmer_smaller, rev_smaller, equals };
	uint32 bin_no;
	CBinPartQueue *bin_part_queue;
	CBinDesc *bd;
	uint32 kmer_len;
	uchar* buffer;
	uint32 buffer_size;
	uint32 buffer_pos;

	uint32 super_kmer_no = 0;
	const uint32 max_super_kmers_expander_pack = 1ul << 12; 
	
	CMemoryPool *pmm_bins;
	uint32 n_recs;
	uint32 n_plus_x_recs;
	uint32 n_super_kmers;
	int lowest_quality;
	uint32 max_x;

	uint32 kmer_bytes;
	bool both_strands;

	template<unsigned DIVIDE_FACTOR> void update_n_plus_x_recs(char* seq, uint32 n);

public:
	CKmerBinCollector(CKMCQueues& Queues, CKMCParams& Params, uint32 _buffer_size, uint32 _bin_no);
	void PutExtendedKmer(char* seq, uint32 n);
	void PutExtendedKmer(char* seq, char* quals, uint32 n);//for quake mode
	inline void Flush();
};

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

CKmerBinCollector::CKmerBinCollector(CKMCQueues& Queues, CKMCParams& Params, uint32 _buffer_size, uint32 _bin_no)
{
	bin_part_queue	= Queues.bpq;
	kmer_len		= Params.kmer_len;
	bd				= Queues.bd;
	buffer_size		= _buffer_size;	
	pmm_bins		= Queues.pmm_bins;
	lowest_quality  = Params.lowest_quality;
	max_x				= Params.max_x;
	
	bin_no			= _bin_no;
	n_recs			= 0;
	n_super_kmers		= 0;
	n_plus_x_recs   = 0;
	buffer_pos		= 0;	
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
	if(buffer_pos + bytes > buffer_size)
	{		
		//send current buff
		Flush();

		pmm_bins->reserve(buffer);
		buffer_pos		= 0;
		n_recs			= 0;
		n_super_kmers		= 0;
		n_plus_x_recs	= 0;
	}
	
	
	buffer[buffer_pos++] = n - kmer_len;		
	for(uint32 i = 0, j = 0 ; i < n / 4 ; ++i,j+=4)
		buffer[buffer_pos++] = (seq[j] << 6) + (seq[j + 1] << 4) + (seq[j + 2] << 2) + seq[j + 3];
	switch (n%4)
	{
	case 1:
		buffer[buffer_pos++] = (seq[n-1] << 6);
		break;
	case 2:
		buffer[buffer_pos++] = (seq[n-2] << 6) + (seq[n-1] << 4);
		break;
	case 3:
		buffer[buffer_pos++] = (seq[n-3] << 6) + (seq[n-2] << 4) + (seq[n-1] << 2);
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
template<unsigned DIVIDE_FACTOR> void CKmerBinCollector::update_n_plus_x_recs(char* seq, uint32 n)
{
	uchar kmer, rev;
	uint32 kmer_pos = 4;
	uint32 rev_pos = kmer_len;
	uint32 x;

	kmer = (seq[0] << 6) + (seq[1] << 4) + (seq[2] << 2) + seq[3];
	rev = ((3 - seq[kmer_len - 1]) << 6) + ((3 - seq[kmer_len - 2]) << 4) + ((3 - seq[kmer_len - 3]) << 2) + (3 - seq[kmer_len - 4]);

	x = 0;
	comparision_state current_state, new_state;
	if (kmer < rev)
		current_state = kmer_smaller;
	else if (rev < kmer)
		current_state = rev_smaller;
	else
		current_state = equals;


	for (uint32 i = 0; i < n - kmer_len; ++i)
	{
		rev >>= 2;
		rev += (3 - seq[rev_pos++]) << 6;
		kmer <<= 2;
		kmer += seq[kmer_pos++];

		if (kmer < rev)
			new_state = kmer_smaller;
		else if (rev < kmer)
			new_state = rev_smaller;
		else
			new_state = equals;

		if (new_state == current_state)
		{
			if (current_state == equals)
				++n_plus_x_recs;
			else
				++x;
		}
		else
		{
			current_state = new_state;
			n_plus_x_recs += 1 + x / DIVIDE_FACTOR;

			x = 0;
		}
	}
	n_plus_x_recs += 1 + x / DIVIDE_FACTOR;
}

//---------------------------------------------------------------------------------
void CKmerBinCollector::PutExtendedKmer(char* seq, char* quals, uint32 n)
{
	uint32 bytes = n + 1;
	if (buffer_pos + bytes > buffer_size)
	{
		Flush();

		pmm_bins->reserve(buffer);
		buffer_pos	= 0;
		n_recs		= 0;
		n_super_kmers	= 0;
	}

	n_recs += n - kmer_len + 1;
	++n_super_kmers;
	buffer[buffer_pos++] = n - kmer_len;
	char qual;
	for (uint32 i = 0; i < n; ++i)
	{
		qual = quals[i] - lowest_quality;
		if (qual > 63)
			qual = 63;
		buffer[buffer_pos++] = (seq[i] << 6) + qual;
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
	bd->insert(bin_no, NULL, "", buffer_pos, n_recs, n_plus_x_recs, n_super_kmers);
}

#endif

// ***** EOF
