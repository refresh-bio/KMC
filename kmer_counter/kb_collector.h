/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
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
	uint32 max_x;

	uint32 kmer_bytes;
	bool both_strands;

	template<unsigned DIVIDE_FACTOR> void update_n_plus_x_recs(char* seq, uint32 n);

public:
	CKmerBinCollector(CKMCQueues& Queues, CKMCParams& Params, uint32 _buffer_size, uint32 _bin_no);
	void PutExtendedKmer(char* seq, uint32 n);	
	void Flush();
};

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

#endif

// ***** EOF
