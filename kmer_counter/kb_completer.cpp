#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/
#include <algorithm>
#include <numeric>
#include <iostream>
#include "kb_completer.h"

using namespace std;

extern uint64 total_reads;



//************************************************************************************************************
// CKmerBinCompleter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Assign queues and monitors
CKmerBinCompleter::CKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues) 
{
	mm		       = Queues.mm;
	file_name      = Params.output_file_name;
	kq             = Queues.kq;
	bd		       = Queues.bd;
	s_mapper	   = Queues.s_mapper;
	memory_bins    = Queues.memory_bins;

	bbkpq		   = Queues.bbkpq;
	use_strict_mem = Params.use_strict_mem;
	kmer_file_name = file_name + ".kmc_suf";
	lut_file_name  = file_name + ".kmc_pre";

	kmer_len       = Params.kmer_len;
	signature_len  = Params.signature_len;

	cutoff_min     = Params.cutoff_min;
	cutoff_max     = (uint32)Params.cutoff_max;
	counter_max    = (uint32)Params.counter_max;
	lut_prefix_len = Params.lut_prefix_len;
	both_strands   = Params.both_strands;
	without_output = Params.without_output;

	kmer_t_size    = Params.KMER_T_size;
		
}

//----------------------------------------------------------------------------------
CKmerBinCompleter::~CKmerBinCompleter()
{
}


//----------------------------------------------------------------------------------
// Store sorted and compacted bins to the output file (stage first)
void CKmerBinCompleter::ProcessBinsFirstStage()
{
	int32 bin_id = 0;
	uchar *data = nullptr;
	//uint64 data_size = 0;
	list<pair<uint64, uint64>> data_packs;
	uchar *lut = nullptr;
	uint64 lut_size = 0;
	counter_size = 0;
	sig_map_size = (1 << (signature_len * 2)) + 1;
	sig_map = new uint32[sig_map_size];
	fill_n(sig_map, sig_map_size, 0);
	lut_pos = 0;
	
	counter_size = min(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));	
	
	if (!without_output)
	{
		// Open output file
		out_kmer = fopen(kmer_file_name.c_str(), "wb");
		if (!out_kmer)
		{
			cerr << "Error: Cannot create " << kmer_file_name << "\n";
			exit(1);
			return;
		}

		out_lut = fopen(lut_file_name.c_str(), "wb");
		if (!out_lut)
		{
			cerr << "Error: Cannot create " << lut_file_name << "\n";
			fclose(out_kmer);
			exit(1);
			return;
		}
	}
	
	n_recs = 0;

	_n_unique = _n_cutoff_min = _n_cutoff_max = _n_total = 0;
	n_unique  = n_cutoff_min  = n_cutoff_max  = n_total  = 0;

	char s_kmc_pre[] = "KMCP";
	char s_kmc_suf[] = "KMCS";

	if (!without_output)
	{
		// Markers at the beginning
		fwrite(s_kmc_pre, 1, 4, out_lut);
		fwrite(s_kmc_suf, 1, 4, out_kmer);
	}

	// Process priority queue of ready-to-output bins
	while (!kq->empty())
	{
		// Get the next bin
		if (!kq->pop(bin_id, data, data_packs, lut, lut_size, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total))
			continue;

		// Decrease memory size allocated by stored bin
		string name;
		uint64 n_rec;
		uint64 n_plus_x_recs;
		uint64 n_super_kmers;
		uint64 raw_size;
		CMemDiskFile *file;

		bd->read(bin_id, file, name, raw_size, n_rec, n_plus_x_recs, n_super_kmers);

		uint64 lut_recs = lut_size / sizeof(uint64);


		if (!without_output)
		{
			for (auto& e : data_packs)
			{
				// Write bin data to the output file
#ifdef WIN32 //fwrite bug https://connect.microsoft.com/VisualStudio/feedback/details/755018/fwrite-hangs-with-large-size-count
				uint64 write_offset = e.first;
				uint64 left_to_write = e.second - e.first;
				while (left_to_write)
				{
					uint64 current_to_write = MIN(left_to_write, (4ull << 30) - 1);
					fwrite(data + write_offset, 1, current_to_write, out_kmer);
					write_offset += current_to_write;
					left_to_write -= current_to_write;
				}
#else
				fwrite(data + e.first, 1, e.second - e.first, out_kmer);
#endif
			}
		}

		memory_bins->free(bin_id, CMemoryBins::mba_suffix);

		if (!without_output)
		{
			uint64 *ulut = (uint64*)lut;
			for (uint64 i = 0; i < lut_recs; ++i)
			{
				uint64 x = ulut[i];
				ulut[i] = n_recs;
				n_recs += x;
			}
			fwrite(lut, lut_recs, sizeof(uint64), out_lut);
		}
		//fwrite(&n_rec, 1, sizeof(uint64), out_lut);
		memory_bins->free(bin_id, CMemoryBins::mba_lut);

		n_unique	 += _n_unique;
		n_cutoff_min += _n_cutoff_min;
		n_cutoff_max += _n_cutoff_max;
		n_total      += _n_total;
		for (uint32 i = 0; i < sig_map_size; ++i)
		{
			if (s_mapper->get_bin_id(i) == bin_id)
			{
				sig_map[i] = lut_pos;
			}
		}
		++lut_pos;
	}		
}

//----------------------------------------------------------------------------------
// Store sorted and compacted bins to the output file (stage second)
void CKmerBinCompleter::ProcessBinsSecondStage()
{
	char s_kmc_pre[] = "KMCP";
	char s_kmc_suf[] = "KMCS";
	if (use_strict_mem)
	{
		int32 bin_id;
		uchar *data = nullptr;
		uint64 data_size = 0;
		uchar *lut = nullptr;
		uint64 lut_size = 0;		
		bool last_in_bin = false;
		while (bbkpq->pop(bin_id, data, data_size, lut, lut_size, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, last_in_bin))
		{
			if (data_size)
			{
				if(!without_output)
					fwrite(data, 1, data_size, out_kmer);				
				sm_pmm_merger_suff->free(data);
			}
			if (lut_size)
			{
				uint64 lut_recs = lut_size / sizeof(uint64);
				uint64* ulut = (uint64*)lut;
				for (uint64 i = 0; i < lut_recs; ++i)
				{
					uint64 x = ulut[i];
					ulut[i] = n_recs;
					n_recs += x;
				}
				if(!without_output)
					fwrite(lut, lut_recs, sizeof(uint64), out_lut);				
				sm_pmm_merger_lut->free(lut);
			}
			if(last_in_bin)
			{
				n_unique += _n_unique;
				n_cutoff_min += _n_cutoff_min;
				n_cutoff_max += _n_cutoff_max;
				n_total += _n_total;
				for (uint32 i = 0; i < sig_map_size; ++i)
				{
					if (s_mapper->get_bin_id(i) == bin_id)
					{
						sig_map[i] = lut_pos;
					}
				}
				++lut_pos;
			}
		}
	}

	if (!without_output)
	{
		// Marker at the end
		fwrite(s_kmc_suf, 1, 4, out_kmer);
		fclose(out_kmer);

		fwrite(&n_recs, 1, sizeof(uint64), out_lut);

		//store signature mapping 
		fwrite(sig_map, sizeof(uint32), sig_map_size, out_lut);

		// Store header
		uint32 offset = 0;

		store_uint(out_lut, kmer_len, 4);				offset += 4;
		store_uint(out_lut, (uint32)0, 4);				offset += 4;	// mode: 0 (counting), 1 (Quake-compatibile counting) which is now not supported
		store_uint(out_lut, counter_size, 4);			offset += 4;
		store_uint(out_lut, lut_prefix_len, 4);			offset += 4;
		store_uint(out_lut, signature_len, 4);			offset += 4;
		store_uint(out_lut, cutoff_min, 4);				offset += 4;
		store_uint(out_lut, cutoff_max, 4);				offset += 4;
		store_uint(out_lut, n_unique - n_cutoff_min - n_cutoff_max, 8);		offset += 8;

		store_uint(out_lut, both_strands ? 0 : 1, 1);			offset++;

		// Space for future use
		for (int32 i = 0; i < 27; ++i)
		{
			store_uint(out_lut, 0, 1);
			offset++;
		}

		store_uint(out_lut, 0x200, 4);
		offset += 4;

		store_uint(out_lut, offset, 4);

		// Marker at the end
		fwrite(s_kmc_pre, 1, 4, out_lut);
		fclose(out_lut);
	}
	cerr << "\n";

	delete[] sig_map;
}

//----------------------------------------------------------------------------------
// Return statistics
void CKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total)
{
	_n_unique	  = n_unique;
	_n_cutoff_min = n_cutoff_min;
	_n_cutoff_max = n_cutoff_max;
	_n_total      = n_total;
}

//----------------------------------------------------------------------------------
// Store single unsigned integer in LSB fashion
bool CKmerBinCompleter::store_uint(FILE *out, uint64 x, uint32 size)
{
	for(uint32 i = 0; i < size; ++i)
		putc((x >> (i * 8)) & 0xFF, out);

	return true;
}

//----------------------------------------------------------------------------------
//Init memory pools for 2nd stage
void CKmerBinCompleter::InitStage2(CKMCParams& /*Params*/, CKMCQueues& Queues)
{
	sm_pmm_merger_lut = Queues.sm_pmm_merger_lut;
	sm_pmm_merger_suff = Queues.sm_pmm_merger_suff;
}


//************************************************************************************************************
// CWKmerBinCompleter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Constructor
CWKmerBinCompleter::CWKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues)
{
	kbc = new CKmerBinCompleter(Params, Queues);
}

void CWKmerBinCompleter::InitStage2(CKMCParams& Params, CKMCQueues& Queues)
{
	kbc->InitStage2(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
CWKmerBinCompleter::~CWKmerBinCompleter()
{
	delete kbc;
}

//----------------------------------------------------------------------------------
// Execution
void CWKmerBinCompleter::operator()(bool first_stage)
{
	if(first_stage)
		kbc->ProcessBinsFirstStage();
	else
		kbc->ProcessBinsSecondStage();
}

//----------------------------------------------------------------------------------
// Return statistics
void CWKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total)
{
	if(kbc)
		kbc->GetTotal(_n_unique, _n_cutoff_min, _n_cutoff_max, _n_total);
}

// ***** EOF
