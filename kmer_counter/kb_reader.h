/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.2.0
  Date   : 2015-04-15
*/

#ifndef _KB_READER_H
#define _KB_READER_H

#include "../kmc/definitions.h"
#include "params.h"
#include "kmer.h"
#include "s_mapper.h"
#include "radix.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <vector>
#include <cstdio>


//************************************************************************************************************
// CKmerBinReader - reader of bins from distribution phase
//************************************************************************************************************
template <typename KMER_T, unsigned SIZE> class CKmerBinReader {
	CMemoryMonitor *mm;
	CSignatureMapper* s_mapper;

	CBinDesc *bd;
	CBinQueue *bq;
	CTooLargeBinsQueue *tlbq;

	CMemoryBins *memory_bins;
	CDiskLogger* disk_logger;

	int32 cutoff_min, cutoff_max;
	int32 counter_max;
	int32 kmer_len;
	int32 lut_prefix_len;
	uint32 max_x;

	bool both_strands;
	bool use_quake;

	int64 round_up_to_alignment(int64 x)
	{
		return (x + ALIGNMENT-1) / ALIGNMENT * ALIGNMENT;
	}

public:
	CKmerBinReader(CKMCParams &Params, CKMCQueues &Queues);
	~CKmerBinReader();

	void ProcessBins();
};


//----------------------------------------------------------------------------------
// Assign monitors and queues
template <typename KMER_T, unsigned SIZE> CKmerBinReader<KMER_T, SIZE>::CKmerBinReader(CKMCParams &Params, CKMCQueues &Queues)
{
	mm = Queues.mm;
//	dm = Queues.dm;
	bd = Queues.bd;
	bq = Queues.bq;
	tlbq = Queues.tlbq;
	disk_logger = Queues.disk_logger;

	memory_bins = Queues.memory_bins;

	kmer_len       = Params.kmer_len;	
	cutoff_min     = Params.cutoff_min;
	cutoff_max     = Params.cutoff_max;
	counter_max    = Params.counter_max;
	both_strands   = Params.both_strands;
	use_quake = Params.use_quake;
	max_x = Params.max_x;
	s_mapper	   = Queues.s_mapper;
	lut_prefix_len = Params.lut_prefix_len;
}

//----------------------------------------------------------------------------------
template <typename KMER_T, unsigned SIZE> CKmerBinReader<KMER_T, SIZE>::~CKmerBinReader()
{
	
}

//----------------------------------------------------------------------------------
// Read all bins from temporary HDD
template <typename KMER_T, unsigned SIZE> void CKmerBinReader<KMER_T, SIZE>::ProcessBins()
{
	uchar *data;
	uint64 readed;

	int32 bin_id;
	CMemDiskFile *file;
	string name;
	uint64 size;
	uint64 n_rec;
	uint64 n_plus_x_recs;
	uint32 buffer_size;
	uint32 kmer_len;

	bd->init_random();
	while((bin_id = bd->get_next_random_bin()) >= 0)		// Get id of the next bin to read
	{
		bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, buffer_size, kmer_len);
		fflush(stdout);


		// Reserve memory necessary to process the current bin at all next stages
		uint64 input_kmer_size;
		uint64 kxmer_counter_size;
		uint32 kxmer_symbols;
		if (max_x && !use_quake)
		{
			input_kmer_size = n_plus_x_recs * sizeof(KMER_T);
			kxmer_counter_size = n_plus_x_recs * sizeof(uint32);
			kxmer_symbols = kmer_len + max_x + 1;
		}
		else
		{
			input_kmer_size = n_rec * sizeof(KMER_T); 
			kxmer_counter_size = 0;
			kxmer_symbols = kmer_len;
		}
		uint64 max_out_recs    = (n_rec+1) / max(cutoff_min, 1);	
		
		uint64 counter_size    = min(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));
		if(KMER_T::QUALITY_SIZE > counter_size)
			counter_size = KMER_T::QUALITY_SIZE;

		uint32 kmer_symbols = kmer_len - lut_prefix_len;
		uint64 kmer_bytes = kmer_symbols / 4;
		uint64 out_buffer_size = max_out_recs * (kmer_bytes + counter_size);
			
		uint32 rec_len         = (kxmer_symbols + 3) / 4;

		uint64 lut_recs = 1 << (2 * lut_prefix_len);
		uint64 lut_size = lut_recs * sizeof(uint64);

		if (!memory_bins->init(bin_id, rec_len, round_up_to_alignment(size), round_up_to_alignment(input_kmer_size), round_up_to_alignment(out_buffer_size), round_up_to_alignment(kxmer_counter_size), round_up_to_alignment(lut_size)))
		{
			tlbq->insert(bin_id);
			continue;
		}

#ifdef DEBUG_MODE
		cout << bin_id << ":  " << name << "  " << c_disk << "  " << size << "  " << n_rec << "\n";
#else
		cout << "*";
#endif

		// Process the bin if it is not empty
		if(size > 0)
		{
			if (file == NULL)
			{
				cout << "Error: Cannot open temporary file: " << name << "\n"; fflush(stdout);
				exit(1);
			}
			else
				file->Rewind();

			memory_bins->reserve(bin_id, data, CMemoryBins::mba_input_file);
			//readed = fread(data, 1, size, file);
			readed = file->Read(data, 1, size);
			if(readed != size)
			{
				cout << "Error: Corrupted file: " << name << "   " << "Real size : " << readed << "   " << "Should be : " << size << "\n";
				fflush(stdout);
				exit(1);
			}

			// Push bin data to a queue of bins to process
			bq->push(bin_id, data, size, n_rec);
		}
		else
			// Push empty bin to process (necessary, since all bin ids must be processed)
			bq->push(bin_id, NULL, 0, 0);

		file->Close();
		
		//Remove temporary file
#ifndef DEVELOP_MODE
		file->Remove();		
#endif
		disk_logger->log_remove(size);
	}


	bq->mark_completed();

	fflush(stdout);
}


//************************************************************************************************************
// CWKmerBinReader - wrapper for multithreading purposes
//************************************************************************************************************

//----------------------------------------------------------------------------------
template <typename KMER_T, unsigned SIZE> class CWKmerBinReader {
	CKmerBinReader<KMER_T, SIZE> *kbr;

public:
	CWKmerBinReader(CKMCParams &Params, CKMCQueues &Queues);
	~CWKmerBinReader();

	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template <typename KMER_T, unsigned SIZE> CWKmerBinReader<KMER_T, SIZE>::CWKmerBinReader(CKMCParams &Params, CKMCQueues &Queues)
{
	kbr = new CKmerBinReader<KMER_T, SIZE>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
template <typename KMER_T, unsigned SIZE> CWKmerBinReader<KMER_T, SIZE>::~CWKmerBinReader()
{
	delete kbr;
}

//----------------------------------------------------------------------------------
// Execution
template <typename KMER_T, unsigned SIZE> void CWKmerBinReader<KMER_T, SIZE>::operator()()
{
	kbr->ProcessBins();
}

#endif

// ***** EOF


