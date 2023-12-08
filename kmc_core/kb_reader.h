/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.3
  Date   : 2023-12-08
*/

#ifndef _KB_READER_H
#define _KB_READER_H

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include "s_mapper.h"
#include "radix.h"
#include "percent_progress.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <vector>
#include <stdio.h>


//************************************************************************************************************
// CKmerBinReader - reader of bins from distribution phase
//************************************************************************************************************
template <unsigned SIZE> class CKmerBinReader
{
	CSignatureMapper* s_mapper;

	CBinDesc *bd;
	CBinQueue *bq;
	CSortersManager* sorters_manager;
	CTooLargeBinsQueue *tlbq;

	CMemoryBins *memory_bins;
	CDiskLogger* disk_logger;

	uint32 cutoff_min, cutoff_max;
	uint32 counter_max;
	uint32 kmer_len;
	int32 lut_prefix_len;
	uint32 max_x;

	bool both_strands;	

	KMC::IPercentProgressObserver* percentProgressObserver;
#ifdef DEVELOP_MODE
	bool verbose_log;
#endif
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
template <unsigned SIZE> CKmerBinReader<SIZE>::CKmerBinReader(CKMCParams &Params, CKMCQueues &Queues)
{
	bd = Queues.bd.get();
	bq = Queues.bq.get();
	sorters_manager = Queues.sorters_manager.get();
	tlbq = Queues.tlbq.get();
	disk_logger = Queues.disk_logger.get();

	memory_bins = Queues.memory_bins.get();

	kmer_len       = (uint32)Params.kmer_len;
	cutoff_min     = Params.cutoff_min;
	cutoff_max	   = (uint32)Params.cutoff_max;
	counter_max    = (uint32)Params.counter_max;
	both_strands   = Params.both_strands;
	max_x = Params.max_x;
	s_mapper	   = Queues.s_mapper.get();
	lut_prefix_len = Params.lut_prefix_len;

#ifdef DEVELOP_MODE
	verbose_log = Params.verbose_log;
#endif

	percentProgressObserver = Params.percentProgressObserver;
}

//----------------------------------------------------------------------------------
template <unsigned SIZE> CKmerBinReader<SIZE>::~CKmerBinReader()
{
	
}

//----------------------------------------------------------------------------------
// Read all bins from temporary HDD
template <unsigned SIZE> void CKmerBinReader<SIZE>::ProcessBins()
{
	uchar *data;
	uint64 readed;

	int32 bin_id;
	CMemDiskFile *file;
	string name;
	uint64 size;
	uint64 n_rec;
	uint64 n_plus_x_recs;

	CPercentProgress percent_progress("Stage 2: ", true, percentProgressObserver);
	percent_progress.SetMaxVal(bd->get_n_rec_sum());
	percent_progress.NotifyProgress(0);

	while ((bin_id = bd->get_next_sort_bin()) >= 0)		// Get id of the next bin to read
	{				
		bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs);
		fflush(stdout);

		// Reserve memory necessary to process the current bin at all next stages
		uint64 input_kmer_size;
		uint64 kxmer_counter_size;
		uint32 kxmer_symbols;
		if (max_x)
		{
			input_kmer_size = n_plus_x_recs * sizeof(CKmer<SIZE>);
			kxmer_counter_size = n_plus_x_recs * sizeof(uint32);
			kxmer_symbols = kmer_len + max_x + 1;
		}
		else
		{
			input_kmer_size = n_rec * sizeof(CKmer<SIZE>); 
			kxmer_counter_size = 0;
			kxmer_symbols = kmer_len;
		}
		uint64 max_out_recs    = (n_rec+1) / max(cutoff_min, 1u);

		uint64 counter_size = calc_counter_size(cutoff_max, counter_max);

		uint32 kmer_symbols = kmer_len - lut_prefix_len;
		uint64 kmer_bytes = kmer_symbols / 4;
		if (lut_prefix_len == 0) //do not split data to prefix and sufix (for example when storying result in KFF)
			kmer_bytes = (kmer_symbols + 3) / 4;

		uint64 out_buffer_size = max_out_recs * (kmer_bytes + counter_size);
			
		uint32 rec_len         = (kxmer_symbols + 3) / 4;

		uint64 lut_recs = 1ull << (2 * lut_prefix_len);
		if (lut_prefix_len == 0)
			lut_recs = 0;
		uint64 lut_size = lut_recs * sizeof(uint64);

		// Reserve memory only for the file data
		if (!memory_bins->init(bin_id, rec_len, round_up_to_alignment(size), round_up_to_alignment(input_kmer_size), round_up_to_alignment(out_buffer_size), round_up_to_alignment(kxmer_counter_size), round_up_to_alignment(lut_size)))
		{
			tlbq->insert(bin_id);
			continue;
		}

		// Process the bin if it is not empty
		if(size > 0)
		{
			if (file == nullptr)
			{
				std::ostringstream ostr;
				ostr << "Error: Cannot open temporary file: " << name;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			else
				file->Rewind();

			memory_bins->reserve(bin_id, data, CMemoryBins::mba_input_file);
			//readed = fread(data, 1, size, file);

			readed = file->Read(data, 1, size);
			if(readed != size)
			{
				std::ostringstream ostr;
				ostr << "Error: Corrupted file: " << name << "   " << "Real size : " << readed << "   " << "Should be : " << size;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}

			// Reserve memory necessary to process the whole bin
			memory_bins->extend(bin_id, rec_len, round_up_to_alignment(size), round_up_to_alignment(input_kmer_size), round_up_to_alignment(out_buffer_size), round_up_to_alignment(kxmer_counter_size), round_up_to_alignment(lut_size));
			memory_bins->reserve(bin_id, data, CMemoryBins::mba_input_file);

			// Push bin data to a queue of bins to process			
			bq->push(bin_id, data, size, n_rec);
			sorters_manager->NotifyBQPush();

		}
		else
		{
			//reserve is allowed also for empty bins
			memory_bins->extend(bin_id, rec_len, round_up_to_alignment(size), round_up_to_alignment(input_kmer_size), round_up_to_alignment(out_buffer_size), round_up_to_alignment(kxmer_counter_size), round_up_to_alignment(lut_size));
			// Push empty bin to process (necessary, since all bin ids must be processed)
			bq->push(bin_id, nullptr, 0, 0);
			sorters_manager->NotifyBQPush();
		}

		file->Close();
		
		//Remove temporary file
#ifdef DEVELOP_MODE
		if(!verbose_log)
			file->Remove();	
#else
		file->Remove();	
#endif
		disk_logger->log_remove(size);

		percent_progress.NotifyProgress(n_rec);
	}

	bq->mark_completed();
	sorters_manager->NotifyQueueCompleted();
	fflush(stdout);
}


//************************************************************************************************************
// CWKmerBinReader - wrapper for multithreading purposes
//************************************************************************************************************

//----------------------------------------------------------------------------------
template <unsigned SIZE> class CWKmerBinReader {
	std::unique_ptr<CKmerBinReader<SIZE>> kbr;

public:
	CWKmerBinReader(CKMCParams &Params, CKMCQueues &Queues);

	void operator()();
};

//----------------------------------------------------------------------------------
// Constructor
template <unsigned SIZE> CWKmerBinReader<SIZE>::CWKmerBinReader(CKMCParams &Params, CKMCQueues &Queues)
{
	kbr = std::make_unique<CKmerBinReader<SIZE>>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
template <unsigned SIZE> void CWKmerBinReader<SIZE>::operator()()
{
	kbr->ProcessBins();
}

#endif

// ***** EOF
