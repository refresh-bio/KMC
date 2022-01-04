/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.1
  Date   : 2022-01-04
*/

#include "bkb_reader.h"


//************************************************************************************************************
// CBigKmerBinReader 
//************************************************************************************************************

//----------------------------------------------------------------------------------
CBigKmerBinReader::CBigKmerBinReader(CKMCParams& Params, CKMCQueues& Queues)
{
	tlbq = Queues.tlbq.get();
	disk_logger = Queues.disk_logger.get();
	bd   = Queues.bd.get();
	bbpq = Queues.bbpq.get();
	sm_pmm_input_file = Queues.sm_pmm_input_file.get();
	sm_mem_part_input_file = Params.sm_mem_part_input_file;
	progressObserver = Params.progressObserver;

	kmer_len = (uint32)Params.kmer_len;
}

//----------------------------------------------------------------------------------
void CBigKmerBinReader::ProcessBigBin()
{
	int32 bin_id;
	CMemDiskFile *file;
	string name;
	uint64 size, n_rec, n_plus_x_recs, in_buffer, end_pos;

	uchar *file_buff, *tmp;
	
	progressObserver->Start("Big bins");
	while (tlbq->get_next(bin_id))
	{
		bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs);
		progressObserver->Step();
		file->Rewind();
		end_pos = 0;
		sm_pmm_input_file->reserve(file_buff);
		while ( (in_buffer = end_pos + file->Read(file_buff + end_pos, 1, sm_mem_part_input_file - end_pos)) )
		{
			end_pos = 0;
			for (; end_pos + 1 + (file_buff[end_pos] + kmer_len + 3) / 4 <= in_buffer; end_pos += 1 + (file_buff[end_pos] + kmer_len + 3) / 4);
			uint64 rest = in_buffer - end_pos;
			sm_pmm_input_file->reserve(tmp);
			memcpy(tmp, file_buff + end_pos, rest);
			bbpq->push(bin_id, file_buff, end_pos);
			file_buff = tmp;
			end_pos = rest;
		}
		sm_pmm_input_file->free(file_buff);
		file->Close();

		//Remove file
		file->Remove();
		disk_logger->log_remove(size);
	}
	bbpq->mark_completed();
	progressObserver->End();
}

//----------------------------------------------------------------------------------
CBigKmerBinReader::~CBigKmerBinReader()
{

}

//************************************************************************************************************
// CWBigKmerBinReader - wrapper for multithreading purposes
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Constructor
CWBigKmerBinReader::CWBigKmerBinReader(CKMCParams& Params, CKMCQueues& Queues)
{
	bkb_reader = std::make_unique<CBigKmerBinReader>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
void CWBigKmerBinReader::operator()()
{
	bkb_reader->ProcessBigBin();
}

// ***** EOF