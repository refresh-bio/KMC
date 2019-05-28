/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#include "stdafx.h"
#include "bkb_reader.h"


//************************************************************************************************************
// CBigKmerBinReader 
//************************************************************************************************************

//----------------------------------------------------------------------------------
CBigKmerBinReader::CBigKmerBinReader(CKMCParams& Params, CKMCQueues& Queues)
{
	tlbq = Queues.tlbq;
	disk_logger = Queues.disk_logger;
	bd   = Queues.bd;
	bbpq = Queues.bbpq;
	sm_pmm_input_file = Queues.sm_pmm_input_file;
	sm_mem_part_input_file = Params.sm_mem_part_input_file;
}

//----------------------------------------------------------------------------------
void CBigKmerBinReader::ProcessBigBin()
{
	int32 bin_id;
	CMemDiskFile *file;
	string name;
	uint64 size, n_rec, n_plus_x_recs, in_buffer, end_pos;
	uint32 buffer_size, kmer_len;		
	uchar *file_buff, *tmp;
	
	while (tlbq->get_next(bin_id))
	{
		bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, buffer_size, kmer_len);
		cerr << "*";
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
	bkb_reader = new CBigKmerBinReader(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
CWBigKmerBinReader::~CWBigKmerBinReader()
{
	delete bkb_reader;
}

//----------------------------------------------------------------------------------
// Execution
void CWBigKmerBinReader::operator()()
{
	bkb_reader->ProcessBigBin();
}

// ***** EOF