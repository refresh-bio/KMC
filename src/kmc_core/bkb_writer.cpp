/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/

#include "bkb_writer.h"
#include "critical_error_handler.h"
#include <sstream>

//************************************************************************************************************
// CBigKmerBinWriter
//************************************************************************************************************

//----------------------------------------------------------------------------------
CBigKmerBinWriter::CBigKmerBinWriter(CKMCParams& Params, CKMCQueues& Queues)
{
	disk_logger = Queues.disk_logger.get();
	bbspq = Queues.bbspq.get();
	sm_pmm_sorter_suffixes = Queues.sm_pmm_sorter_suffixes.get();
	sm_pmm_sorter_lut = Queues.sm_pmm_sorter_lut.get();
	working_directory = Params.working_directory;
	bbd = Queues.bbd.get();
	sm_cbc = Queues.sm_cbc.get();
}

//----------------------------------------------------------------------------------
void CBigKmerBinWriter::Process()
{	
	int32 curr_bin_id = -1;
	uchar* suff_buff = nullptr;
	uint64 suff_buff_size;
	uint64* lut = nullptr;
	uint64 lut_size = 0;
	bool last_one_in_sub_bin;
	bool first_in_sub_bin = true;
	FILE* file = nullptr;
	string name;
	uint64 file_size = 0;
	while (bbspq->pop(bin_id, sub_bin_id, suff_buff, suff_buff_size, lut, lut_size, last_one_in_sub_bin))
	{
		if (curr_bin_id != bin_id)
		{
			if (curr_bin_id != -1)			
				sm_cbc->push(curr_bin_id);			
			curr_bin_id = bin_id;
		}
		if (first_in_sub_bin)
		{
			file_size = 0;			
			name = GetName();
			file = fopen(name.c_str(), "wb+");
			if (!file)
			{
				std::ostringstream ostr;
				ostr << "Error: can not open file : " << name;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			setbuf(file, nullptr);
		}
		first_in_sub_bin = false;

		if (suff_buff_size)
		{
			disk_logger->log_write(suff_buff_size);
			if (fwrite(suff_buff, 1, suff_buff_size, file) != suff_buff_size)
			{
				std::ostringstream ostr;
				ostr << "Error while writing to file : " << name;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			file_size += suff_buff_size;
			sm_pmm_sorter_suffixes->free(suff_buff);
		}

		if (lut_size)
		{
			disk_logger->log_write(lut_size * sizeof(uint64));
			if (fwrite(lut, sizeof(uint64), lut_size, file) != lut_size)
			{
				std::ostringstream ostr;
				ostr << "Error while writing to file : " << name;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			file_size += lut_size * sizeof(uint64);
			sm_pmm_sorter_lut->free((uchar*)lut);
		}
		
		if (last_one_in_sub_bin)
		{			
			bbd->push(bin_id, sub_bin_id, 0, 0, file, name, file_size);
			first_in_sub_bin = true;
		}
	}
	if(curr_bin_id != -1)
		sm_cbc->push(curr_bin_id);
	sm_cbc->mark_completed();
}

//----------------------------------------------------------------------------------
string CBigKmerBinWriter::GetName()
{	
	string s_tmp = std::to_string(bin_id);
	while (s_tmp.length() < 5)
		s_tmp = string("0") + s_tmp;	
	string s1 = std::to_string(sub_bin_id);
	while (s1.length() < 3)
		s1 = string("0") + s1;

	if (*working_directory.rbegin() != '/' && *working_directory.rbegin() != '\\')
		working_directory += "/";
	return working_directory + "kmc_" + s_tmp + "_" + s1 + "_" + s1 + ".bin";
}


//************************************************************************************************************
// CWBigKmerBinWriter - wrapper for multithreading purposes
//************************************************************************************************************


//----------------------------------------------------------------------------------
// Constructor
CWBigKmerBinWriter::CWBigKmerBinWriter(CKMCParams& Params, CKMCQueues& Queues)
{
	bkb_writer = std::make_unique<CBigKmerBinWriter>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
void CWBigKmerBinWriter::operator()()
{
	bkb_writer->Process();
}

// ***** EOF 