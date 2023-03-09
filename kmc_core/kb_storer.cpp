/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.2
  Date   : 2023-03-09
*/
#include <algorithm>
#include <numeric>
#include "kb_storer.h"
#include "critical_error_handler.h"
#include <sstream>

using namespace std;

//************************************************************************************************************
// CKmerBinStorer - storer for bins
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Constructor
CKmerBinStorer::CKmerBinStorer(CKMCParams &Params, CKMCQueues &Queues)
{
	pmm_bins			= Queues.pmm_bins.get();
	n_bins			    = Params.n_bins;
	q_part			    = Queues.bpq.get();
	bd                  = Queues.bd.get();
	epd					= Queues.epd.get();
	working_directory	= Params.working_directory;

	tmp_files_owner		= Queues.tmp_files_owner.get();

	s_mapper			= Queues.s_mapper.get();
	disk_logger			= Queues.disk_logger.get();
	buffer_size_bytes      = 0;
	max_buf_size		   = 0;
	max_buf_size_id		   = 0;
	max_mem_buffer         = Params.max_mem_storer;

	max_mem_single_package = Params.max_mem_storer_pkg;
	tmp_buff = std::unique_ptr<uchar[]>(new uchar[max_mem_single_package*2]); //no std::make_unique<uchar[]>(max_mem_single_package*2), because it clears memory which I don't want here
	
	buffer.resize(n_bins);

	total_size = 0 ; 

}

//----------------------------------------------------------------------------------
// Destructor
CKmerBinStorer::~CKmerBinStorer()
{
	Release();
}

//----------------------------------------------------------------------------------
// Write ends of bins and release memory
void CKmerBinStorer::Release()
{
	buffer.clear();
	buf_sizes.clear();
	tmp_buff.reset();
}

//----------------------------------------------------------------------------------
// Put buffer items to the queue
void CKmerBinStorer::ReleaseBuffer()
{
	for(int i = 0; i < n_bins; ++i)
		if(buffer[i])
			PutBinToTmpFile(i);

	buffer.clear();
}

//----------------------------------------------------------------------------------
// Return name of a file related to a kmer of given id.
string CKmerBinStorer::GetName(int n)
{
	string s_tmp = std::to_string(n);
	while(s_tmp.length() < 5)
		s_tmp = string("0") + s_tmp;
	
	if (*working_directory.rbegin() != '/' && *working_directory.rbegin() != '\\')
		working_directory += "/";
	return working_directory + "kmc_" + s_tmp + ".bin";
}

//----------------------------------------------------------------------------------
// Check wheter it is necessary to store some bin to a HDD
void CKmerBinStorer::CheckBuffer()
{
	int32 i;

	if(buffer_size_bytes < max_mem_buffer && max_buf_size < max_mem_single_package)
		return;

	PutBinToTmpFile(max_buf_size_id);

	buf_sizes[max_buf_size_id] = 0;

	max_buf_size    = buf_sizes[0];
	max_buf_size_id = 0;
	for(i = 1; i < n_bins; ++i)
	{
		if(buf_sizes[i] > max_buf_size)
		{
			max_buf_size    = buf_sizes[i];
			max_buf_size_id = i;
		}

	}
}

//----------------------------------------------------------------------------------
// Send bin to temp file
void CKmerBinStorer::PutBinToTmpFile(uint32 n)
{
	if(buf_sizes[n])
	{
		uint64 w;
		uint64 tmp_buff_pos = 0;
		uint32 size;
		uchar* buf;
		for(auto p = buffer[n]->begin() ; p != buffer[n]->end() ; ++p)
		{
			buf = get<0>(*p);
			size = get<1>(*p);
			memcpy(tmp_buff.get() + tmp_buff_pos, buf, size);
			tmp_buff_pos += size;
			pmm_bins->free(buf);
		}

		disk_logger->log_write(tmp_buff_pos);
		w = tmp_files_owner->Get(n)->Write(tmp_buff.get(), 1, tmp_buff_pos);
		if(w != tmp_buff_pos)
		{
			std::ostringstream ostr;
			ostr << "Error while writing to temporary file " << n;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		total_size += w;		
		buffer_size_bytes -= buf_sizes[n];
	}
	buffer[n]->clear();
}
//


//----------------------------------------------------------------------------------
// Open temporary files for all bins
bool CKmerBinStorer::OpenFiles()
{
	string f_name;

	tmp_files_owner->CreateInstances();

	buf_sizes.resize(n_bins);

	for(int i = 0; i < n_bins; ++i)
	{
		f_name = GetName(i);
		buf_sizes[i] = 0;

		auto file = tmp_files_owner->Get(i);
		file->Open(f_name);

		bd->insert(i, file, f_name);
	}

	return true;
}

//----------------------------------------------------------------------------------
// 
void CKmerBinStorer::ProcessQueue()
{
	// Process the queue
	while(!q_part->completed())
	{
		int32 bin_id;
		uchar *part;
		uint32 true_size;
		uint32 alloc_size;

		list<pair<uint64, uint64>> expander_parts;
		if (q_part->pop(bin_id, part, true_size, alloc_size, expander_parts))
		{
			epd->push(bin_id, expander_parts);
			expander_parts.clear();
			if(!buffer[bin_id])
				buffer[bin_id] = std::make_unique<elem_t>();
			buffer[bin_id]->push_back(make_tuple(part, true_size, alloc_size));
			buffer_size_bytes += alloc_size;
			buf_sizes[bin_id] += alloc_size;

			if(buf_sizes[bin_id] > max_buf_size)
			{
				max_buf_size    = buf_sizes[bin_id];
				max_buf_size_id = bin_id;
			}

			CheckBuffer();
		}
	}

	// Move all remaining parts to queue
	ReleaseBuffer();


}


//************************************************************************************************************
// CWKmerBinStorer - wrapper
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Constructor
CWKmerBinStorer::CWKmerBinStorer(CKMCParams &Params, CKMCQueues &Queues)
{
	kbs = std::make_unique<CKmerBinStorer>(Params, Queues);
	kbs->OpenFiles();
}

//----------------------------------------------------------------------------------
// Execution
void CWKmerBinStorer::operator()()
{
	kbs->ProcessQueue();
}

// ***** EOF
