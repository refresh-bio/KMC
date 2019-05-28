#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/
#include <algorithm>
#include <numeric>
#include <iostream>
#include "kb_storer.h"

using namespace std;

extern uint64 total_reads;

//************************************************************************************************************
// CKmerBinStorer - storer for bins
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Constructor
CKmerBinStorer::CKmerBinStorer(CKMCParams &Params, CKMCQueues &Queues)
{
	pmm_bins			= Queues.pmm_bins;
	mm					= Queues.mm;
	n_bins			    = Params.n_bins;
	q_part			    = Queues.bpq;
	bd                  = Queues.bd;
	epd					= Queues.epd;
	working_directory = Params.working_directory;

	mem_mode			= Params.mem_mode;

	s_mapper			= Queues.s_mapper;
	disk_logger			= Queues.disk_logger;
	files                  = nullptr;
	buf_sizes		       = nullptr;
	buffer_size_bytes      = 0;
	max_buf_size		   = 0;
	max_buf_size_id		   = 0;
	max_mem_buffer         = Params.max_mem_storer;

	max_mem_single_package = Params.max_mem_storer_pkg;
	tmp_buff = new uchar[max_mem_single_package*2]; 
	
	buffer = new elem_t*[n_bins];
	for(int i = 0; i < n_bins; ++i)
		buffer[i] = nullptr;


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
	if(!files)
		return;
	for(int i = 0; i < n_bins; ++i)
		if(buffer[i])
			delete buffer[i];

	delete[] buffer;
	buffer = nullptr;

	delete[] files;
	files = nullptr;

	delete[] buf_sizes;
	buf_sizes = nullptr;
	
	delete [] tmp_buff;

	cerr << "\n";
}

//----------------------------------------------------------------------------------
// Put buffer items to the queue
void CKmerBinStorer::ReleaseBuffer()
{
	for(int i = 0; i < n_bins; ++i)
		if(buffer[i])
			PutBinToTmpFile(i);

	for(int i = n_bins-1; i >= 0; --i)
		if(buffer[i])
		{
			delete buffer[i];
			buffer[i] = nullptr;
		}
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
			memcpy(tmp_buff + tmp_buff_pos, buf, size);
			tmp_buff_pos += size;
			pmm_bins->free(buf);
		}

		disk_logger->log_write(tmp_buff_pos);
		w = files[n]->Write(tmp_buff, 1, tmp_buff_pos);
		if(w != tmp_buff_pos)
		{
			cerr << "Error while writing to temporary file " << n;
			exit(1);
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

	files     = new CMemDiskFile*[n_bins];
	for (int i = 0 ; i < n_bins ; ++i)
	{
		files[i] = new CMemDiskFile(mem_mode);
	}
	buf_sizes = new uint64[n_bins];

	for(int i = 0; i < n_bins; ++i)
	{
		f_name = GetName(i);
		buf_sizes[i] = 0;
		
		files[i]->Open(f_name);

		bd->insert(i, files[i], f_name, 0, 0, 0, 0);
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
				buffer[bin_id] = new elem_t;
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
	kbs = new CKmerBinStorer(Params, Queues);
	kbs->OpenFiles();
}

//----------------------------------------------------------------------------------
// Destructore
CWKmerBinStorer::~CWKmerBinStorer()
{
	delete kbs;
}

//----------------------------------------------------------------------------------
// Execution
void CWKmerBinStorer::operator()()
{
	kbs->ProcessQueue();
}

// ***** EOF
