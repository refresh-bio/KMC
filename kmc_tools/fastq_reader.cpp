/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"

#include <algorithm>
#include <cstring>

#include "defs.h"
#include "fastq_reader.h"

//************************************************************************************************************
// CFastqReader	- reader class
//************************************************************************************************************

uint64 CFastqReader::OVERHEAD_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
// Constructor of FASTA/FASTQ reader
// Parameters:
//    * _mm - pointer to memory monitor (to check the memory limits)
CFastqReader::CFastqReader(CMemoryPool *_pmm_fastq, CFilteringParams::file_type _file_type, uint32 _gzip_buffer_size, uint32 _bzip2_buffer_size, int _kmer_len)
{
	pmm_fastq = _pmm_fastq;

	file_type  = _file_type;
	kmer_len = _kmer_len;
	// Input file mode (default: uncompressed)
	mode      = m_plain;

	// Pointers to input files in various formats (uncompressed, gzip-compressed, bzip2-compressed)
	in		  = NULL;
	in_gzip   = NULL;
	in_bzip2  = NULL;
	bzerror   = BZ_OK;

	// Size and pointer for the buffer
	part_size = 1 << 23;
	part      = NULL;

	gzip_buffer_size  = _gzip_buffer_size;
	bzip2_buffer_size = _bzip2_buffer_size;


}

//----------------------------------------------------------------------------------
// Destructor - close the files
CFastqReader::~CFastqReader()
{
	if(mode == m_plain)
	{
		if(in)
			fclose(in);
	}
	else if(mode == m_gzip)
	{
		if(in_gzip)
			gzclose(in_gzip);
	}
	else if(mode == m_bzip2)
	{
		if(in)
		{
			BZ2_bzReadClose(&bzerror, in_bzip2);
			fclose(in);
		}
	}

	if(part)
		pmm_fastq->free(part);
}

//----------------------------------------------------------------------------------
// Set the name of the file to process
bool CFastqReader::SetNames(string _input_file_name)
{
	input_file_name = _input_file_name;

	// Set mode according to the extension of the file name
	if(input_file_name.size() > 3 && string(input_file_name.end()-3, input_file_name.end()) == ".gz")
		mode = m_gzip;
	else if(input_file_name.size() > 4 && string(input_file_name.end()-4, input_file_name.end()) == ".bz2")
		mode = m_bzip2;
	else
		mode = m_plain;

	return true;
}

//----------------------------------------------------------------------------------
// Set part size of the buffer
bool CFastqReader::SetPartSize(uint64 _part_size)
{
	if(in || in_gzip || in_bzip2)
		return false;

	if(_part_size < (1 << 20) || _part_size > (1 << 30))
		return false;

	part_size = _part_size;

	return true;
}

//----------------------------------------------------------------------------------
// Open the file
bool CFastqReader::OpenFiles()
{
	if(in || in_gzip || in_bzip2)
		return false;

	// Uncompressed file
	if(mode == m_plain)	
	{
		if((in = fopen(input_file_name.c_str(), "rb")) == NULL)
			return false;
	}
	// Gzip-compressed file
	else if(mode == m_gzip)
	{
		if((in_gzip = gzopen(input_file_name.c_str(), "rb")) == NULL)
			return false;
		gzbuffer(in_gzip, gzip_buffer_size);
	}
	// Bzip2-compressed file
	else if(mode == m_bzip2)
	{
		in = fopen(input_file_name.c_str(), "rb");
		if(!in)
			return false;
		setvbuf(in, NULL, _IOFBF, bzip2_buffer_size);
		if((in_bzip2 = BZ2_bzReadOpen(&bzerror, in, 0, 0, NULL, 0)) == NULL)
		{
			fclose(in);
			return false;
		}
	}
	
	// Reserve via PMM
	pmm_fastq->reserve(part);

	part_filled = 0;

	return true;
}

//----------------------------------------------------------------------------------
// Read a part of the file
bool CFastqReader::GetPart(uchar *&_part, uint64 &_size)
{	
	if(!in && !in_gzip && !in_bzip2)
		return false;


	if(IsEof())
		return false;
	uint64 readed;
	
	// Read data
	if(mode == m_plain)
		readed = fread(part+part_filled, 1, part_size - part_filled, in);
	else if(mode == m_gzip)
		readed = gzread(in_gzip, part+part_filled, (int) (part_size - part_filled));
	else if(mode == m_bzip2)
		readed = BZ2_bzRead(&bzerror, in_bzip2, part+part_filled, (int) (part_size - part_filled));
	else
		readed = 0;				// Never should be here

	int64 total_filled = part_filled + readed;
	int64 i;

	if(IsEof())
	{
		_part = part;
		_size = total_filled;

		part = NULL;
		return true;
	}
	
	// Look for the end of the last complete record in a buffer
	if(file_type == CFilteringParams::file_type::fasta)			// FASTA files
	{
		// Looking for a FASTA record at the end of the area
		i = total_filled - 1;
		int64 start, end;
		int64 line_start[4], line_end[4];
		int readed_lines = 0;
		bool success = false;
		int k;
		while (i >= 0 && readed_lines < 4)
		{
			GetFullLineFromEnd(start, end, part, i);

			line_start[4 - readed_lines - 1] = start;
			line_end[4 - readed_lines - 1] = end;
			++readed_lines;

			if (readed_lines >= 2)
			{
				k = 4 - readed_lines;
				if (part[line_start[k]] == '>')
				{
					success = true;
					break;

				}
			}
		}
		// Looking for a FASTQ record at the end of the area
		if (!success)
		{
			cerr << "Error: Wrong input file!\n";
			exit(1);
		}

		_part = part;
		_size = line_end[k + 1];
	}
	else
	{
		i = total_filled - 1;
		int64 start, end;
		int64 line_start[8], line_end[8];
		int readed_lines = 0;
		bool success = false;
		int k;
		while (i >= 0 && readed_lines < 8)
		{
			GetFullLineFromEnd(start, end, part, i);
			line_start[8 - readed_lines - 1] = start;
			line_end[8 - readed_lines - 1] = end;
			++readed_lines;

			if (readed_lines >= 4)
			{
				k = 8 - readed_lines;
				if (part[line_start[k]] == '@' && part[line_start[k + 2]] == '+')
				{
					if (part[line_start[k + 2] + 1] == '\n' || part[line_start[k + 2] + 1] == '\r')
					{
						success = true;
						break;
					}
					if (line_start[k + 1] - line_start[k] == line_start[k + 3] - line_start[k + 2] &&
						memcmp(part + line_start[k] + 1, part + line_start[k + 2] + 1, line_start[k + 3] - line_start[k + 2] - 1) == 0)
					{
						success = true;
						break;
					}
				}
			}
		}

		if (!success)
		{
			cerr << "Error: Wrong input file!\n";
			exit(1);
		}
		_part = part;
		_size = line_end[k + 3];
	}
	// Allocate new memory for the buffer

	pmm_fastq->reserve(part);
	copy(_part+_size, _part+total_filled, part);
	part_filled = total_filled - _size;

	return true;
}

//----------------------------------------------------------------------------------
// Skip to next EOL from the current position in a buffer
bool CFastqReader::SkipNextEOL(uchar *part, int64 &pos, int64 max_pos)
{
	int64 i;
	for(i = pos; i < max_pos-2; ++i)
		if((part[i] == '\n' || part[i] == '\r') && !(part[i+1] == '\n' || part[i+1] == '\r'))
			break;

	if(i >= max_pos-2)
		return false;

	pos = i+1;

	return true;
}

void CFastqReader::GetFullLineFromEnd(int64& line_sart, int64& line_end, uchar* buff, int64& pos)
{
	while (pos >= 0 && buff[pos] != '\n' && buff[pos] != '\r')
		--pos;
	line_end = pos + 1;
	while (pos >= 0 && (buff[pos] == '\n' || buff[pos] == '\r'))
		--pos;
	while (pos >= 0 && buff[pos] != '\n' && buff[pos] != '\r')
		--pos;
	line_sart = pos + 1;
};
//----------------------------------------------------------------------------------
// Check whether there is an EOF
bool CFastqReader::IsEof()
{
	if(mode == m_plain)
		return feof(in) != 0;
	else if(mode == m_gzip)
		return gzeof(in_gzip) != 0;
	else if(mode == m_bzip2)
		return bzerror == BZ_STREAM_END;

	return true;
}



//************************************************************************************************************
// CWFastqReader - wrapper for multithreading purposes
//************************************************************************************************************
CWFastqReader::CWFastqReader(CFilteringParams &Params, CFilteringQueues &Queues)
{	
	pmm_fastq = Queues.pmm_fastq_reader;

	input_files_queue = Queues.input_files_queue;
	part_size		  = Params.fastq_buffer_size;
	part_queue		  = Queues.input_part_queue;
	file_type         = Params.input_file_type;
	kmer_len		  = Params.kmer_len;

	gzip_buffer_size  = Params.gzip_buffer_size;
	bzip2_buffer_size = Params.bzip2_buffer_size;

	fqr = nullptr;
}

//----------------------------------------------------------------------------------
CWFastqReader::~CWFastqReader()
{
}

//----------------------------------------------------------------------------------
void CWFastqReader::operator()()
{
	uchar *part;
	uint64 part_filled;
	
	while(input_files_queue->pop(file_name))
	{
		fqr = new CFastqReader(pmm_fastq, file_type, gzip_buffer_size, bzip2_buffer_size, kmer_len);
		fqr->SetNames(file_name);
		fqr->SetPartSize(part_size);

		if(fqr->OpenFiles())
		{
			// Reading Fastq parts
			while(fqr->GetPart(part, part_filled))
				part_queue->push(part, part_filled);
		}
		else
			cerr << "Error: Cannot open file " << file_name << "\n";
		delete fqr;
	}
	part_queue->mark_completed();
}

// ***** EOF
