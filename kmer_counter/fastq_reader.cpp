#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#include <algorithm>
#include "defs.h"
#include "fastq_reader.h"
#include "asmlib_wrapper.h"
//************************************************************************************************************
// CFastqReader	- reader class
//************************************************************************************************************

uint64 CFastqReader::OVERHEAD_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
// Constructor of FASTA/FASTQ reader
// Parameters:
//    * _mm - pointer to memory monitor (to check the memory limits)
CFastqReader::CFastqReader(CMemoryMonitor *_mm, CMemoryPool *_pmm_fastq, input_type _file_type, int _kmer_len, CBinaryPackQueue* _binary_pack_queue, CMemoryPool* _pmm_binary_file_reader)
{
	binary_pack_queue = _binary_pack_queue;

	mm		  = _mm;
	pmm_fastq = _pmm_fastq;
	pmm_binary_file_reader = _pmm_binary_file_reader;

	file_type  = _file_type;
	kmer_len = _kmer_len;
	
	// Size and pointer for the buffer
	part_size = 1 << 23;
	part      = NULL;

	containsNextChromosome = false;

	data_src.SetQueue(binary_pack_queue, pmm_binary_file_reader);
}

//----------------------------------------------------------------------------------
// Destructor - close the files
CFastqReader::~CFastqReader()
{
	if(part)
		pmm_fastq->free(part);
}

//----------------------------------------------------------------------------------
// Set part size of the buffer
bool CFastqReader::SetPartSize(uint64 _part_size)
{
	if(_part_size < (1 << 20) || _part_size > (1 << 30))
		return false;

	part_size = _part_size;

	return true;
}


//----------------------------------------------------------------------------------
// Read a part of the file in multi line fasta format
bool CFastqReader::GetPartFromMultilneFasta(uchar *&_part, uint64 &_size)
{
		uint64 readed = 0;

		if(!containsNextChromosome)
		{
			if (data_src.Finished())
				return false;
		}

		readed = data_src.read(part + part_filled, part_size - part_filled);
		
		int64 total_filled = part_filled + readed;
		int64 last_header_pos = 0;
		int64 pos = 0;
		for(int64 i = 0 ; i < total_filled ;++i )//find last '>' and remove EOLs
		{
			if(part[i] == '>')
			{
				int64 tmp = i;
				bool next_line = SkipNextEOL(part, i, total_filled);
				if (!next_line)
					i = total_filled;
				copy(part+tmp, part+i, part+pos);
				last_header_pos = pos;
				pos += i - tmp;
				if (!next_line)
					break;
			}
			if(part[i] != '\n' && part[i] != '\r')
			{
				part[pos++] = part[i];
			}
 		}

		_part = part;
		if(last_header_pos == 0)//data in block belong to one seq
		{
			part_filled = kmer_len - 1;
			_size = pos;
			pmm_fastq->reserve(part);
			copy(_part+_size-part_filled, _part+_size, part);
			containsNextChromosome = false;
		}
		else//next seq starts at last_header_pos
		{
			_size = last_header_pos;
			part_filled = pos - last_header_pos;
			pmm_fastq->reserve(part);
			copy(_part + last_header_pos, _part + pos, part);
			containsNextChromosome = true;
		}
		return true;
}

bool CFastqReader::GetPartNew(uchar *&_part, uint64 &_size)
{
	if (file_type == multiline_fasta)
		return GetPartFromMultilneFasta(_part, _size);

	if (data_src.Finished())
		return false;
	uint64 readed;

	// Read data
	readed = data_src.read(part + part_filled, part_size);

	//std::cerr.write((char*)part + part_filled, readed);

	int64 total_filled = part_filled + readed;
	int64 i;

	if (part_filled >= OVERHEAD_SIZE)
	{
		cerr << "Error: Wrong input file!\n";
		exit(1);
	}

	if (data_src.Finished())
	{
		_part = part;
		_size = total_filled;

		part = NULL;
		return true;
	}

	// Look for the end of the last complete record in a buffer
	if (file_type == fasta)			// FASTA files
	{
		// Looking for a FASTA record at the end of the area
		int64 line_start[3];
		int32 j;

		i = total_filled - OVERHEAD_SIZE / 2;
		for (j = 0; j < 3; ++j)
		{
			if (!SkipNextEOL(part, i, total_filled))
				break;
			line_start[j] = i;
		}

		_part = part;
		if (j < 3)
			_size = 0;
		else
		{
			int k;
			for (k = 0; k < 2; ++k)
			if (part[line_start[k] + 0] == '>')
				break;

			if (k == 2)
				_size = 0;
			else
				_size = line_start[k];
		}
	}
	else			// FASTQ file
	{
		// Looking for a FASTQ record at the end of the area
		int64 line_start[9];
		int32 j;

		i = total_filled - OVERHEAD_SIZE / 2;
		for (j = 0; j < 9; ++j)
		{
			if (!SkipNextEOL(part, i, total_filled))
				break;
			line_start[j] = i;
		}

		_part = part;
		if (j < 9)
			_size = 0;
		else
		{
			int k;
			for (k = 0; k < 4; ++k)
			{
				if (part[line_start[k] + 0] == '@' && part[line_start[k + 2] + 0] == '+')
				{
					if (part[line_start[k + 2] + 1] == '\n' || part[line_start[k + 2] + 1] == '\r')
						break;
					if (line_start[k + 1] - line_start[k] == line_start[k + 3] - line_start[k + 2] &&
						memcmp(part + line_start[k] + 1, part + line_start[k + 2] + 1, line_start[k + 3] - line_start[k + 2] - 1) == 0)
						break;
				}
			}

			if (k == 4)
				_size = 0;
			else
				_size = line_start[k];
		}
	}
	// Allocate new memory for the buffer

	pmm_fastq->reserve(part);
	copy(_part + _size, _part + total_filled, part);
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




//************************************************************************************************************
// CFastqReaderDataSrc - data source for FastqReader
//************************************************************************************************************
void CFastqReaderDataSrc::init_stream()
{
	switch (compression_type)
	{
	case CompressionType::plain: 
		in_data_pos = 0;
		break;
	case CompressionType::gzip:
		stream.zalloc = Z_NULL;
		stream.zfree = Z_NULL;
		stream.opaque = Z_NULL;
		stream.avail_in = 0;
		stream.next_in = Z_NULL;
		if (inflateInit2(&stream, 31) != Z_OK)
		{
			cerr << "Error while reading gz file\n";
			exit(1);
		}
		stream.avail_in = (uint32)in_data_size;
		stream.next_in = in_data;
		break;
	case CompressionType::bzip2:
		_bz_stram.bzalloc = NULL;
		_bz_stram.bzfree = NULL;
		_bz_stram.opaque = NULL;
		_bz_stram.avail_in = 0;
		_bz_stram.next_in = NULL;
		if (BZ2_bzDecompressInit(&_bz_stram, 0, 0) != BZ_OK)
		{
			cerr << "Error while reading bz2 file\n";
			exit(1);
		}
		_bz_stram.avail_in = (uint32)in_data_size;
		_bz_stram.next_in = (char*)in_data;
		break;
	default:
		break;
	}
}

//----------------------------------------------------------------------------------
void CFastqReaderDataSrc::SetQueue(CBinaryPackQueue* _binary_pack_queue, CMemoryPool *_pmm_binary_file_reader)
{
	binary_pack_queue = _binary_pack_queue;
	pmm_binary_file_reader = _pmm_binary_file_reader;
}

//----------------------------------------------------------------------------------
bool CFastqReaderDataSrc::Finished()
{
	return end_reached;
}

//----------------------------------------------------------------------------------
uint64 CFastqReaderDataSrc::read(uchar* buff, uint64 size)
{
	if (!in_progress)
	{
		if (!binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type))
		{
			end_reached = true;
			return 0;
		}
		in_progress = true;		
		init_stream();
	}

	//reading
	if (compression_type == CompressionType::gzip)
	{
		stream.next_out = buff;
		stream.avail_out = (uint32)size;
		int ret;
		do
		{
			if (!stream.avail_in)
			{
				pmm_binary_file_reader->free(in_data);
				in_data = nullptr;
				if(!binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type))
					return 0;
				stream.avail_in = (uint32)in_data_size;
				stream.next_in = in_data;				
			}

			ret = inflate(&stream, Z_NO_FLUSH);

			switch (ret)
			{
			case Z_NEED_DICT:
				ret = Z_DATA_ERROR;     /* and fall through */
			case Z_DATA_ERROR:
				cerr << "Some error while reading gzip file\n"; 
				exit(1);
			case Z_MEM_ERROR:
				inflateEnd(&stream);
				return ret;
			}

			if (ret == Z_STREAM_END)
			{				
				bool multistream = stream.avail_in || !binary_pack_queue->is_next_last();
				if (!multistream)
				{
					pmm_binary_file_reader->free(in_data);
					in_data = nullptr;
					inflateEnd(&stream);
					in_progress = false;
					//pull end
					bool queue_end = !binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type);
					if (!queue_end && file_part != FilePart::End)
					{
						cerr << "Error: An internal error occurred. Please contact authors\n";
					}
					break;
				}
				else //multiple streams in one file
				{
					//equivalent of inflateReset 
				/*	inflateEnd(&stream);
					if (inflateInit2(&stream, 31) != Z_OK) 
					{
						cerr << "Error while reading gzip file\n";
						exit(1);
					}*/
					if (inflateReset(&stream) != Z_OK)
					{
						cerr << "Error while reading gzip file\n";
						exit(1);
					}
				}
			}
		} while (stream.avail_out);
		return size - stream.avail_out;
	}
	else if (compression_type == CompressionType::bzip2)
	{
		_bz_stram.next_out = (char*)buff;
		_bz_stram.avail_out = (uint32)size;
		int ret;
		do
		{
			if (!_bz_stram.avail_in)
			{
				pmm_binary_file_reader->free(in_data);
				in_data = nullptr;
				binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type);
				_bz_stram.avail_in = (uint32)in_data_size;
				_bz_stram.next_in = (char*)in_data;				
			}
			ret = BZ2_bzDecompress(&_bz_stram);
			if (ret == BZ_PARAM_ERROR || ret == BZ_DATA_ERROR || ret == BZ_DATA_ERROR_MAGIC || ret == BZ_MEM_ERROR)
			{
				BZ2_bzDecompressEnd(&_bz_stram);
				cerr << "bz2 reading error\n"; 
			}
			if (ret == BZ_STREAM_END)
			{
				bool multistream = _bz_stram.avail_in || !binary_pack_queue->is_next_last();
				if (!multistream)
				{
					pmm_binary_file_reader->free(in_data);
					in_data = nullptr;
					BZ2_bzDecompressEnd(&_bz_stram);
					in_progress = false;
					//pull end
					bool queue_end = !binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type);
					if (!queue_end && file_part != FilePart::End)
					{
						cerr << "Error: An internal error occurred. Please contact authors\n";
					}
					break;
				}
				else
				{
					BZ2_bzDecompressEnd(&_bz_stram);
					if (BZ2_bzDecompressInit(&_bz_stram, 0, 0) != BZ_OK)
					{
						cerr << "Error while reading bz2 file\n";
						exit(1);
					}
				}
			}

		} while (_bz_stram.avail_out);
		return size - _bz_stram.avail_out;
	}
	else if (compression_type == CompressionType::plain) 
	{
		uint64 out_pos = 0;
		do
		{
			if (in_data_pos >= in_data_size)
			{
				pmm_binary_file_reader->free(in_data);
				in_data = nullptr;
				binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type);
				if (file_part == FilePart::End)
				{
					in_progress = false;
					break;
				}
				in_data_pos = 0;
			}
			uint64 in_left = in_data_size - in_data_pos;
			uint64 out_left = size - out_pos;
			uint64 n_copy = min(in_left, out_left);
			A_memcpy(buff + out_pos, in_data + in_data_pos, n_copy);
			in_data_pos += n_copy;
			out_pos += n_copy;
		} while (out_pos < size);
		return out_pos;
	}
	else
	{
		cerr << "Error: unknown compression\n";
		exit(1);
	}
}




//************************************************************************************************************
// CWFastqReader - wrapper for multithreading purposes
//************************************************************************************************************
CWFastqReader::CWFastqReader(CKMCParams &Params, CKMCQueues &Queues, CBinaryPackQueue* _binary_pack_queue)
{
	mm = Queues.mm;
	pmm_fastq = Queues.pmm_fastq;
	pmm_binary_file_reader = Queues.pmm_binary_file_reader;
	//input_files_queue = Queues.input_files_queue;
	binary_pack_queue = _binary_pack_queue;
	part_size		  = Params.fastq_buffer_size;
	part_queue		  = Queues.part_queue;
	file_type         = Params.file_type;
	kmer_len		  = Params.p_k;


	fqr = NULL;
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

	fqr = new CFastqReader(mm, pmm_fastq, file_type, kmer_len, binary_pack_queue, pmm_binary_file_reader);
	fqr->SetPartSize(part_size);
	fqr->Init();
	while (fqr->GetPartNew(part, part_filled))
		part_queue->push(part, part_filled);
	
	delete fqr;
	part_queue->mark_completed();
}



//************************************************************************************************************
// CWStatsFastqReader - wrapper for multithreading purposes
//************************************************************************************************************
CWStatsFastqReader::CWStatsFastqReader(CKMCParams &Params, CKMCQueues &Queues, CBinaryPackQueue* _binary_pack_queue)
{
	mm = Queues.mm;
	pmm_fastq = Queues.pmm_fastq;
	pmm_binary_file_reader = Queues.pmm_binary_file_reader;

	binary_pack_queue = _binary_pack_queue;
	
	part_size = Params.fastq_buffer_size;
	stats_part_queue = Queues.stats_part_queue;
	file_type = Params.file_type;
	kmer_len = Params.p_k;

	fqr = NULL;
}

//----------------------------------------------------------------------------------
CWStatsFastqReader::~CWStatsFastqReader()
{
}

//----------------------------------------------------------------------------------
void CWStatsFastqReader::operator()()
{
	uchar *part;
	uint64 part_filled;

	fqr = new CFastqReader(mm, pmm_fastq, file_type, kmer_len, binary_pack_queue, pmm_binary_file_reader);
	fqr->SetPartSize(part_size);
	fqr->Init();
	bool finished = false;
	while (fqr->GetPartNew(part, part_filled) && !finished)
	{
		if (!stats_part_queue->push(part, part_filled))
		{
			finished = true;
			pmm_fastq->free(part);
			binary_pack_queue->ignore_rest();
			fqr->IgnoreRest();			
			break;
		}
	}

	delete fqr;
	stats_part_queue->mark_completed();
}


// ***** EOF
