/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.4
  Date   : 2024-02-09
*/
#include <algorithm>
#include "defs.h"
#include "fastq_reader.h"
#include "bam_utils.h"
#include "critical_error_handler.h"
#include <sstream>

#include "bkb_sorter.h"
//************************************************************************************************************
// CFastqReader	- reader class
//************************************************************************************************************

uint64 CFastqReader::OVERHEAD_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
// Constructor of FASTA/FASTQ reader
// Parameters:
//    * _mm - pointer to memory monitor (to check the memory limits)
CFastqReader::CFastqReader(CMemoryPoolWithBamSupport *_pmm_fastq, InputType _file_type, int _kmer_len,
	CBinaryPackQueue* _binary_pack_queue, CMemoryPool* _pmm_binary_file_reader, CBamTaskManager* _bam_task_manager, 
	CPartQueue* _part_queue, CPartQueue* _stats_part_queue, CMissingEOL_at_EOF_counter* _missingEOL_at_EOF_counter)
{
	binary_pack_queue = _binary_pack_queue;
	missingEOL_at_EOF_counter = _missingEOL_at_EOF_counter;

	bam_task_manager = _bam_task_manager;
	part_queue = _part_queue;
	stats_part_queue = _stats_part_queue;

	pmm_fastq = _pmm_fastq;
	pmm_binary_file_reader = _pmm_binary_file_reader;

	file_type = _file_type;
	kmer_len = _kmer_len;

	// Size and pointer for the buffer
	part_size = 1 << 23;	
	part = nullptr;

	containsNextChromosome = false;

	data_src.SetQueue(binary_pack_queue, pmm_binary_file_reader);
}

//----------------------------------------------------------------------------------
// Destructor - close the files
CFastqReader::~CFastqReader()
{
	if (part)
		pmm_fastq->free(part);
}

//----------------------------------------------------------------------------------
// Set part size of the buffer
bool CFastqReader::SetPartSize(uint64 _part_size)
{	
	if (_part_size < (1 << 20) || _part_size >(1 << 30))
		return false;

	part_size = _part_size;

	return true;
}

void CFastqReader::ProcessBamBinaryPart(uchar* data, uint64 size, uint32 id, uint32 file_no)
{
	pmm_fastq->bam_reserve_gunzip(part, id);

	part_filled = 0;

	uint64_t inpos = 0;
	while (inpos < size)
	{
		inpos += 16;

		uint16_t BGZF_block_size_tmp;
		read_uint16_t(BGZF_block_size_tmp, data, inpos);

		uint64_t BGZF_block_size = BGZF_block_size_tmp + 1ull;  //This integer gives the size of the containing BGZF block minus one.

		uint64_t ISIZE_pos = inpos - 18 + BGZF_block_size - 4;
		uint32_t ISIZE;
		read_uint32_t(ISIZE, data, ISIZE_pos);

		if (part_filled + ISIZE > part_size)
		{
			bam_task_manager->PushGunzippedPart(part, part_filled, id, file_no);

			uchar* prepare_for_splitter_data;
			uint64 prepare_for_splitter_size;
			uint32 prepare_for_splitter_id;
			uint32 prepare_for_splitter_file_no;
			while (bam_task_manager->TakeNextPrepareForSplitterTaskIfExists(prepare_for_splitter_data, prepare_for_splitter_size, prepare_for_splitter_id, prepare_for_splitter_file_no))
			{
				PreparePartForSplitter(prepare_for_splitter_data, prepare_for_splitter_size, prepare_for_splitter_id, prepare_for_splitter_file_no);
				bam_task_manager->NotifySplitterPrepareTaskDone();
			}

			pmm_fastq->bam_reserve_gunzip(part, id);
			part_filled = 0;
		}

		z_stream stream;
		stream.zalloc = Z_NULL;
		stream.zfree = Z_NULL;
		stream.opaque = Z_NULL;
		stream.avail_in = BGZF_block_size - 18;
		stream.next_in = data + inpos;
		stream.avail_out = part_size - part_filled;
		stream.next_out = part + part_filled;
		if (inflateInit2(&stream, -15) != Z_OK)
		{
			std::ostringstream ostr;
			ostr << "inflateInit2 error";
			ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		int ret = inflate(&stream, Z_NO_FLUSH);
		switch (ret)
		{
		case Z_NEED_DICT:
			ret = Z_DATA_ERROR;     /* and fall through */
		case Z_DATA_ERROR:
		{
			std::ostringstream ostr;
			ostr << "Some error (Z_DATA_ERROR) while reading bam file, please contact authors, CODE: FastqReader_" << __LINE__;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		case Z_MEM_ERROR:
		{
			inflateEnd(&stream);

			std::ostringstream ostr;
			ostr << "Some error (Z_MEM_ERROR) while reading bam file, please contact authors, CODE: FastqReader_" << __LINE__;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		}
		if (ret == Z_STREAM_END)
			if (inflateEnd(&stream) != Z_OK) {
				std::ostringstream ostr;
				ostr << "Some error while reading gzip file (inflateEnd) in (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
		else
		{
			std::ostringstream ostr;
			ostr << "Error: Some error in bam file decompression, please contact authors, CODE: FastqReader_" << __LINE__;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		if (stream.avail_in != 8)
		{
			std::ostringstream ostr;
			ostr << "Error: Some error in bam file decompression, please contact authors, CODE: FastqReader_" << __LINE__;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		part_filled += ISIZE;

		inpos += BGZF_block_size_tmp - 18 + 1;

		//TODO: check CRC? in the future I may want to use libdeflate for bam decompression and crc checks

		uint64_t CRC_pos = inpos - 8;
		uint32_t CRC;
		read_uint32_t(CRC, data, CRC_pos);

		uLong crc = crc32(0L, Z_NULL, 0);
		crc = crc32(crc, part + part_filled - ISIZE, ISIZE);
		if (crc != CRC)
		{
			std::ostringstream ostr;
			ostr << "Error: BGZF CRC check error. Input file may be corrupted";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

	}
	bam_task_manager->PushGunzippedPart(part, part_filled, id, file_no);
	part = nullptr;
}

void CFastqReader::PreparePartForSplitter(uchar* data, uint64 size, uint32 /*id*/, uint32 file_no)
{
	auto& state = bam_task_manager->splitter_prepare_state;
	uint64_t bpos = 0;
	//if first in file, skip header
	if (file_no != state.current_file_no)
	{
		state.current_file_no = file_no;

		if (strncmp((const char*)data, "BAM\1", 4) != 0)
		{
			std::ostringstream ostr;
			ostr << "BAM\\1 magic string missed";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		int32_t l_text = 0;
		bpos = 4;
		read_int32_t(l_text, data, bpos);

		bpos += l_text;
		int32_t n_ref = 0;

		read_int32_t(n_ref, data, bpos);
		for (int32_t i = 0; i < n_ref; ++i)
		{
			int32_t l_name = 0;
			read_int32_t(l_name, data, bpos);
			bpos += l_name;
			bpos += 4;
		}

		uint64 skip_header_offset = bpos;
		//if contain header, remove it

		memmove(data, data + skip_header_offset, size - skip_header_offset);
		size -= skip_header_offset;
		skip_header_offset = 0;
		bpos = 0;
	}

	if (state.prev_part_size)
	{
		//assure block size present in prev_part
		if (state.prev_part_size < 4)
		{
			memcpy(state.prev_part_data + state.prev_part_size, data, 4 - state.prev_part_size);
			bpos += 4 - state.prev_part_size;
			state.prev_part_size += 4 - state.prev_part_size;
		}
		int32_t prev_block_size;
		uint64_t prev_bpos = 0;
		read_int32_t(prev_block_size, state.prev_part_data, prev_bpos);


		uint64_t missing = prev_block_size - (state.prev_part_size - 4);
		memcpy(state.prev_part_data + state.prev_part_size, data + bpos, missing);
		state.prev_part_size += missing;
		bpos += missing;

		if (part_queue)
		{
			part_queue->push(state.prev_part_data, state.prev_part_size, ReadType::na);
		}
		else if (stats_part_queue)
		{
			stats_part_queue->push(state.prev_part_data, state.prev_part_size, ReadType::na);
		}
		else
		{
			std::ostringstream ostr;
			ostr << "Error: Should never be here, please contact authors, CODE: FastqReader_" << __LINE__;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		//TODO: this memmove could be avoided, the part queue could store starting pos, apropriate changes will be required in splitters
		//The amount of moved data here may be quite large (close to size)
		//on the other hand bam records should rare span multiple bgzf blocks and in case of newer BAM files it should not happen
		//so for now I will let it be this way
		memmove(data, data + bpos, size - bpos);
		size -= bpos;
		bpos = 0;
		state.prev_part_data = nullptr;
		state.prev_part_size = 0;
	}

	while (true)
	{
		if (bpos == size) // this pack contains full records, may be passed 
		{

			if (part_queue)
			{
				part_queue->push(data, size, ReadType::na);
			}
			else if (stats_part_queue)
			{
				stats_part_queue->push(data, size, ReadType::na);
			}
			else
			{
				std::ostringstream ostr;
				ostr << "Error: Should never be here, please contact authors, CODE: FastqReader_" << __LINE__;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			break;
		}

		bool span_multiple_bgzf_blocks = false;

		int32_t block_size;
		if (bpos + 4 > size) //part of block_size is in the next gunzipped part
		{
			span_multiple_bgzf_blocks = true;
		}
		else
		{
			read_int32_t(block_size, data, bpos);

			if (bpos + block_size > size)
			{
				span_multiple_bgzf_blocks = true;
				bpos -= 4;
			}
		}

		if (span_multiple_bgzf_blocks)
		{
			pmm_fastq->reserve(state.prev_part_data);
			memcpy(state.prev_part_data, data + bpos, size - bpos);
			state.prev_part_size = size - bpos;
			if (part_queue)
			{
				part_queue->push(data, bpos, ReadType::na);
			}
			else if (stats_part_queue)
			{
				stats_part_queue->push(data, bpos, ReadType::na);
			}
			else
			{
				std::ostringstream ostr;
				ostr << "Error: Should never be here, please contact authors, CODE: FastqReader_" << __LINE__;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			break;
		}

		bpos += block_size;
	}
}

//----------------------------------------------------------------------------------
// Read a part of the file in bam file format
void CFastqReader::ProcessBam()
{
	CBamTaskManager::TaskType taskType;
	uchar* data = nullptr;
	uint64 size = 0;
	uint32 id = (uint32)(-1);
	uint32 file_no = (uint32)(-1);

	while (bam_task_manager->PopTask(taskType, data, size, id, file_no))
	{
		if (taskType == CBamTaskManager::TaskType::Gunzip)
		{
			ProcessBamBinaryPart(data, size, id, file_no);
			pmm_binary_file_reader->free(data);
			bam_task_manager->NotifyIdFinished(id);
			pmm_fastq->bam_notify_id_finished(id);
		}
		else if (taskType == CBamTaskManager::TaskType::PrepareForSplitter)
		{
			PreparePartForSplitter(data, size, id, file_no);
			bam_task_manager->NotifySplitterPrepareTaskDone();
		}
		else
		{
			std::ostringstream ostr;
			ostr << "Error: should never be here, plase contact authors, CODE: FastqReader_" << __LINE__;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
	}
}

//----------------------------------------------------------------------------------
// Read a part of the file in multi line fasta format
bool CFastqReader::GetPartFromMultilneFasta(bool allow_unexpected_end_of_gzip_stream, uchar *&_part, uint64 &_size)
{
	uint64 readed = 0;

	if (!containsNextChromosome)
	{
		if (data_src.Finished())
			return false;
	}

	bool last_in_file;
	bool first_in_file;
	readed = data_src.read(part + part_filled, (part_size - 1) - part_filled, allow_unexpected_end_of_gzip_stream, last_in_file, first_in_file); //part_size - 1 to eventually append EOL if not present at EOF

	int64 total_filled = part_filled + readed;

	if (readed && last_in_file) //mkokot_TODO: because all EOLs was removed for prev part, the reading module really needs to be refactored
		FixEOLIfNeeded(part, total_filled);
	
	if(first_in_file && total_filled)
		if (part[0] != '>')
		{
			std::ostringstream ostr;
			ostr << "Wrong input file";
			ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

	int64 last_header_pos = 0;
	int64 pos = 0;
	for (int64 i = 0; i < total_filled; ++i)//find last '>' and remove EOLs
	{
		while (part[i] == '>') //issue 116, (if -> while) USE CASE: when there is empty seqence (only header, no data)
		{
			int64 tmp = i;
			bool next_line = SkipNextEOL(part, i, total_filled);
			if (!next_line)
				i = total_filled; 
			copy(part + tmp, part + i, part + pos);
			last_header_pos = pos;
			pos += i - tmp;
			if (!next_line)
				break;
		}
		if (i < total_filled && part[i] != '\n' && part[i] != '\r')
		{
			part[pos++] = part[i];
		}
	}

	_part = part;
	if (last_header_pos == 0)//data in block belong to one seq ('>' is the first symbol or do not occur at all)
	{
		part_filled = kmer_len - 1;
		_size = pos;
		if (_size < part_filled || last_in_file) //fixes 96 (first part) and prevents copying if it was the last part of a single input file
			part_filled = 0;
		pmm_fastq->reserve(part);
		copy(_part + _size - part_filled, _part + _size, part);
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

FORCE_INLINE void CFastqReader::FixEOLIfNeeded(uchar* part, int64& size)
{
	//Just in case. In fact this method should never be called because if size is 0, the last_in_file flag should not be set
	//But I am not sure if it is true in all cases
	if(size == 0)
		return;
	uchar c = part[size - 1];
	if (c != '\n' && c != '\r')
	{
		missingEOL_at_EOF_counter->RegisterMissingEOL();
		part[size++] = '\n'; //append fake EOL
	}
}

FORCE_INLINE bool CFastqReader::GetNextSymbOfLongReadRecord(bool allow_unexpected_end_of_gzip_stream, uchar& res, int64& p, int64& size)
{
	if (p == size)
	{
		bool last_in_file;
		size = data_src.read(part, part_size - 1, allow_unexpected_end_of_gzip_stream, last_in_file);//part_size - 1 to eventually append EOL if not present at EOF
		if(last_in_file)
			FixEOLIfNeeded(part, size);
		p = 0;
		if (size == 0)
			return false;
	}
	res = part[p++];
	return true;
}

void CFastqReader::CleanUpAfterLongFastaRead(bool allow_unexpected_end_of_gzip_stream)
{
	pmm_fastq->reserve(part);

	uchar symb;
	int64 in_part = 0;
	int64 skip_pos = 0;	
	while (GetNextSymbOfLongReadRecord(allow_unexpected_end_of_gzip_stream, symb, skip_pos, in_part))
	{
		if (symb == '\n' || symb == '\r')
			;
		else
		{				
			if (symb != '>')
			{
				std::ostringstream ostr;
				ostr << "Wrong input file";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			std::copy(part + skip_pos - 1, part + in_part, part);
			part_filled = in_part - (skip_pos - 1);
			return;				
		}				
	}

	//the file has ended
	part_filled = 0;
}
void CFastqReader::CleanUpAfterLongFastqRead(bool allow_unexpected_end_of_gzip_stream, uint32 number_of_lines_to_skip)
{
	pmm_fastq->reserve(part);

	uchar symb;
	int64 in_part = 0;
	int64 skip_pos = 0;
	uint32 state = 0; // 0 - skipping remaining EOLs, 1 - skipping line
	while (GetNextSymbOfLongReadRecord(allow_unexpected_end_of_gzip_stream, symb, skip_pos, in_part))
	{
		if (state == 0)
		{
			if (symb == '\n' || symb == '\r')
				;
			else
			{
				if (number_of_lines_to_skip)
					state = 1;
				else
				{
					if (symb != '@')
					{
						std::ostringstream ostr;
						ostr << "Wrong input file";
						ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
						CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
					}
					std::copy(part + skip_pos - 1, part + in_part, part);
					part_filled = in_part - (skip_pos - 1);
					return;
				}
			}
		}
		else if (state == 1)
		{
			if (symb == '\n' || symb == '\r')
			{
				--number_of_lines_to_skip;
				state = 0;
			}
		}
	}

	//the file has ended
	part_filled = 0;
}

bool CFastqReader::GetPartNew(bool allow_unexpected_end_of_gzip_stream, uchar *&_part, uint64 &_size, ReadType& read_type)
{
	if (file_type == InputType::MULTILINE_FASTA)
	{
		read_type = ReadType::na;
		return GetPartFromMultilneFasta(allow_unexpected_end_of_gzip_stream, _part, _size);
	}

	if (data_src.Finished())
		return false;
	uint64 readed;

	// Read data
	bool last_in_file;
	readed = data_src.read(part + part_filled, (part_size - 1) - part_filled, allow_unexpected_end_of_gzip_stream, last_in_file);//part_size - 1 to eventually append EOL if not present at EOF
	
	//std::cerr.write((char*)part + part_filled, readed);

	int64 total_filled = part_filled + readed;

	if (last_in_file)
		FixEOLIfNeeded(part, total_filled);

	int64 i;

	if (data_src.Finished() && !long_read_in_progress)
	{
		read_type = ReadType::normal_read;
		_part = part;
		_size = total_filled;

		part = nullptr;
		return true;
	}

	if (!total_filled) //better fix for issue 95
	{
		_part = part;
		_size = 0;
	}									// Look for the end of the last complete record in a buffer
	else if (file_type == InputType::FASTA || file_type == InputType::KMC)			// FASTA files
	{
		if (long_read_in_progress)
		{			
			//check if there is EOL in the data
			int64 pos = 0;
			for (; pos < total_filled; ++pos)
			{
				if (part[pos] == '\n' || part[pos] == '\r')
				{
					long_read_in_progress = false;
					break;
				}
			}

			if (!long_read_in_progress)
			{
				_part = part;
				_size = pos;

				//all from this part was readed, maybe there is another EOL character in the file
				if(pos == total_filled)
					CleanUpAfterLongFastaRead(allow_unexpected_end_of_gzip_stream);
				else //there is still some important data in the part!!
				{
					//skip possible eol
					for (; pos < total_filled; ++pos)
						if (part[pos] != '\n' && part[pos] != '\r')
							break;

					pmm_fastq->reserve(part);
					std::copy(_part + pos, _part + total_filled, part);
					part_filled = total_filled - pos;
				}
			}
			else
			{
				_part = part;
				_size = total_filled;
				pmm_fastq->reserve(part);
				std::copy(_part + total_filled - kmer_len + 1, _part + total_filled, part);
				part_filled = kmer_len - 1;
			}
			return true;
		}
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
			if (readed_lines == 2) //because if successfully readed full 2 lines in the worst case there is only one read that begins at the buffer start, and there is only onle fasta record
			{
				std::ostringstream ostr;
				ostr << "some error while reading fasta file, please contact authors";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			k = 4 - readed_lines;
			if (line_start[k] != 0)
			{
				std::ostringstream ostr;
				ostr << "some error while reading fasta file, please contact authors";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			if (part[0] != '>')
			{
				std::ostringstream ostr;
				ostr << "Wrong input file";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			if (last_in_file)
			{
				part_filled = 0;
				_size = 0;
				return true;
			}
			if (readed_lines == 1)
			{
				long_read_in_progress = true;				

				_part = part;
				_size = total_filled;
				read_type = ReadType::long_read;

				//copy last k-1 symbols
				pmm_fastq->reserve(part);
				copy(_part + total_filled - kmer_len + 1, _part + total_filled, part);
				part_filled = kmer_len - 1;

				return true;
			}
			else 
			{
				std::ostringstream ostr;
				ostr << "some error while reading fasta file, please contact authors";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}

			return true;
		}

		_part = part;
		_size = line_end[k + 1];
		read_type = ReadType::normal_read;
	}
	else
	{
		if (long_read_in_progress)
		{
			//check if there is EOL in the data
			int64 pos = 0;
			for (; pos < total_filled; ++pos)
			{
				if (part[pos] == '\n' || part[pos] == '\r')
				{
					long_read_in_progress = false;
					break;
				}
			}

			if (!long_read_in_progress)
			{
				_part = part;
				_size = pos;

				//count the number of lines to skip after the part that is currently in the buffer (after read)
				uint32 no_lines = 2;//qual header and qual
				uint32 state = 0; //0 - skipping remaining EOLs, 1 - skipping line content
				for (; pos < total_filled; ++pos)
				{
					if (state == 0)
					{
						if (part[pos] == '\n' || part[pos] == '\r')
							state = 1;
					}
					else if (state == 1)
					{
						if (part[pos] == '\n' || part[pos] == '\r')
						{
							--no_lines;
							state = 0;
						}
					}
				}
				CleanUpAfterLongFastqRead(allow_unexpected_end_of_gzip_stream, no_lines);
			}
			else
			{
				_part = part;
				_size = total_filled;
				pmm_fastq->reserve(part);
				std::copy(_part + total_filled - kmer_len + 1, _part + total_filled, part);
				part_filled = kmer_len - 1;
			}
			return true;
		}

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

		if (!success) //wrong input file or very long read 
		{
			if (readed_lines == 4) //because if successfully reade full 4 lines in the worst case there is only one read that begins at the buffer start, and there is only onle fastq record
			{
				std::ostringstream ostr;
				ostr << "some error while reading fastq file, please contact authors";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			k = 8 - readed_lines;
			if (line_start[k] != 0)
			{
				std::ostringstream ostr;
				ostr << "some error while reading fastq file, please contact authors";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			if (part[0] != '@')
			{
				std::ostringstream ostr;
				ostr << "Wrong input file";
				ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			if (last_in_file)
			{
				part_filled = 0;
				_size = 0;
				return true;
			}
			if (readed_lines == 1)
			{
				long_read_in_progress = true;

				_part = part;
				_size = total_filled;
				read_type = ReadType::long_read;

				//copy last k-1 symbols
				pmm_fastq->reserve(part);
				copy(_part + total_filled - kmer_len + 1, _part + total_filled, part);
				part_filled = kmer_len - 1;

				return true;
			}
			else //whole read had beed read
			{
				//readed_lines == 1 -> read header, and part of read
				//readed_lines == 2 -> read header + read, and (possibly empty) part of qual header
				//readed_lines == 3 -> read header + read + qual header, and (possibly empty) part of qual 
				//readed_lines == 4 -> should not be possible here
				//so here only read_lines == 2 or 3 is possible

				long_read_in_progress = false;
				_part = part;
				_size = line_end[k + 1];
				read_type = ReadType::long_read;

				//Some clean up to assure in next function call the buffer starts with @ or nothing (file finished)

				uint32 number_of_lines_to_skip;
				if (readed_lines == 2)
					number_of_lines_to_skip = 2; // means inside qual header
				else if (readed_lines == 3)
					number_of_lines_to_skip = 1; //means inside qual
				else //shold not be possible
				{
					std::ostringstream ostr;
					ostr << "some error while reading fastq file, please contact authors";
					ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
					CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
				}

				CleanUpAfterLongFastqRead(allow_unexpected_end_of_gzip_stream, number_of_lines_to_skip);
				return true;

			}

		}
		_part = part;
		_size = line_end[k + 3];
		read_type = ReadType::normal_read;
	}
	// Allocate new memory for the buffer

	pmm_fastq->reserve(part);
	copy(_part + _size, _part + total_filled, part);
	part_filled = total_filled - _size;

	return true;
}

//----------------------------------------------------------------------------------
// Skip to next EOL from the current position in a buffer
bool CFastqReader::SkipNextEOL(uchar *part, int64 &pos, int64 size)
{
	int64 i;
	for (i = pos; i < size - 1; ++i)
		if ((part[i] == '\n' || part[i] == '\r') && !(part[i + 1] == '\n' || part[i + 1] == '\r'))
			break;

	if (i >= size - 1)
		return false;

	pos = i + 1;

	return true;
}

void CFastqReader::GetFullLineFromEnd(int64& line_sart, int64& line_end, uchar* buff, int64& pos)
{
	while (pos >= 0 && buff[pos] != '\n' && buff[pos] != '\r')
		--pos;
	line_end = pos + 1;
	if (pos >= 0 && (buff[pos] == '\n' || buff[pos] == '\r'))
	{
		--pos;
		if (pos >= 0 && buff[pos] != buff[pos + 1] && (buff[pos] == '\n' || buff[pos] == '\r'))
			--pos;
	}	
	while (pos >= 0 && buff[pos] != '\n' && buff[pos] != '\r')
		--pos;
	line_sart = pos + 1;
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
			std::ostringstream ostr;
			ostr << "Error while reading gz file";
			ostr << " (" << __FILE__ << ": " << __LINE__ << ")";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		stream.avail_in = (uint32)in_data_size;
		stream.next_in = in_data;
		break;
	default:
		break;
	}
}

//----------------------------------------------------------------------------------
//returns false if there is
// a) end of a single file (or if binary reader interrupted reading (in stats mode)), or
// b) there is nothing more and will never be anything more
bool CFastqReaderDataSrc::pop_pack(uchar*& data, uint64& size, FilePart& file_part, CompressionType& mode, bool& last_in_file)
{
	end_reached = !binary_pack_queue->pop(data, size, file_part, mode);
	if (file_part == FilePart::End)
	{
		last_in_file = true;
		in_progress = false;
		return false;
	}
	return !end_reached;
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
uint64 CFastqReaderDataSrc::read(uchar* buff, uint64 size, bool allow_unexpected_end_of_gzip_stream, bool& last_in_file, bool& first_in_file)
{
	last_in_file = false;
	first_in_file = false;
	if (!in_progress)
	{
		if (!pop_pack(in_data, in_data_size, file_part, compression_type, last_in_file))
			return 0;
		in_progress = true;
		first_in_file = true;
		init_stream();
	}

	//reading
	if (compression_type == CompressionType::gzip)
	{
		stream.next_out = buff;
		stream.avail_out = (uint32)size;
		int ret = Z_OK;
		do
		{
			if (!stream.avail_in)
			{
				pmm_binary_file_reader->free(in_data);
				in_data = nullptr;
				if (!pop_pack(in_data, in_data_size, file_part, compression_type, last_in_file)) {
					auto ret_val = size - stream.avail_out;

					if (inflateEnd(&stream) != Z_OK) {
						std::ostringstream ostr;
						ostr << "Some error while reading gzip file (inflateEnd) in (" << __FILE__ << ": " << __LINE__ << ")";
						CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
					}
					assert(last_in_file);

					//there is no more data but we don't have Z_STREAM_END
					if (!allow_unexpected_end_of_gzip_stream && ret != Z_STREAM_END) {
						std::ostringstream ostr;
						ostr << "Unexpected end of gzip file";
						CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
					}
					return ret_val; //mkokot_TODO: is this correct now?
				}
				stream.avail_in = (uint32)in_data_size;
				stream.next_in = in_data;
			}

			ret = inflate(&stream, Z_NO_FLUSH);

			switch (ret)
			{
			case Z_NEED_DICT:
				ret = Z_DATA_ERROR;     /* and fall through */
			case Z_DATA_ERROR:
			{
				std::ostringstream ostr;
				ostr << "Some error while reading gzip file in (" << __FILE__ << ": " << __LINE__ << ")";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			case Z_MEM_ERROR:
			{
				inflateEnd(&stream);
				return ret;
			}
			}

			if (ret == Z_STREAM_END)
			{
				uchar* tmp_data = nullptr;
				uint64 tmp_size = 0;
				bool multistream = stream.avail_in;

				if (stream.avail_in < 2) {
					bool peek_res = binary_pack_queue->peek_next_pack(tmp_data, tmp_size);
					multistream = stream.avail_in || peek_res;
				}
				bool garbage = false;
				if (multistream)
				{
					if (stream.avail_in + tmp_size < 2)
					{
						std::ostringstream ostr;
						ostr << "Some error while reading gzip file in (" << __FILE__ << ": " << __LINE__ << ")";
						CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
					}
					uchar b1, b2;
					if (stream.avail_in >= 2)
					{
						b1 = stream.next_in[0];
						b2 = stream.next_in[1];
					}
					else if (stream.avail_in == 1)
					{
						b1 = stream.next_in[0];
						b2 = tmp_data[0];
					}
					else
					{
						b1 = tmp_data[0];
						b2 = tmp_data[1];
					}
					garbage = b1 != 0x1f || b2 != 0x8b;
				}

				//bool multistream = stream.avail_in || !binary_pack_queue->is_next_last();				
				if (!multistream || garbage)
				{
					pmm_binary_file_reader->free(in_data);
					in_data = nullptr;
					if (inflateEnd(&stream) != Z_OK) {
						std::ostringstream ostr;
						ostr << "Some error while reading gzip file (inflateEnd) in (" << __FILE__ << ": " << __LINE__ << ")";
						CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
					}
					in_progress = false;
					//pull end
					bool queue_end = !pop_pack(in_data, in_data_size, file_part, compression_type, last_in_file);
					if (!queue_end && file_part != FilePart::End && !garbage)
					{
						std::ostringstream ostr;
						ostr << "An internal error occurred. Please contact authors. File: " << __FILE__ << ", line: " << __LINE__;
						CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
					}
					if (garbage)
					{
						//ignore rest (garbage) data of a current file
						FilePart tmpfilepart = file_part;
						while (tmpfilepart != FilePart::End)
						{
							uchar* tmp;
							uint64 tmpsize;
							CompressionType tmpcomptype;
							pop_pack(tmp, tmpsize, tmpfilepart, tmpcomptype, last_in_file);
						}
					}
					last_in_file = true;
					break;
				}
				else //multiple streams in one file
				{
					if (inflateReset(&stream) != Z_OK)
					{
						std::ostringstream ostr;
						ostr << "Error while reading gzip file";
						CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
					}
				}
			}
		} while (stream.avail_out);
		return size - stream.avail_out;
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
				//may be false even if file_part != FilePart::End in stats mode
				auto pop_res = pop_pack(in_data, in_data_size, file_part, compression_type, last_in_file);
				if (!pop_res || file_part == FilePart::End)
				{
					in_progress = false;
					last_in_file = true;
					break;
				}
				in_data_pos = 0;
			}
			uint64 in_left = in_data_size - in_data_pos;
			uint64 out_left = size - out_pos;
			uint64 n_copy = min(in_left, out_left);
			memcpy(buff + out_pos, in_data + in_data_pos, n_copy);
			in_data_pos += n_copy;
			out_pos += n_copy;
		} while (out_pos < size);
		return out_pos;
	}
	else
	{
		std::ostringstream ostr;
		ostr << "unknown compression";
		CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
	}
	assert(false); //should never be here
}
//----------------------------------------------------------------------------------
uint64 CFastqReaderDataSrc::read(uchar* buff, uint64 size, bool allow_unexpected_end_of_gzip_stream, bool& last_in_file)
{
	bool ignore;
	return read(buff, size, allow_unexpected_end_of_gzip_stream, last_in_file, ignore);
}




//************************************************************************************************************
// CWFastqReader - wrapper for multithreading purposes
//************************************************************************************************************
CWFastqReader::CWFastqReader(CKMCParams &Params, CKMCQueues &Queues, CBinaryPackQueue* _binary_pack_queue)
{
	pmm_fastq = Queues.pmm_fastq.get();
	pmm_binary_file_reader = Queues.pmm_binary_file_reader.get();
	//input_files_queue = Queues.input_files_queue;
	binary_pack_queue = _binary_pack_queue;
	missingEOL_at_EOF_counter = Queues.missingEOL_at_EOF_counter.get();
	bam_task_manager = Queues.bam_task_manager.get();
	part_size = Params.fastq_buffer_size; 
	part_queue = Queues.part_queue.get();
	file_type = Params.file_type;
	kmer_len = Params.kmer_len;
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

	CFastqReader fqr(pmm_fastq, file_type, kmer_len, binary_pack_queue, pmm_binary_file_reader, bam_task_manager, part_queue, nullptr, missingEOL_at_EOF_counter);
	fqr.SetPartSize(part_size);
	if (file_type == InputType::BAM)
	{
		fqr.ProcessBam();
	}
	else
	{
		fqr.Init();
		ReadType read_type;
		while (fqr.GetPartNew(false, part, part_filled, read_type))
			part_queue->push(part, part_filled, read_type);
	}
	part_queue->mark_completed();
}



//************************************************************************************************************
// CWStatsFastqReader - wrapper for multithreading purposes
//************************************************************************************************************
CWStatsFastqReader::CWStatsFastqReader(CKMCParams &Params, CKMCQueues &Queues, CBinaryPackQueue* _binary_pack_queue)
{
	pmm_fastq = Queues.pmm_fastq.get();
	pmm_binary_file_reader = Queues.pmm_binary_file_reader.get();

	binary_pack_queue = _binary_pack_queue;
	bam_task_manager = Queues.bam_task_manager.get();

	part_size = Params.fastq_buffer_size;
	stats_part_queue = Queues.stats_part_queue.get();
	file_type = Params.file_type;
	kmer_len = Params.kmer_len;

	missingEOL_at_EOF_counter = Queues.missingEOL_at_EOF_counter.get();
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


	CFastqReader fqr(pmm_fastq, file_type, kmer_len, binary_pack_queue, pmm_binary_file_reader, bam_task_manager, nullptr, stats_part_queue, missingEOL_at_EOF_counter);
	fqr.SetPartSize(part_size);
	if (file_type == InputType::BAM)
	{
		fqr.ProcessBam();
	}
	else
	{
		fqr.Init();
		ReadType read_type;
		while (fqr.GetPartNew(true, part, part_filled, read_type))
		{
			if (part_filled) //if part_filled == 0 we are probably in stats computing mode
				stats_part_queue->push(part, part_filled, read_type);
			else
				pmm_fastq->free(part);
		}
	}
	stats_part_queue->mark_completed();
}


// ***** EOF
