#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.1.0
  Date   : 2018-05-10
*/

#include <algorithm>
#include "defs.h"
#include "fastq_reader.h"
#include "bam_utils.h"
//************************************************************************************************************
// CFastqReader	- reader class
//************************************************************************************************************

uint64 CFastqReader::OVERHEAD_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
// Constructor of FASTA/FASTQ reader
// Parameters:
//    * _mm - pointer to memory monitor (to check the memory limits)
CFastqReader::CFastqReader(CMemoryMonitor *_mm, CMemoryPoolWithBamSupport *_pmm_fastq, input_type _file_type, int _kmer_len, CBinaryPackQueue* _binary_pack_queue, CMemoryPool* _pmm_binary_file_reader, CBamTaskManager* _bam_task_manager, CPartQueue* _part_queue, CStatsPartQueue* _stats_part_queue)
{
	binary_pack_queue = _binary_pack_queue;

	bam_task_manager = _bam_task_manager;
	part_queue = _part_queue;
	stats_part_queue = _stats_part_queue;

	mm = _mm;
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
			if (!bam_task_manager->PushGunzippedPart(part, part_filled, id, file_no))
			{
				pmm_fastq->free(part);
				part = nullptr;
				return;
			}

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
			cerr << "Error: inflateInit2\n";
			exit(1);
		}

		int ret = inflate(&stream, Z_NO_FLUSH);
		switch (ret)
		{
		case Z_NEED_DICT:
			ret = Z_DATA_ERROR;     /* and fall through */
		case Z_DATA_ERROR:
			cerr << "Some error (Z_DATA_ERROR) while reading bam file, please contact authors, CODE: FastqReader_" << __LINE__ << "\n";
			exit(1);
		case Z_MEM_ERROR:
			inflateEnd(&stream);
			exit(ret);
		}
		if (ret == Z_STREAM_END)
			inflateEnd(&stream);
		else
		{
			cerr << "Error: Some error in bam file decompression, please contact authors, CODE: FastqReader_" << __LINE__ << "\n";
			exit(1);
		}
		if (stream.avail_in != 8)
		{
			//cerr << "Err avail_in != 8:" << stream.avail_in << "\n";
			cerr << "Error: Some error in bam file decompression, please contact authors, CODE: FastqReader_" << __LINE__ << "\n";
			exit(1);
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
			cerr << "Error: BGZF CRC check error. Input file may be corrupted\n";
			exit(1);
		}

	}
	if (!bam_task_manager->PushGunzippedPart(part, part_filled, id, file_no))
		pmm_fastq->free(part);
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
			cerr << "BAM\\1 magic string missed\n";
			exit(1);
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
			part_queue->push(state.prev_part_data, state.prev_part_size);
		}
		else if (stats_part_queue)
		{
			if (!stats_part_queue->push(state.prev_part_data, state.prev_part_size))
			{
				pmm_fastq->free(state.prev_part_data);
				pmm_fastq->free(data);
				state.prev_part_data = nullptr;
				state.prev_part_size = 0;
				bam_task_manager->IgnoreRest(pmm_fastq, pmm_binary_file_reader);
				return;
			}
		}
		else
		{
			cerr << "Error: Should never be here, please contact authors, CODE: FastqReader_" << __LINE__ << "\n";
			exit(1);
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
				part_queue->push(data, size);
			}
			else if (stats_part_queue)
			{
				if (!stats_part_queue->push(data, size))
				{
					pmm_fastq->free(data);
					bam_task_manager->IgnoreRest(pmm_fastq, pmm_binary_file_reader);
					return;
				}
			}
			else
			{
				cerr << "Error: Should never be here, please contact authors, CODE: FastqReader_" << __LINE__ << "\n";
				exit(1);
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
				part_queue->push(data, bpos);
			}
			else if (stats_part_queue)
			{
				if (!stats_part_queue->push(data, bpos))
				{
					pmm_fastq->free(data);
					pmm_fastq->free(state.prev_part_data);
					state.prev_part_data = nullptr;
					state.prev_part_size = 0;
					bam_task_manager->IgnoreRest(pmm_fastq, pmm_binary_file_reader);
					return;
				}
			}
			else
			{
				cerr << "Error: Should never be here, please contact authors, CODE: FastqReader_" << __LINE__ << "\n";
				exit(1);
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
			cerr << "Error: should never be here, plase contact authors, CODE: FastqReader_" << __LINE__ << "\n";
		}
	}
}

//----------------------------------------------------------------------------------
// Read a part of the file in multi line fasta format
bool CFastqReader::GetPartFromMultilneFasta(uchar *&_part, uint64 &_size)
{
	uint64 readed = 0;

	if (!containsNextChromosome)
	{
		if (data_src.Finished())
			return false;
	}

	readed = data_src.read(part + part_filled, part_size - part_filled);

	int64 total_filled = part_filled + readed;
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
		if (part[i] != '\n' && part[i] != '\r')
		{
			part[pos++] = part[i];
		}
	}

	_part = part;
	if (last_header_pos == 0)//data in block belong to one seq ('>' is the first symbol or do not occur at all)
	{
		part_filled = kmer_len - 1;
		_size = pos;
		if (_size < part_filled) //fixes 96
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

bool CFastqReader::GetPartNew(uchar *&_part, uint64 &_size)
{
	if (file_type == multiline_fasta)
		return GetPartFromMultilneFasta(_part, _size);

	if (data_src.Finished())
		return false;
	uint64 readed;

	// Read data
	readed = data_src.read(part + part_filled, part_size - part_filled);

	//std::cerr.write((char*)part + part_filled, readed);

	int64 total_filled = part_filled + readed;
	int64 i;

	if (data_src.Finished())
	{
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
	else if (file_type == fasta)			// FASTA files
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
	copy(_part + _size, _part + total_filled, part);
	part_filled = total_filled - _size;

	return true;
}

//----------------------------------------------------------------------------------
// Skip to next EOL from the current position in a buffer
bool CFastqReader::SkipNextEOL(uchar *part, int64 &pos, int64 max_pos)
{
	int64 i;
	for (i = pos; i < max_pos - 2; ++i)
		if ((part[i] == '\n' || part[i] == '\r') && !(part[i + 1] == '\n' || part[i + 1] == '\r'))
			break;

	if (i >= max_pos - 2)
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
			cerr << "Error while reading gz file\n";
			exit(1);
		}
		stream.avail_in = (uint32)in_data_size;
		stream.next_in = in_data;
		break;
	case CompressionType::bzip2:
		_bz_stram.bzalloc = nullptr;
		_bz_stram.bzfree = nullptr;
		_bz_stram.opaque = nullptr;
		_bz_stram.avail_in = 0;
		_bz_stram.next_in = nullptr;
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
				if (!binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type))
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
				uchar* tmp_data = nullptr;
				uint64 tmp_size = 0;
				bool multistream = stream.avail_in || binary_pack_queue->peek_next_pack(tmp_data, tmp_size);
				bool garbage = false;
				if (multistream)
				{
					if (stream.avail_in + tmp_size < 2)
					{
						cerr << "Some error while reading gzip file\n";
						exit(1);
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
					inflateEnd(&stream);
					in_progress = false;
					//pull end
					bool queue_end = !binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type);
					if (!queue_end && file_part != FilePart::End && !garbage)
					{
						cerr << "Error: An internal error occurred. Please contact authors\n";
					}
					if (garbage)
						binary_pack_queue->ignore_rest();
					break;
				}
				else //multiple streams in one file
				{
					//equivalent of inflateReset 
					//inflateEnd(&stream);
					//if (inflateInit2(&stream, 31) != Z_OK)
					//{
					//	cerr << "Error while reading gzip file\n";
					//	exit(1);
					//}
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
			memcpy(buff + out_pos, in_data + in_data_pos, n_copy);
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
	bam_task_manager = Queues.bam_task_manager;
	part_size = Params.fastq_buffer_size;
	part_queue = Queues.part_queue;
	file_type = Params.file_type;
	kmer_len = Params.p_k;


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

	fqr = new CFastqReader(mm, pmm_fastq, file_type, kmer_len, binary_pack_queue, pmm_binary_file_reader, bam_task_manager, part_queue, nullptr);
	fqr->SetPartSize(part_size);
	if (file_type == bam)
	{
		fqr->ProcessBam();
	}
	else
	{
		fqr->Init();
		while (fqr->GetPartNew(part, part_filled))
			part_queue->push(part, part_filled);
	}
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
	bam_task_manager = Queues.bam_task_manager;

	part_size = Params.fastq_buffer_size;
	stats_part_queue = Queues.stats_part_queue;
	file_type = Params.file_type;
	kmer_len = Params.p_k;

	fqr = nullptr;
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

	fqr = new CFastqReader(mm, pmm_fastq, file_type, kmer_len, binary_pack_queue, pmm_binary_file_reader, bam_task_manager, nullptr, stats_part_queue);
	fqr->SetPartSize(part_size);
	if (file_type == bam)
	{
		fqr->ProcessBam();
	}
	else
	{
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
	}
	delete fqr;
	stats_part_queue->mark_completed();
}


// ***** EOF
