/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _FASTQ_READER_H
#define _FASTQ_READER_H

#include "defs.h"
#include "params.h"
#include <stdio.h>
#include <iostream>

#include "libs/zlib.h"
#include "libs/bzlib.h"

using namespace std;

//************************************************************************************************************
// data source for FASTA/FASTQ reader
//************************************************************************************************************
class CFastqReaderDataSrc
{
	z_stream stream;	
	bz_stream _bz_stram;
	uchar* in_buffer;
	CBinaryPackQueue* binary_pack_queue;
	CMemoryPool *pmm_binary_file_reader;
	bool in_progress = false;
	bool end_reached = false;
	FilePart file_part;
	CompressionType compression_type;
	uchar* in_data;
	uint64 in_data_size;
	uint64 in_data_pos; //for plain
	void init_stream();
public:
	inline void SetQueue(CBinaryPackQueue* _binary_pack_queue, CMemoryPool *_pmm_binary_file_reader);
	inline bool Finished();
	uint64 read(uchar* buff, uint64 size, bool& last_in_file);
	void IgnoreRest()
	{
		if (in_data)
			pmm_binary_file_reader->free(in_data);
		in_data = nullptr;
		//clean queue
		while (binary_pack_queue->pop(in_data, in_data_size, file_part, compression_type))
		{
			if(in_data_size)
				pmm_binary_file_reader->free(in_data);
			in_data = nullptr;
		}
		switch (compression_type)
		{
		case CompressionType::plain:
			break;
		case CompressionType::gzip:
			inflateEnd(&stream);
			break;
		case CompressionType::bzip2:
			BZ2_bzDecompressEnd(&_bz_stram);
			break;
		default:
			break;
		}

	}
};


//************************************************************************************************************
// FASTA/FASTQ reader class
//************************************************************************************************************
class CFastqReader {	
	CBinaryPackQueue* binary_pack_queue;

	CBamTaskManager* bam_task_manager = nullptr; //only for bam input

	CMemoryMonitor *mm;	
	CMemoryPoolWithBamSupport* pmm_fastq;
	CMissingEOL_at_EOF_counter* missingEOL_at_EOF_counter;

	CMemoryPool *pmm_binary_file_reader;
	CPartQueue *part_queue;
	CStatsPartQueue *stats_part_queue;

	string input_file_name;
	input_type file_type;
	int kmer_len;
	
	CFastqReaderDataSrc data_src;

	uint64 part_size;
	
	uchar *part;
	uint64 part_filled;

	bool long_read_in_progress = false;
	
	bool containsNextChromosome; //for multiline_fasta processing

	bool SkipNextEOL(uchar *part, int64 &pos, int64 size);

	void GetFullLineFromEnd(int64& line_sart, int64& line_end, uchar* buff, int64& pos);
	
	void ProcessBamBinaryPart(uchar* data, uint64 size, uint32 id, uint32 file_no);
	void PreparePartForSplitter(uchar* data, uint64 size, uint32 id, uint32 file_no);

	bool GetNextSymbOfLongReadRecord(uchar& res, int64& p, int64& size);

	void CleanUpAfterLongFastqRead(uint32 number_of_lines_to_skip);

	void CleanUpAfterLongFastaRead();
	void FixEOLIfNeeded(uchar* part, int64& size);	
public:
	CFastqReader(CMemoryMonitor *_mm, CMemoryPoolWithBamSupport *_pmm_fastq, input_type _file_type, int _kmer_len, 
		CBinaryPackQueue* _binary_pack_queue, CMemoryPool* _pmm_binary_file_reader, CBamTaskManager* _bam_task_manager, 
		CPartQueue* _part_queue, CStatsPartQueue* _stats_part_queue, CMissingEOL_at_EOF_counter* _missingEOL_at_EOF_counter);
	~CFastqReader();

	static uint64 OVERHEAD_SIZE;

	bool SetNames(string _input_file_name);
	bool SetPartSize(uint64 _part_size);
	bool OpenFiles();

	bool GetPartFromMultilneFasta(uchar *&_part, uint64 &_size);
	
	void ProcessBam();	

	bool GetPart(uchar *&_part, uint64 &_size);

	bool GetPartNew(uchar *&_part, uint64 &_size, ReadType& read_type);
	void Init()
	{
		pmm_fastq->reserve(part);
		part_filled = 0;
	}

	void IgnoreRest()
	{
		data_src.IgnoreRest();
	}

};

//************************************************************************************************************
// Wrapper for FASTA/FASTQ reader class - for multithreading purposes
//************************************************************************************************************
class CWFastqReader {
	CMemoryMonitor *mm;
	CMemoryPoolWithBamSupport *pmm_fastq;
	CMemoryPool *pmm_binary_file_reader;

	CFastqReader *fqr;
	uint64 part_size;
	CBinaryPackQueue* binary_pack_queue;
	CBamTaskManager* bam_task_manager = nullptr; //only for bam input
	CPartQueue *part_queue;
	CStatsPartQueue *stats_part_queue;

	input_type file_type;	
	int kmer_len;

	CMissingEOL_at_EOF_counter* missingEOL_at_EOF_counter;

public:
	CWFastqReader(CKMCParams &Params, CKMCQueues &Queues, CBinaryPackQueue* _binary_pack_queue);
	~CWFastqReader();

	void operator()();
};



//************************************************************************************************************
// Wrapper for FASTA/FASTQ reader class (stats mode) - for multithreading purposes
//************************************************************************************************************
class CWStatsFastqReader {
	CMemoryMonitor *mm;
	CMemoryPoolWithBamSupport *pmm_fastq;
	CMemoryPool *pmm_binary_file_reader;
	CFastqReader *fqr;
	uint64 part_size;
	CBamTaskManager* bam_task_manager = nullptr; //only for bam input
	CStatsPartQueue *stats_part_queue;
	input_type file_type;	
	int kmer_len;
	CBinaryPackQueue* binary_pack_queue;
	CMissingEOL_at_EOF_counter* missingEOL_at_EOF_counter;
public:
	CWStatsFastqReader(CKMCParams &Params, CKMCQueues &Queues, CBinaryPackQueue* _binary_pack_queue);
	~CWStatsFastqReader();

	void operator()();
};
#endif

// ***** EOF
