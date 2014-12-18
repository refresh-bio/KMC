/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.1
  Date   : 2014-12-18
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
// FASTA/FASTQ reader class
//************************************************************************************************************
class CFastqReader {
	typedef enum {m_plain, m_gzip, m_bzip2} t_mode;

	CMemoryMonitor *mm;
	CMemoryPool *pmm_fastq;

	string input_file_name;
	input_type file_type;
	int kmer_len;
	t_mode mode;

	FILE *in;
	gzFile_s *in_gzip;
	BZFILE *in_bzip2;
	int bzerror;

	uint64 part_size;
	
	uchar *part;
	uint64 part_filled;
	
	uint32 gzip_buffer_size;
	uint32 bzip2_buffer_size;

	bool containsNextChromosome; //for multiline_fasta processing

	bool SkipNextEOL(uchar *part, int64 &pos, int64 max_pos);

	bool IsEof();

public:
	CFastqReader(CMemoryMonitor *_mm, CMemoryPool *_pmm_fastq, input_type _file_type, uint32 _gzip_buffer_size, uint32 _bzip2_buffer_size, int _kmer_len);
	~CFastqReader();

	static uint64 OVERHEAD_SIZE;

	bool SetNames(string _input_file_name);
	bool SetPartSize(uint64 _part_size);
	bool OpenFiles();

	bool GetPartFromMultilneFasta(uchar *&_part, uint64 &_size);
	bool GetPart(uchar *&_part, uint64 &_size);
};

//************************************************************************************************************
// Wrapper for FASTA/FASTQ reader class - for multithreading purposes
//************************************************************************************************************
class CWFastqReader {
	CMemoryMonitor *mm;
	CMemoryPool *pmm_fastq;

	CFastqReader *fqr;
	string file_name;
	uint64 part_size;
	CInputFilesQueue *input_files_queue;
	CPartQueue *part_queue;
	input_type file_type;
	uint32 gzip_buffer_size;
	uint32 bzip2_buffer_size;
	int kmer_len;

public:
	CWFastqReader(CKMCParams &Params, CKMCQueues &Queues);
	~CWFastqReader();

	void operator()();
};



//************************************************************************************************************
// Wrapper for FASTA/FASTQ reader class (stats mode) - for multithreading purposes
//************************************************************************************************************
class CWStatsFastqReader {
	CMemoryMonitor *mm;
	CMemoryPool *pmm_fastq;

	CFastqReader *fqr;
	string file_name;
	uint64 part_size;
	CInputFilesQueue *input_files_queue;
	CStatsPartQueue *stats_part_queue;
	input_type file_type;
	uint32 gzip_buffer_size;
	uint32 bzip2_buffer_size;
	int kmer_len;

public:
	CWStatsFastqReader(CKMCParams &Params, CKMCQueues &Queues);
	~CWStatsFastqReader();

	void operator()();
};
#endif

// ***** EOF
