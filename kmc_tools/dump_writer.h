/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _DUMP_WRITER_H
#define _DUMP_WRITER_H
#include "defs.h"
#include "kmer.h"
#include "nc_utils.h"
#include "config.h"
#include <fstream>


//wrapper to simplify interface

//For kmc1 input and kmc2 input without -s parameter
template<typename KMCDB, unsigned SIZE, bool SORTED>
class CKMCDBForDump
{
	KMCDB kmcdb;
public:
	CKMCDBForDump() :
		kmcdb(CConfig::GetInstance().headers.front(), CConfig::GetInstance().input_desc.front(), CConfig::GetInstance().percent_progress, KMCDBOpenMode::sequential){}
	bool NextKmer(CKmer<SIZE>& kmer, uint32& counter)
	{
		return kmcdb.NextKmerSequential(kmer, counter);
	}
};


//specialization for -s parameter nad kmc2 input
template<unsigned SIZE>
class CKMCDBForDump<CKMC2DbReader<SIZE>, SIZE, true>
{
	CKMC2DbReader<SIZE>* kmcdb;
	CBundle<SIZE> bundle;
public:
	CKMCDBForDump() :
		kmcdb(new CKMC2DbReader<SIZE>(CConfig::GetInstance().headers.front(), CConfig::GetInstance().input_desc.front(), CConfig::GetInstance().percent_progress, KMCDBOpenMode::sorted)), 
		bundle(kmcdb){}
	bool NextKmer(CKmer<SIZE>& kmer, uint32& counter)
	{
		if(!bundle.Finished())
		{
			kmer = bundle.TopKmer();
			counter = bundle.TopCounter();
			bundle.Pop();
			return true;
		}
		return false;
	}
};

template<typename KMCDB, unsigned SIZE>
class CDumpWriter
{
	static const uint32 OVERHEAD_SIZE = 1000;
	KMCDB& kmcdb;
	COutputDesc output_desc;
	uint32 kmer_len;
	uint32 kmer_bytes;
	CConfig& config;
	uint32 in_first_byte;
	char* buf;
	uint32 buf_size;
	uint32 buf_pos;
	struct DumpOpt
	{
		char* opt_ACGT;
		DumpOpt()
		{
			opt_ACGT = new char[1024];
			char codes[] = { 'A', 'C', 'G', 'T' };
			uint32 pos = 0;
			for (uint32 kmer = 0; kmer < 256; ++kmer)
			{
				opt_ACGT[pos++] = codes[(kmer >> 6) & 3];
				opt_ACGT[pos++] = codes[(kmer >> 4) & 3];
				opt_ACGT[pos++] = codes[(kmer >> 2) & 3];
				opt_ACGT[pos++] = codes[kmer & 3];
			}

		}
		~DumpOpt()
		{
			delete[]opt_ACGT;
		}
		
	}opt;

	void kmerToStr(CKmer<SIZE>& kmer, char* kmer_str)
	{
		//first byte
		char* base = opt.opt_ACGT + 4 * kmer.get_byte(kmer_bytes - 1) + 4 - in_first_byte;
		for(uint32 i = 0 ; i < in_first_byte ; ++i)
			*kmer_str++ = *base++;
		//rest
		for (int pos = kmer_bytes - 2; pos >= 0; --pos)
		{
			base = opt.opt_ACGT + 4 * kmer.get_byte(pos);
			*kmer_str++ = *base++;
			*kmer_str++ = *base++;
			*kmer_str++ = *base++;
			*kmer_str++ = *base++;
		}
	}

public:
	CDumpWriter(KMCDB& kmcdb) :kmcdb(kmcdb), output_desc(CConfig::GetInstance().output_desc), config(CConfig::GetInstance())
	{
		kmer_len = config.headers.front().kmer_len;
		kmer_bytes = (kmer_len + 3) / 4;
		in_first_byte = kmer_len % 4;
		if (in_first_byte == 0)
			in_first_byte = 4;
	}

	bool Process()
	{
		CKmer<SIZE> kmer;
		uint32 counter;
	
		uint32 counter_len;
		FILE* file = fopen(output_desc.file_src.c_str(), "wb");
		if (!file)
		{
			std::cout << "Error: cannot open file: " << output_desc.file_src << "\n";
			exit(1);
		}
		buf_pos = 0;
		buf_size = DUMP_BUF_SIZE;
		buf = new char[buf_size];

		//while (kmcdb.NextKmerSequential(kmer, counter))
		while (kmcdb.NextKmer(kmer, counter))
		{
			if (counter >= output_desc.cutoff_min && counter <= output_desc.cutoff_max)
			{
				kmerToStr(kmer, buf + buf_pos);
				buf[buf_pos + kmer_len] = '\t';
				counter_len = CNumericConversions::Int2PChar(counter, (uchar*)(buf + buf_pos + kmer_len + 1));
				buf[buf_pos + kmer_len + 1 + counter_len] = '\n';				
				buf_pos += kmer_len + 2 + counter_len;
				if (buf_pos + OVERHEAD_SIZE > buf_size)
				{
					fwrite(buf, 1, buf_pos, file);
					buf_pos = 0;
				}

			}
		}

		//save rest if necessary
		if (buf_pos)
		{
			fwrite(buf, 1, buf_pos, file);
			buf_pos = 0;
		}

		fclose(file);

		delete[] buf;

		return true;
	}
};

#endif