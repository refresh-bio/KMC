/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _KMC1_DB_READER_H
#define _KMC1_DB_READER_H
#include "kmer.h"
#include "defs.h"
#include "config.h"
#include "bundle.h"
#include "kmc_header.h"
#include "queues.h"
#include <iostream>
#include <cstring>
#include <thread>

enum class KMCDBOpenMode { sequential, sorted, counters_only };

//************************************************************************************************************
// CKMC1DbReader - reader of KMC1 database
//************************************************************************************************************
template<unsigned SIZE> class CKMC1DbReader : public CInput<SIZE>
{
public:
	CKMC1DbReader(const CKMC_header& header, const CInputDesc& desc, CPercentProgress& percent_progress, KMCDBOpenMode open_mode);

	void NextBundle(CBundle<SIZE>& bundle) override
	{
		bool exists = circular_queue->pop(bundle.Data());

		percent_progress.UpdateItem(progress_id, bundle.Size());
		
		if (exists)
			return;
		
		percent_progress.Complete(progress_id);

		this->finished = true;
		this->sorted_access_thread.join();
		delete this->circular_queue;		
	}

	void IgnoreRest() override
	{
		circular_queue->force_finish();
		this->finished = true;
		this->sorted_access_thread.join();
		delete this->circular_queue;
	}

	~CKMC1DbReader()
	{
		if(prefix_file != nullptr)
			fclose(prefix_file);
		if(sufix_file != nullptr)
			fclose(sufix_file);
		delete[] prefix_buff;
		delete[] sufix_buff;
	}
	
	bool NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter)
	{		
		if (next_kmer_sorted(kmer, counter))
		{
			percent_progress.UpdateItem(progress_id);
			return true;
		}
		percent_progress.Complete(progress_id);
		return false;
	}

	bool NextCounter(uint32& counter);

private:
	static const uint32 PREFIX_BUFF_BYTES = KMC1_DB_READER_PREFIX_BUFF_BYTES; 
	static const uint32 SUFIX_BUFF_BYTES = KMC1_DB_READER_SUFIX_BUFF_BYTES;
	const CKMC_header& header;
	const CInputDesc& desc;

	CPercentProgress& percent_progress;
	KMCDBOpenMode open_mode;

	uint32 progress_id;

	FILE* prefix_file;
	FILE* sufix_file;

	uint32 record_size; //of sufix, in bytes	
	uint32 current_preffix;
	uint32 sufix_bytes;
	uint64* prefix_buff = nullptr;
	uchar* sufix_buff = nullptr;

	uint32 prefix_bytes;
	uint32 kmer_bytes;

	uint64 prefix_buff_size;
	uint64 sufix_buff_size;

	uint64 prefix_buff_pos;
	uint64 sufix_buff_pos;

	uint64 prefix_left_to_read;
	uint64 sufix_left_to_read;

	std::string prefix_file_name;
	std::string sufix_file_name;

	uint64 sufix_number;

	CCircularQueue<SIZE>* circular_queue = nullptr; //for sorted access only
	std::thread sorted_access_thread;

	void reload_pref_buff();

	bool reload_suf_buff();

	bool next_kmer_sorted(CKmer<SIZE>& kmer, uint32& counter);

	void open_files();

	void allocate_buffers()
	{
		sufix_buff = new uchar[sufix_buff_size];
		if (open_mode == KMCDBOpenMode::sequential || open_mode == KMCDBOpenMode::sorted)
			prefix_buff = new uint64[prefix_buff_size];
	}
};

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

template<unsigned SIZE> CKMC1DbReader<SIZE>::CKMC1DbReader(const CKMC_header& header, const CInputDesc& desc, CPercentProgress& percent_progress, KMCDBOpenMode open_mode) :
	header(header), desc(desc), percent_progress(percent_progress), open_mode(open_mode)
{
	progress_id = percent_progress.RegisterItem(header.total_kmers);

	prefix_file = sufix_file = nullptr;
	sufix_bytes = (header.kmer_len - header.lut_prefix_len) / 4;
	record_size = sufix_bytes + header.counter_size;
	sufix_buff_size = SUFIX_BUFF_BYTES / record_size * record_size;
	prefix_buff_size = PREFIX_BUFF_BYTES / sizeof(uint64);

	sufix_left_to_read = header.total_kmers * record_size;

	if (sufix_left_to_read < sufix_buff_size)
		sufix_buff_size = sufix_left_to_read;

	prefix_left_to_read = (1 << header.lut_prefix_len * 2) - 1;

	if (prefix_left_to_read < prefix_buff_size)
		prefix_buff_size = prefix_left_to_read;

	prefix_bytes = (header.lut_prefix_len + 3) / 4;

	kmer_bytes = prefix_bytes + sufix_bytes;

	open_files();
	allocate_buffers();
	if(open_mode == KMCDBOpenMode::sequential || open_mode == KMCDBOpenMode::sorted)
		reload_pref_buff();
	reload_suf_buff();

	current_preffix = 0;
	sufix_number = 0;

	if (open_mode == KMCDBOpenMode::sorted)
	{
		circular_queue = new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY);
		sorted_access_thread = std::thread([this]{

			CKmer<SIZE> kmer;
			uint32 counter;
			CBundleData<SIZE> bundle_data;

			while (next_kmer_sorted(kmer, counter))
			{
				bundle_data.Insert(kmer, counter);
				if (bundle_data.Full())
				{
					if (!this->circular_queue->push(bundle_data))
						break;
				}
			}
			if (!bundle_data.Empty())
				this->circular_queue->push(bundle_data);
			this->circular_queue->mark_completed();		
		});
	}
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbReader<SIZE>::NextCounter(uint32& counter)
{
	while (true)
	{
		if (sufix_number >= header.total_kmers)
			return false;

		uchar* record = sufix_buff + sufix_buff_pos + sufix_bytes;

		counter = 0;
		for (int32 i = header.counter_size - 1; i >= 0; --i)
		{
			counter <<= 8;
			counter += record[i];
		}

		++sufix_number;
		sufix_buff_pos += record_size;

		if (sufix_buff_pos >= sufix_buff_size)
			reload_suf_buff();

		if (counter >= desc.cutoff_min && counter <= desc.cutoff_max)
			return true;
	}
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbReader<SIZE>::reload_pref_buff()
{
	uint64 to_read = MIN(prefix_left_to_read, prefix_buff_size);
	prefix_buff_pos = 0;
	if (to_read == 0)
	{
		prefix_buff[0] = header.total_kmers;//guard		
		return;
	}

	if (fread(prefix_buff, sizeof(uint64), to_read, prefix_file) != to_read)
	{
		std::cout << "Error: some error while reading " << prefix_file_name << "\n";
		exit(1);
	}	
	prefix_left_to_read -= to_read;
	if (to_read < prefix_buff_size)
	{
		prefix_buff[to_read] = header.total_kmers;//guard
	}
}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbReader<SIZE>::reload_suf_buff()
{
	uint64 to_read = MIN(sufix_left_to_read, sufix_buff_size);
	if (to_read == 0)
		return false;
	uint64 readed = fread(sufix_buff, 1, to_read, sufix_file);
	if (readed != to_read)
	{
		std::cout << "Error: some error while reading " << sufix_file_name << "\n";
		exit(1);
	}
	sufix_buff_pos = 0;
	sufix_left_to_read -= to_read;
	return true;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbReader<SIZE>::open_files()
{
	
	sufix_file_name = desc.file_src + ".kmc_suf";
	
	sufix_file = fopen(sufix_file_name.c_str(), "rb");
	setvbuf(sufix_file, NULL, _IONBF, 0);
	
	if (!sufix_file)
	{
		std::cout << "Error: cannot open file: " << sufix_file_name << "\n";
		exit(1);
	}

	char marker[4];
	if (fread(marker, 1, 4, sufix_file) != 4)
	{
		std::cout << "Error: while reading start marker in file: " << sufix_file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCS", 4) != 0)
	{
		std::cout << "Error: wrong start marker in file: " << sufix_file_name << "\n";
		exit(1);
	}


	my_fseek(sufix_file, -4, SEEK_END);
	if (fread(marker, 1, 4, sufix_file) != 4)
	{
		std::cout << "Error: while reading end marker in file: " << sufix_file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCS", 4) != 0)
	{
		std::cout << "Error: wrong end marker in file: " << sufix_file_name << "\n";
		exit(1);
	}
	my_fseek(sufix_file, 4, SEEK_SET); //skip KMCS

	if (open_mode == KMCDBOpenMode::sequential || open_mode == KMCDBOpenMode::sorted)
	{
		prefix_file_name = desc.file_src + ".kmc_pre";
		
		prefix_file = fopen(prefix_file_name.c_str(), "rb");
		setvbuf(prefix_file, NULL, _IONBF, 0);
		
		if (!prefix_file)
		{
			std::cout << "Error: cannot open file: " << prefix_file_name << "\n";
			exit(1);
		}
		my_fseek(prefix_file, 4 + sizeof(uint64), SEEK_SET);//skip KMCP and first value as it must be 0
	}
}


/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbReader<SIZE>::next_kmer_sorted(CKmer<SIZE>& kmer, uint32& counter)
{
	while (true)
	{
		if (sufix_number >= header.total_kmers)
			return false;

		while (prefix_buff[prefix_buff_pos] <= sufix_number)
		{
			++current_preffix;
			++prefix_buff_pos;
			if (prefix_buff_pos >= prefix_buff_size)
				reload_pref_buff();
		}

		uchar* record = sufix_buff + sufix_buff_pos;
		uint32 pos = kmer_bytes - 1;

		kmer.load(record, sufix_bytes);
		for (int32 i = prefix_bytes - 1; i >= 0; --i)
			kmer.set_byte(pos--, current_preffix >> (i << 3));

		counter = 0;
		for (int32 i = header.counter_size - 1; i >= 0; --i)
		{
			counter <<= 8;
			counter += record[i];
		}

		++sufix_number;
		sufix_buff_pos += record_size;
		
		if (sufix_buff_pos >= sufix_buff_size)
			reload_suf_buff();

		if (counter >= desc.cutoff_min && counter <= desc.cutoff_max)
			return true;
	}
}


#endif

// ***** EOF