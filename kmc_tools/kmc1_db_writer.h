/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
*/

#ifndef _KMC1_DB_WRITER_H
#define _KMC1_DB_WRITER_H

#include "defs.h"
#include "config.h"
#include "queues.h"

#include <string>
#include <vector>

//************************************************************************************************************
// CKMC1SuffixFileWriter - thread for writing suffixes' parts
//************************************************************************************************************
class CKMC1SuffixFileWriter
{
public:
	CKMC1SuffixFileWriter(CSufWriteQueue& input_queue, FILE* kmc_suf) :
		input_queue(input_queue),
		kmc_suf(kmc_suf)
	{
	}
	void operator()()
	{
		uchar* buf;
		uint32 size;
		while (input_queue.pop(buf, size))
		{
			if (fwrite(buf, 1, size, kmc_suf) != size)
			{
				std::cerr << "Error while writing to kmc_suf file\n";
				exit(1);
			}
			delete[] buf;
		}
	}
private:
	CSufWriteQueue& input_queue;
	FILE* kmc_suf;
};

//************************************************************************************************************
// CKMC1DbWriter - writer of KMC1 database
//************************************************************************************************************
template<unsigned SIZE> class CKMC1DbWriter
{
public:
	CKMC1DbWriter(CBundle<SIZE>* bundle, COutputDesc& output_desc);
	~CKMC1DbWriter();
	bool Process();

private:
	static const uint32 PRE_BUFF_SIZE_BYTES = KMC1_DB_WRITER_PREFIX_BUFF_BYTES;
	static const uint32 SUF_BUFF_SIZE_BYTES = KMC1_DB_WRITER_SUFFIX_BUFF_BYTES;

	CConfig& config;
	COutputDesc& output_desc;
	CBundle<SIZE>* bundle = nullptr;
	FILE* kmc_pre, *kmc_suf;
	uint32 lut_prefix_len;
	uint32 current_prefix;
	uint32 counter_size;
	uint32 pre_buff_size;
	uint32 suf_buff_size;
	uint64* pre_buff;
	uchar* suf_buff;
	uint64 added_kmers;
	uint32 suffix_rec_bytes;
	uint32 suf_pos, pre_pos;


	void store_pre_buf();
	void send_suf_buf_to_queue();
	void start_writing();
	inline void add_kmer(CKmer<SIZE>& kmer, uint32 counter);
	void finish_writing();

	template<typename T> void write_header_part(T data);
	void calc_lut_prefix_len();


	CCircularQueue<SIZE> bundles_queue;
	CSufWriteQueue suf_buf_queue;


	//for simple and transform operations
	std::thread preparing_thread;
	CKMC1SuffixFileWriter* suffix_writer = nullptr;
	std::thread suf_buf_writing_thread;
public:
	void MultiOptputInit();
	void MultiOptputAddResultPart(COutputBundle<SIZE>& bundle);
	void MultiOptputAddResultPart(CBundle<SIZE>& bundle);
	void MultiOptputFinish();

};

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template <unsigned SIZE> CKMC1DbWriter<SIZE>::CKMC1DbWriter(CBundle<SIZE>* bundle, COutputDesc& output_desc) :
config(CConfig::GetInstance()),
output_desc(output_desc),
bundle(bundle),
bundles_queue(DEFAULT_CIRCULAL_QUEUE_CAPACITY)
{
	kmc_pre = NULL;
	kmc_suf = NULL;
	pre_buff = NULL;
	suf_buff = NULL;
	std::string kmc_pre_file_name = output_desc.file_src + ".kmc_pre";
	std::string kmc_suf_file_name = output_desc.file_src + ".kmc_suf";

	kmc_pre = fopen(kmc_pre_file_name.c_str(), "wb");
	setvbuf(kmc_pre, NULL, _IONBF, 0);

	if (!kmc_pre)
	{
		std::cerr << "Error: cannot open file : " << kmc_pre_file_name << "\n";
		exit(1);
	}
	kmc_suf = fopen(kmc_suf_file_name.c_str(), "wb");
	setvbuf(kmc_suf, NULL, _IONBF, 0);

	if (!kmc_suf)
	{
		fclose(kmc_pre);
		std::cerr << "Error: cannot open file : " << kmc_suf_file_name << "\n";
		exit(1);
	}

	setvbuf(kmc_pre, NULL, _IONBF, 0);
	setvbuf(kmc_suf, NULL, _IONBF, 0);
	// Calculate LUT size



	calc_lut_prefix_len();

	counter_size = MIN(BYTE_LOG(output_desc.counter_max), BYTE_LOG(output_desc.cutoff_max));
	if (output_desc.counter_value)
		counter_size = BYTE_LOG(output_desc.counter_value);
	suffix_rec_bytes = (config.kmer_len - lut_prefix_len) / 4 + counter_size;
	current_prefix = 0;
	added_kmers = 0;
	pre_buff_size = PRE_BUFF_SIZE_BYTES / sizeof(uint64);
	suf_buff_size = SUF_BUFF_SIZE_BYTES / suffix_rec_bytes;
	suf_pos = pre_pos = 0;

	pre_buff = new uint64[pre_buff_size];
	pre_buff[pre_pos++] = 0;
	suf_buff = new uchar[suf_buff_size * suffix_rec_bytes];


	suf_buf_queue.init(suf_buff_size * suffix_rec_bytes, SUFFIX_WRITE_QUEUE_CAPACITY);

}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbWriter<SIZE>::Process()
{

	start_writing();

	//Converts bundles to output buffers, suffix buffer is placed to another queue and write in separate thread (suffix_writer)
	preparing_thread = std::thread([this]{
		CBundleData<SIZE> bundle_data;
		while (bundles_queue.pop(bundle_data))
		{
			while (!bundle_data.Empty())
			{
				add_kmer(bundle_data.TopKmer(), bundle_data.TopCounter());
				bundle_data.Pop();
			}
		}
		suf_buf_queue.push(suf_buff, suffix_rec_bytes * suf_pos);
		suf_buf_queue.mark_completed();
	});



	suffix_writer = new CKMC1SuffixFileWriter(suf_buf_queue, kmc_suf);
	suf_buf_writing_thread = std::thread(std::ref(*suffix_writer));

#ifdef ENABLE_LOGGER
	CTimer timer;

#endif
	while (!bundle->Finished())
	{
#ifdef ENABLE_LOGGER
		timer.start();
#endif
		bundles_queue.push(bundle->Data());
#ifdef ENABLE_LOGGER
		CLoger::GetLogger().log_operation("dodawanie do kolejki wyjsciowej bundla", this, timer.get_time());
#endif
	}

	bundles_queue.mark_completed();

	preparing_thread.join();
	suf_buf_writing_thread.join();

	finish_writing();

	delete suffix_writer;
	return true;
}


template<unsigned SIZE> void CKMC1DbWriter<SIZE>::MultiOptputInit()
{
	start_writing();
	//Converts bundles to output buffers, suffix buffer is placed to another queue and write in separate thread (suffix_writer)
	preparing_thread = std::thread([this]{
		CBundleData<SIZE> bundle_data;
		while (bundles_queue.pop(bundle_data))
		{
			while (!bundle_data.Empty())
			{
				add_kmer(bundle_data.TopKmer(), bundle_data.TopCounter());
				bundle_data.Pop();
			}
		}
		suf_buf_queue.push(suf_buff, suffix_rec_bytes * suf_pos);
		suf_buf_queue.mark_completed();
	});

	suffix_writer = new CKMC1SuffixFileWriter(suf_buf_queue, kmc_suf);
	suf_buf_writing_thread = std::thread(std::ref(*suffix_writer));
}

template<unsigned SIZE> void CKMC1DbWriter<SIZE>::MultiOptputAddResultPart(COutputBundle<SIZE>& bundle)
{
	bundles_queue.push(bundle.Data());
}

template<unsigned SIZE> void CKMC1DbWriter<SIZE>::MultiOptputAddResultPart(CBundle<SIZE>& bundle)
{
	bundles_queue.push(bundle.Data());
}

template<unsigned SIZE> void CKMC1DbWriter<SIZE>::MultiOptputFinish()
{
	bundles_queue.mark_completed();
	preparing_thread.join();
	suf_buf_writing_thread.join();
	delete suffix_writer;
	finish_writing();
}

/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC1DbWriter<SIZE>::~CKMC1DbWriter()
{
	delete[] suf_buff;
	delete[] pre_buff;
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template <unsigned SIZE> template <typename T> void CKMC1DbWriter<SIZE>::write_header_part(T data)
{
	for (uint32 i = 0; i < sizeof(T); ++i)
	{
		char c = (data >> (i << 3)) & 0xff;
		if (putc(c, kmc_pre) == EOF)
		{
			std::cerr << "Error while writing header of kmc1\n";
			exit(1);
		}
	}
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbWriter<SIZE>::start_writing()
{
	if (fwrite("KMCP", 1, 4, kmc_pre) != 4)
	{
		std::cerr << "Error while writing starting KMCP marker";
		exit(1);
	}
	if (fwrite("KMCS", 1, 4, kmc_suf) != 4)
	{
		std::cerr << "Error while writing starting KMCS marker";
		exit(1);
	}
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbWriter<SIZE>::finish_writing()
{
	uint32 max_prefix = (1 << 2 * lut_prefix_len);
	while (current_prefix < max_prefix - 1)
	{
		pre_buff[pre_pos++] = added_kmers;
		++current_prefix;
		if (pre_pos == pre_buff_size)
			store_pre_buf();
	}
	store_pre_buf();
	send_suf_buf_to_queue();

	//store header
	write_header_part(config.kmer_len);
	write_header_part(config.headers.front().mode);
	write_header_part(counter_size);
	write_header_part(lut_prefix_len);
	write_header_part(output_desc.cutoff_min);
	write_header_part(output_desc.cutoff_max);
	write_header_part(added_kmers);

	bool both_stands = false;
	for (auto& input : config.headers)
		both_stands = both_stands || input.both_strands; //if any input database is in both strands, output is also in both strands

	write_header_part(!both_stands);


	for (uint32 i = 0; i < 31; ++i)
		write_header_part(uchar(0));

	write_header_part((uint32)64);


	if (fwrite("KMCP", 1, 4, kmc_pre) != 4)
	{
		std::cerr << "Error while writing end KMCP marker";
		exit(1);
	}
	if (fwrite("KMCS", 1, 4, kmc_suf) != 4)
	{
		std::cerr << "Error while writing end KMCS marker";
		exit(1);
	}
	fclose(kmc_pre);
	fclose(kmc_suf);
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbWriter<SIZE>::add_kmer(CKmer<SIZE>& kmer, uint32 counter)
{
	//if specific counter value is set use it as counter value (set_counts operation), do not check if counter is valid in term of cutoffs and counter max
	if (output_desc.counter_value)
		counter = static_cast<uint32>(output_desc.counter_value);
	else
	{
		if (counter < output_desc.cutoff_min || counter > output_desc.cutoff_max)
			return;
		if (counter > output_desc.counter_max)
			counter = output_desc.counter_max;
	}
	uint64 kmer_prefix = kmer.remove_suffix((config.kmer_len - lut_prefix_len) * 2);
	while (current_prefix < kmer_prefix)
	{
		pre_buff[pre_pos++] = added_kmers;
		++current_prefix;
		if (pre_pos == pre_buff_size)
			store_pre_buf();
	}
	uchar* rec = suf_buff + suf_pos * suffix_rec_bytes;

	kmer.store(rec, suffix_rec_bytes - counter_size);
	for (uint32 i = 0; i < counter_size; ++i)
		*rec++ = counter >> (i << 3);
	++suf_pos;
	if (suf_pos == suf_buff_size)
		send_suf_buf_to_queue();
	++added_kmers;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbWriter<SIZE>::store_pre_buf()
{
	if (fwrite(pre_buff, sizeof(uint64), pre_pos, kmc_pre) != pre_pos)
	{
		std::cerr << "Error while writing to kmc_pre file\n";
		exit(1);
	}
	pre_pos = 0;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbWriter<SIZE>::send_suf_buf_to_queue()
{
	suf_buf_queue.push(suf_buff, suffix_rec_bytes * suf_pos);
	suf_pos = 0;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbWriter<SIZE>::calc_lut_prefix_len()
{


	std::vector<uint32> best_lut_prefix_len_inputs(config.headers.size());


	for (uint32 i = 0; i < config.headers.size(); ++i)
	{
		uint32 best_lut_prefix_len = 0;
		uint64 best_mem_amount = 1ull << 62;
		for (lut_prefix_len = 1; lut_prefix_len < 16; ++lut_prefix_len)
		{
			uint32 suffix_len = config.headers[i].kmer_len - lut_prefix_len;
			if (suffix_len % 4)
				continue;

			uint64 suf_mem = config.headers[i].total_kmers * suffix_len / 4;
			uint64 lut_mem = (1ull << (2 * lut_prefix_len)) * sizeof(uint64);

			if (suf_mem + lut_mem < best_mem_amount)
			{
				best_lut_prefix_len = lut_prefix_len;
				best_mem_amount = suf_mem + lut_mem;
			}
		}
		best_lut_prefix_len_inputs[i] = best_lut_prefix_len;
	}

	//TODO: poki co jako lut size biore najwieszy z najlepszych dla baz wejsciowych
	lut_prefix_len = *std::max_element(best_lut_prefix_len_inputs.begin(), best_lut_prefix_len_inputs.end());
}

#endif


// ***** EOF