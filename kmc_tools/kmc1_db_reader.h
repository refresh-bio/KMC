/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
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
#include <tuple>

enum class KMCDBOpenMode { sequential, sorted, counters_only };

class CSuffBufQueue
{
	using desc_t = std::tuple<uchar*, uint64>;
	std::vector<desc_t> data;
	int start = 0;
	int end = 0;
	bool is_completed = false;
	bool forced_to_finish = false;
	bool full;
	mutable std::mutex mtx;

	std::condition_variable cv_push;
	std::condition_variable cv_pop;
public:
	CSuffBufQueue(uint32 n_recs, uint64 buff_size) :
		is_completed(false), full(false)
	{
		data.resize(n_recs);
		for (auto& e : data)
		{
			std::get<0>(e) = new uchar[buff_size];
			std::get<1>(e) = 0;
		}
	}
	~CSuffBufQueue()
	{
		for (auto& e : data)
			delete std::get<0>(e);
	}
	bool pop(uchar* &suffix_buff, uint64& size)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] { return start != end || full || is_completed || forced_to_finish; });

		if (forced_to_finish)
			return false;

		if (is_completed && !full && start == end)
			return false;

		bool was_full = full;
		std::swap(suffix_buff, std::get<0>(data[start]));
		size = std::get<1>(data[start]);
		start = (start + 1) % data.size();
		full = false;
		if (was_full)
			cv_push.notify_all();
		return true;
	}

	bool push(uchar* &suffix_buff, uint64 size)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_push.wait(lck, [this] { return !full || forced_to_finish;  });

		if (forced_to_finish)
			return false;

		bool was_empty = start == end;
		std::swap(std::get<0>(data[end]), suffix_buff);
		std::get<1>(data[end]) = size;

		end = (end + 1) % data.size();
		if (end == start)
			full = true;

		if (was_empty)
			cv_pop.notify_all();

		return true;
	}

	void mark_completed()
	{
		std::lock_guard<std::mutex> lck(mtx);
		is_completed = true;
		cv_pop.notify_all();
	}

	void force_finish()
	{
		std::lock_guard<std::mutex> lck(mtx);
		forced_to_finish = true;
		cv_pop.notify_all();
		cv_push.notify_all();
	}

};


class CSuffBuffReader
{
	uint64 suffix_left_to_read = 0;
	uint64 suffix_buf_size = 0;
	CSuffBufQueue& suff_buff_queue;
	uchar* suff_buff;
	FILE* suffix_file;
	const std::string& suffix_file_name;
public:
	CSuffBuffReader(uint64 suffix_left_to_read, uint64 suffix_buff_size, CSuffBufQueue& suff_buff_queue, FILE* suffix_file, const std::string& suffix_file_name) :
		suffix_left_to_read(suffix_left_to_read),
		suffix_buf_size(suffix_buff_size),
		suff_buff_queue(suff_buff_queue),
		suffix_file(suffix_file),
		suffix_file_name(suffix_file_name)
	{
		suff_buff = new uchar[suffix_buf_size];
	}

	~CSuffBuffReader()
	{
		delete[]  suff_buff;
	}

	void operator()()
	{
		while (true)
		{
			uint64 to_read = MIN(suffix_left_to_read, suffix_buf_size);
			if (to_read == 0)
			{
				suff_buff_queue.mark_completed();
				return;
			}
			uint64 readed = fread(suff_buff, 1, to_read, suffix_file);
			if (readed != to_read)
			{
				std::cerr << "Error: some error while reading " << suffix_file_name << "\n";
				exit(1);
			}
			suffix_left_to_read -= readed;
			if (!suff_buff_queue.push(suff_buff, readed))
				break;
		}
		suff_buff_queue.mark_completed();
	}
};

//************************************************************************************************************
// CKMC1DbReader - reader of KMC1 database
//************************************************************************************************************
template<unsigned SIZE>
class CKMC1DbReader : public CInput<SIZE>
{
public:
	CKMC1DbReader(const CKMC_header& header, const CInputDesc& desc, CPercentProgress& percent_progress, KMCDBOpenMode open_mode);

	void NextBundle(CBundle<SIZE>& bundle) override
	{
		bool exists = circular_queue.pop(bundle.Data());

		percent_progress.UpdateItem(progress_id, bundle.Size());

		if (exists)
			return;

		percent_progress.Complete(progress_id);

		this->finished = true;
		this->building_thr.join();
		this->suffix_builder_thread.join();
	}

	void IgnoreRest() override
	{
		if (this->finished)
			return;
		circular_queue.force_finish();
		suffix_queue.force_finish();
		this->suff_buff_queue->force_finish();
		this->finished = true;
		this->building_thr.join();
		this->suffix_builder_thread.join();
	}

	~CKMC1DbReader()
	{
		suff_buff_reader_th.join();
		delete suff_buff_reader;
		delete suff_buff_queue;
		if (prefix_file != nullptr)
			fclose(prefix_file);
		if (suffix_file != nullptr)
			fclose(suffix_file);
		delete[] prefix_buff;
		delete[] suffix_buff;
	}

	bool NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter)
	{
		if (next_kmer_sequential_bundle.Empty())
		{
			bool exists = circular_queue.pop(next_kmer_sequential_bundle);
			if (!exists)
			{
				building_thr.join();
				suffix_builder_thread.join();
				percent_progress.Complete(progress_id);
				return false;
			}
		}
		kmer = next_kmer_sequential_bundle.TopKmer();
		counter = next_kmer_sequential_bundle.TopCounter();
		next_kmer_sequential_bundle.Pop();
		percent_progress.UpdateItem(progress_id);
		return true;
	}

	bool NextCounter(uint32& counter);

private:
	static const uint32 PREFIX_BUFF_BYTES = KMC1_DB_READER_PREFIX_BUFF_BYTES;
	static const uint32 SUFFIX_BUFF_BYTES = KMC1_DB_READER_SUFFIX_BUFF_BYTES;
	const CKMC_header& header;
	uint32 counter_size;
	uint64 total_kmers;

	uint64 kmers_left_for_current_prefix;
	uint64 total_kmers_left;

	const CInputDesc& desc;

	CPercentProgress& percent_progress;
	KMCDBOpenMode open_mode;

	uint32 progress_id;

	FILE* prefix_file;
	FILE* suffix_file;

	uint32 record_size; //of suffix, in bytes	
	CKmer<SIZE> current_prefix;

	uint32 suffix_bytes;
	uint64* prefix_buff = nullptr;
	uchar* suffix_buff = nullptr;
	uchar* record = nullptr;

	uint32 prefix_bytes;
	uint32 kmer_bytes;

	uint64 prefix_buff_size;	
	uint64 suffix_buff_size;
	
	uint64 prefix_buff_pos;
	uint64 suffix_buff_pos;

	uint64 prefix_left_to_read;	

	std::string prefix_file_name;
	std::string suffix_file_name;

	uint64 suffix_number;

	CCircularQueue<SIZE> circular_queue;
	CCircularQueue<SIZE> suffix_queue;
	std::thread building_thr;
	std::thread suffix_builder_thread;

	CBundleData<SIZE> suffix_bundle_get; CBundleData<SIZE> bundle_data;
	CBundleData<SIZE> suffix_bundle_insert;
	CBundleData<SIZE> next_kmer_sequential_bundle;
	CBundleData<SIZE> next_counter_bundle;

	void reload_pref_buff();

	bool reload_suf_buff();

	bool fill_bundle();
	bool fill_suffix();

	void open_files();

	void allocate_buffers()
	{
		suffix_buff = new uchar[suffix_buff_size];
		prefix_buff = new uint64[prefix_buff_size];
	}

	CSuffBuffReader* suff_buff_reader = nullptr;
	CSuffBufQueue* suff_buff_queue = nullptr;
	std::thread suff_buff_reader_th;
};


/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

template<unsigned SIZE> CKMC1DbReader<SIZE>::CKMC1DbReader(const CKMC_header& header, const CInputDesc& desc, CPercentProgress& percent_progress, KMCDBOpenMode open_mode) :
	header(header),
	counter_size(header.counter_size),
	total_kmers(header.total_kmers),
	desc(desc),
	percent_progress(percent_progress),
	open_mode(open_mode),
	circular_queue(DEFAULT_CIRCULAL_QUEUE_CAPACITY),
	suffix_queue(DEFAULT_CIRCULAL_QUEUE_CAPACITY)
{
	progress_id = percent_progress.RegisterItem(header.total_kmers);

	prefix_file = suffix_file = nullptr;
	suffix_bytes = (header.kmer_len - header.lut_prefix_len) / 4;
	record_size = suffix_bytes + header.counter_size;
	suffix_buff_size = SUFFIX_BUFF_BYTES / record_size * record_size;
	prefix_buff_size = PREFIX_BUFF_BYTES / sizeof(uint64);

	uint64 suffix_left_to_read = header.total_kmers * record_size;

	if (suffix_left_to_read < suffix_buff_size)
		suffix_buff_size = suffix_left_to_read;

	

	prefix_left_to_read = (1 << header.lut_prefix_len * 2) - 1;

	if (prefix_left_to_read < prefix_buff_size)
		prefix_buff_size = prefix_left_to_read;

	prefix_bytes = (header.lut_prefix_len + 3) / 4;

	kmer_bytes = prefix_bytes + suffix_bytes;

	open_files();
	allocate_buffers();

	suff_buff_queue = new CSuffBufQueue(4, suffix_buff_size);
	suff_buff_reader = new CSuffBuffReader(suffix_left_to_read, suffix_buff_size, *suff_buff_queue, suffix_file, suffix_file_name);

	suff_buff_reader_th = std::thread(std::ref(*suff_buff_reader));

	reload_pref_buff();
	reload_suf_buff();

	total_kmers_left = header.total_kmers;

	current_prefix.clear();
	suffix_number = 0;


	//run threads
	suffix_builder_thread = std::thread([this]{
		while (fill_suffix())
		{
			if (!suffix_queue.push(suffix_bundle_insert))
				break;
		}
		suffix_queue.mark_completed();
	});

	building_thr = std::thread([this]{

		while (fill_bundle())
		{
			if (!circular_queue.push(bundle_data))
				break;
		}
		circular_queue.mark_completed();
	});
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbReader<SIZE>::NextCounter(uint32& counter)
{
	if (next_counter_bundle.Empty())
	{
		bool exists = circular_queue.pop(next_counter_bundle);
		if (!exists)
		{
			building_thr.join();
			suffix_builder_thread.join();
			percent_progress.Complete(progress_id);
			return false;
		}
	}
	counter = next_counter_bundle.TopCounter();
	next_counter_bundle.Pop();
	percent_progress.UpdateItem(progress_id);

	return true;
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
		std::cerr << "Error: some error while reading " << prefix_file_name << "\n";
		exit(1);
	}
	prefix_left_to_read -= to_read;
	if (to_read < prefix_buff_size)
	{
		prefix_buff[to_read] = header.total_kmers;//guard
	}
	kmers_left_for_current_prefix = prefix_buff[0];
}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbReader<SIZE>::reload_suf_buff()
{
	suffix_buff_pos = 0;
	uint64 size;
	if (!suff_buff_queue->pop(suffix_buff, size))
		return false;

	record = suffix_buff;
	suffix_buff_pos = size / record_size;	
	return true;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC1DbReader<SIZE>::open_files()
{

	suffix_file_name = desc.file_src + ".kmc_suf";

	suffix_file = fopen(suffix_file_name.c_str(), "rb");
	setvbuf(suffix_file, NULL, _IONBF, 0);

	if (!suffix_file)
	{
		std::cerr << "Error: cannot open file: " << suffix_file_name << "\n";
		exit(1);
	}
	setvbuf(suffix_file, NULL, _IONBF, 0);

	char marker[4];
	if (fread(marker, 1, 4, suffix_file) != 4)
	{
		std::cerr << "Error: while reading start marker in file: " << suffix_file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCS", 4) != 0)
	{
		std::cerr << "Error: wrong start marker in file: " << suffix_file_name << "\n";
		exit(1);
	}


	my_fseek(suffix_file, -4, SEEK_END);
	if (fread(marker, 1, 4, suffix_file) != 4)
	{
		std::cerr << "Error: while reading end marker in file: " << suffix_file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCS", 4) != 0)
	{
		std::cerr << "Error: wrong end marker in file: " << suffix_file_name << "\n";
		exit(1);
	}
	my_fseek(suffix_file, 4, SEEK_SET); //skip KMCS


	prefix_file_name = desc.file_src + ".kmc_pre";

	prefix_file = fopen(prefix_file_name.c_str(), "rb");
	setvbuf(prefix_file, NULL, _IONBF, 0);

	if (!prefix_file)
	{
		std::cerr << "Error: cannot open file: " << prefix_file_name << "\n";
		exit(1);
	}
	my_fseek(prefix_file, 4 + sizeof(uint64), SEEK_SET);//skip KMCP and first value as it must be 0

}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbReader<SIZE>::fill_suffix()
{
	uint64 local_suffix_buff_pos = suffix_buff_pos;
	uchar* local_record = record;

	uint32 bundle_size = suffix_bundle_insert.size;
	uint32 bundle_pos = suffix_bundle_insert.insert_pos;

	uint32 counter_bytes = this->counter_size;
	uint32 suffix_bytes = this->suffix_bytes;

	uint32 counter_mask = (uint32)((1ull << (counter_bytes << 3)) - 1);

	bool little_endian = CConfig::GetInstance().IsLittleEndian();

	if (local_suffix_buff_pos)
	{
		while (true)
		{
#ifdef DISABLE_FAST_LOAD
			suffix_bundle_insert.kmers_with_counters[bundle_pos].kmer.load(local_record, suffix_bytes);
#else
			suffix_bundle_insert.kmers_with_counters[bundle_pos].kmer.load_fast(local_record, suffix_bytes, little_endian);
#endif
			CCounterBuilder::build_counter(suffix_bundle_insert.kmers_with_counters[bundle_pos].counter, local_record, counter_bytes, counter_mask, little_endian);

			++bundle_pos;



			if (!--local_suffix_buff_pos)
			{
				if (!reload_suf_buff())
					break;
				local_suffix_buff_pos = suffix_buff_pos;
				local_record = record;
			}

			if (bundle_pos == bundle_size)
				break;

		}
	}
	suffix_buff_pos = local_suffix_buff_pos;
	record = local_record;
	suffix_bundle_insert.insert_pos = bundle_pos;
	if (!suffix_bundle_insert.Empty())
		return true;
	return false;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC1DbReader<SIZE>::fill_bundle()
{
	CKmer<SIZE> kmer;
	uint32 counter;

	uint32 cutoff_min = desc.cutoff_min;
	uint32 cutoff_max = desc.cutoff_max;
	uint32 cutoff_range = cutoff_max - cutoff_min;

	uint64 local_kmers_left_for_current_prefix = kmers_left_for_current_prefix;
	uint64 local_total_kmers_left = total_kmers_left;

	uint32 bundle_pos = 0;
	uint32 size = bundle_data.size;



	if (suffix_bundle_get.Empty())
	{
		suffix_queue.pop(suffix_bundle_get);
	}


	uint32 suffix_bundle_pos = suffix_bundle_get.get_pos;
	uint32 suffix_bundle_size = suffix_bundle_get.insert_pos;

	uint32 left_in_suffix_bundle = suffix_bundle_size - suffix_bundle_pos;

	uint32 suffix_bytes = this->suffix_bytes;

	const uint32 unroll = 8;
	uint32 iter = MIN(suffix_bundle_size - suffix_bundle_pos, size - bundle_pos) / unroll;


	//macro to force inlining
#define FILL_BUNDLE_LOOP_UNROLL_CODE 														\
	while (!local_kmers_left_for_current_prefix)			   								\
	{														   								\
	current_prefix.increment_at(suffix_bytes);			   									\
	uint64 tmp = prefix_buff[prefix_buff_pos];			   									\
	++prefix_buff_pos;									   									\
	if (prefix_buff_pos >= prefix_buff_size)												\
	reload_pref_buff();																		\
																							\
	local_kmers_left_for_current_prefix = prefix_buff[prefix_buff_pos] - tmp;				\
	}																						\
																							\
	counter = suffix_bundle_get.kmers_with_counters[suffix_bundle_pos++].counter;			\
	--local_kmers_left_for_current_prefix;													\
																							\
	if (counter - cutoff_min <= cutoff_range)												\
	{																						\
	kmer = suffix_bundle_get.kmers_with_counters[suffix_bundle_pos - 1].kmer;				\
	kmer.set_prefix(current_prefix, suffix_bytes);											\
	bundle_data.kmers_with_counters[bundle_pos].kmer = kmer;								\
	bundle_data.kmers_with_counters[bundle_pos++].counter = counter;						\
	}
#ifdef WIN32 //VS is better if unroll by hand
	for (uint32 ii = 0; ii < iter; ++ii)
	{
		FILL_BUNDLE_LOOP_UNROLL_CODE
			FILL_BUNDLE_LOOP_UNROLL_CODE
			FILL_BUNDLE_LOOP_UNROLL_CODE
			FILL_BUNDLE_LOOP_UNROLL_CODE
			FILL_BUNDLE_LOOP_UNROLL_CODE
			FILL_BUNDLE_LOOP_UNROLL_CODE
			FILL_BUNDLE_LOOP_UNROLL_CODE
			FILL_BUNDLE_LOOP_UNROLL_CODE
	}
#else //g++ is better when it is not unrolled...
	for(uint32 ii = 0 ; ii < iter * unroll ; ++ii)
	{
		FILL_BUNDLE_LOOP_UNROLL_CODE
	}

#endif

	local_total_kmers_left -= unroll * iter;
	left_in_suffix_bundle -= unroll * iter;

	left_in_suffix_bundle += 1; //for simpler condition
	if (bundle_pos == size)
	{
		kmers_left_for_current_prefix = local_kmers_left_for_current_prefix;
		total_kmers_left = local_total_kmers_left;

		suffix_bundle_get.get_pos = suffix_bundle_pos;
		suffix_bundle_get.insert_pos = suffix_bundle_size;
		bundle_data.insert_pos = bundle_pos;
		return true;
	}

	while (local_total_kmers_left)
	{
		while (!local_kmers_left_for_current_prefix)
		{
			current_prefix.increment_at(suffix_bytes);
			uint64 tmp = prefix_buff[prefix_buff_pos];
			++prefix_buff_pos;
			if (prefix_buff_pos >= prefix_buff_size)
				reload_pref_buff();

			local_kmers_left_for_current_prefix = prefix_buff[prefix_buff_pos] - tmp;
		}


		if (!--left_in_suffix_bundle)
		{
			suffix_queue.pop(suffix_bundle_get);
			suffix_bundle_pos = suffix_bundle_get.get_pos;
			suffix_bundle_size = suffix_bundle_get.insert_pos;
			left_in_suffix_bundle = suffix_bundle_size - suffix_bundle_pos;
		}

		counter = suffix_bundle_get.kmers_with_counters[suffix_bundle_pos++].counter;

		--local_kmers_left_for_current_prefix;
		--local_total_kmers_left;

		if (counter - cutoff_min <= cutoff_range)
		{
			kmer = suffix_bundle_get.kmers_with_counters[suffix_bundle_pos - 1].kmer;
			kmer.set_prefix(current_prefix, suffix_bytes);

			bundle_data.kmers_with_counters[bundle_pos].kmer = kmer;
			bundle_data.kmers_with_counters[bundle_pos].counter = counter;

			if (++bundle_pos == size)
			{
				kmers_left_for_current_prefix = local_kmers_left_for_current_prefix;
				total_kmers_left = local_total_kmers_left;

				suffix_bundle_get.get_pos = suffix_bundle_pos;
				suffix_bundle_get.insert_pos = suffix_bundle_size;
				bundle_data.insert_pos = bundle_pos;
				return true;
			}
		}
	}
	kmers_left_for_current_prefix = local_kmers_left_for_current_prefix;
	total_kmers_left = local_total_kmers_left;

	suffix_bundle_get.get_pos = suffix_bundle_pos;
	suffix_bundle_get.insert_pos = suffix_bundle_size;

	bundle_data.insert_pos = bundle_pos;
	if (!bundle_data.Empty())
		return true;
	return false;
}


#endif

// ***** EOF