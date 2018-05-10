/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
*/

#ifndef _KMC2_DB_READER_H
#define _KMC2_DB_READER_H

#include "config.h"
#include "bundle.h"
#include "queues.h"
#include <vector>
#include <mutex>
#include <memory>
#include <tuple>
//#include <stack>
#include <queue>

#include <algorithm>
#include <condition_variable>

//Forward declaration
template<unsigned SIZE> class CKMC2DbReaderSorted;

template<unsigned SIZE> class CBin;



struct CBinBuff //must be moveble
{
	uchar* buf;
	uint32 size;

	CBinBuff() :
		buf(nullptr), size(0)
	{
	}

	CBinBuff(uchar* buf, uint32 size) :buf(buf), size(size)
	{

	}

#ifdef WIN32
	CBinBuff& operator=(CBinBuff&& rhs) throw()
#else
	CBinBuff& operator=(CBinBuff&& rhs) noexcept
#endif
	{
		if (this != &rhs)
		{
			buf = rhs.buf;
			size = rhs.size;
			rhs.buf = nullptr;
			rhs.size = 0;
		}
		return *this;
	}

#ifdef WIN32
	CBinBuff(CBinBuff&& rhs) throw()
#else
	CBinBuff(CBinBuff&& rhs) noexcept
#endif
	{
		buf = rhs.buf;
		size = rhs.size;
		rhs.buf = nullptr;
		rhs.size = 0;
	}

	CBinBuff(const CBinBuff&) = delete;
	CBinBuff& operator=(const CBinBuff&) = delete;
};

template<unsigned SIZE> class CBinBufProvider
{
	std::vector<CBinBuff> internal_bufs;
	uchar *buf_bins, *buf_internal;
	uint32 bins_left_to_read = 0;
	uint32 max_bin_bytes;
	uint32 rec_size;

	using desc_t = std::tuple<uint64, uint64, bool>;//current_kmer, last_kmer, is_empty
	using to_read_t = std::tuple<uint32, uint64, uchar*, uint32>;//bin_id, file_pos, buffer to read, size to read

	std::vector<desc_t> desc;
	//std::stack<to_read_t, std::vector<to_read_t>> to_read;
	std::queue<to_read_t, std::list<to_read_t>> to_read;

	mutable std::mutex mtx;
	std::condition_variable cv_pop;
	std::condition_variable cv_get_next_to_read;

	bool forced_to_finish = false;

public:
	void init(std::vector<CBin<SIZE>>& bins);

	void pop(uint32 bin_id, CBinBuff& bin_buf)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this, bin_id]{return !std::get<2>(desc[bin_id]); });

		std::swap(bin_buf, internal_bufs[bin_id]);
		std::get<2>(desc[bin_id]) = true;

		uint64 kmers_left = std::get<1>(desc[bin_id]) - std::get<0>(desc[bin_id]);
		if (kmers_left)
		{
			uint32 kmers_to_read = (uint32)MIN(kmers_left, max_bin_bytes / rec_size);
			internal_bufs[bin_id].size = kmers_to_read * rec_size;
			bool was_empty = to_read.empty();
			to_read.push(std::make_tuple(bin_id, 4 + std::get<0>(desc[bin_id]) * rec_size, internal_bufs[bin_id].buf, internal_bufs[bin_id].size));
			std::get<0>(desc[bin_id]) += kmers_to_read;
			if (was_empty)
				cv_get_next_to_read.notify_all();
		}
		else
		{
			--bins_left_to_read;
			if (!bins_left_to_read)
				cv_get_next_to_read.notify_all();
		}
	}

	void notify_bin_filled(uint32 bin_id)
	{
		std::lock_guard<std::mutex> lck(mtx);
		std::get<2>(desc[bin_id]) = false;
		cv_pop.notify_all();
	}

	bool get_next_to_read(uint32& bin_id, uint64& file_pos, uchar* &buf, uint32& size)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_get_next_to_read.wait(lck, [this]{return !to_read.empty() || !bins_left_to_read || forced_to_finish; });
		if (forced_to_finish || (to_read.empty() && !bins_left_to_read))
			return false;

		std::tie(bin_id, file_pos, buf, size) = to_read.front();
		to_read.pop();
		return true;
	}

	void force_to_finish()
	{
		std::lock_guard<std::mutex> lck(mtx);
		forced_to_finish = true;
		cv_get_next_to_read.notify_all();
	}

	~CBinBufProvider()
	{
		delete[] buf_bins;
		delete[] buf_internal;
	}
};

template<unsigned SIZE>
class CSufBinReader
{
	CBinBufProvider<SIZE>& bin_provider;
	FILE* suf_file;
public:
	CSufBinReader(CBinBufProvider<SIZE>& bin_provider, FILE* suf_file) :
		bin_provider(bin_provider),
		suf_file(suf_file)
	{

	}
	void operator()()
	{
		uint32 bin_id;
		uint64 file_pos;
		uchar* buf;
		uint32 size;

		while (bin_provider.get_next_to_read(bin_id, file_pos, buf, size))
		{
			my_fseek(suf_file, file_pos, SEEK_SET);
			if (fread(buf, 1, size, suf_file) != size)
			{
				std::cerr << "Error while reading suffix file\n";
				exit(1);
			}
			bin_provider.notify_bin_filled(bin_id);
		}
	}
};

template<unsigned SIZE> class CKmerPQ;
template<unsigned SIZE> class CBin
{
public:
	CBin(uint32 bin_id, uint64* LUT, CKMC2DbReaderSorted<SIZE>& kmc2_db);
	bool NextKmer(CKmer<SIZE>& kmer, uint32& counter);

	uint64 get_kmer_number_start()
	{
		return kmer_number_start;
	}

	uint64 get_kmers_left()
	{
		return kmers_left;
	}

	uint64 get_kmer_number_end()
	{
		return kmer_number_end;
	}

	uint32 get_record_size()
	{
		return record_size;
	}

	void set_bin_buff(CBinBuff&& _bin_buff)
	{
		bin_buff = std::move(_bin_buff);
		pos = bin_buff.size; //force reload
	}


#ifdef WIN32
	//Because VS2013 does generate default move ctor here
	CBin(CBin&& o) throw() :
		bin_id(o.bin_id),
		bin_buff(std::move(o.bin_buff)),
		LUT(o.LUT),
		pos(o.pos),
		bin_provider(o.bin_provider),
		kmc2_db(o.kmc2_db),
		kmer_number_start(o.kmer_number_start), kmer_number_end(o.kmer_number_end),
		kmers_left(o.kmers_left),
		kmers_left_for_current_prefix(o.kmers_left_for_current_prefix),
		kmer_bytes(o.kmer_bytes), prefix_bytes(o.prefix_bytes), suffix_bytes(o.suffix_bytes), counter_size(o.counter_size), record_size(o.record_size),
		cutoff_range(o.cutoff_range), cutoff_min(o.cutoff_min),
		prefix_pos(o.prefix_pos),
		prefix(o.prefix),
		max_prefix(o.max_prefix),
		is_little_endian(o.is_little_endian),
		counter_mask(o.counter_mask)
	{

	}
#else
	//g++ generate here move ctor automatically
#endif




private:
	uint32 bin_id;
	CBinBuff bin_buff;
	uint64* LUT;
	uint32 pos = 0;
	CBinBufProvider<SIZE>& bin_provider;
	CKMC2DbReaderSorted<SIZE>& kmc2_db;
	void reload_suf_buf();
	uint64 kmer_number_start, kmer_number_end;
	uint64 kmers_left;
	uint64 kmers_left_for_current_prefix;
	uint32 kmer_bytes, prefix_bytes, suffix_bytes, counter_size, record_size;
	uint32 cutoff_range, cutoff_min;
	uint64 prefix_pos = 0;
	CKmer<SIZE> prefix;
	uint64 max_prefix;
	bool is_little_endian;
	uint32 counter_mask;


	friend class CKmerPQ<SIZE>;
};


template<unsigned SIZE> void CBinBufProvider<SIZE>::init(std::vector<CBin<SIZE>>& bins)
{
	uint64 start, end;
	uint64 needed_mem = 0;
	rec_size = bins.front().get_record_size();
	max_bin_bytes = SINGLE_BIN_BUFF_SIZE_FOR_DB2_READER / rec_size * rec_size;
	uint32 mem;

	internal_bufs.resize(bins.size());
	for (uint32 i = 0; i < bins.size(); ++i)
	{
		auto& b = bins[i];
		start = b.get_kmer_number_start();
		end = b.get_kmer_number_end();
		mem = (uint32)MIN((end - start) * rec_size, max_bin_bytes);

		internal_bufs[i] = CBinBuff(nullptr, mem);
		desc.push_back(std::make_tuple(start, end, true));
		needed_mem += mem;
	}

	bins_left_to_read = (uint32)bins.size();
	buf_bins = new uchar[needed_mem];
	buf_internal = new uchar[needed_mem];

	internal_bufs[0].buf = buf_internal;

	uchar* ptr = buf_bins;
	bins[0].set_bin_buff(CBinBuff(ptr, internal_bufs[0].size));

	for (uint32 i = 1; i < internal_bufs.size(); ++i)
	{
		internal_bufs[i].buf = internal_bufs[i - 1].buf + internal_bufs[i - 1].size;
		ptr += internal_bufs[i - 1].size;
		bins[i].set_bin_buff(CBinBuff(ptr, internal_bufs[i].size));
	}

	for (uint32 bin_id = 0; bin_id < desc.size(); ++bin_id)
	{
		uint64 kmers_left = std::get<1>(desc[bin_id]) - std::get<0>(desc[bin_id]);
		if (kmers_left)
		{
			uint32 kmers_to_read = (uint32)MIN(kmers_left, max_bin_bytes / rec_size);
			internal_bufs[bin_id].size = kmers_to_read * rec_size;
			to_read.push(std::make_tuple(bin_id, 4 + std::get<0>(desc[bin_id]) * rec_size, internal_bufs[bin_id].buf, internal_bufs[bin_id].size));
			std::get<0>(desc[bin_id]) += kmers_to_read;
		}
		else
		{
			--bins_left_to_read;
		}
	}
}

//************************************************************************************************************
// CKmerPQ - Priority Queue of k-mers - binary heap. K-mers from bins are processed by this priority queue
//************************************************************************************************************
template<unsigned SIZE> class CKmerPQ
{
public:
	CKmerPQ(uint32 _no_of_bins);
	inline void init_add(CBin<SIZE>* bin);

	void Process(CCircularQueue<SIZE>& output_queue);

	inline void reset();

private:
	using elem_t = std::pair<CKmer<SIZE>, uint32>;//kmer, desc_id
	using desc_t = std::pair<CBin<SIZE>*, uint32>;//bin, counter
	std::vector<elem_t> elems;
	std::vector<desc_t> descs;
	uint32 pos, desc_pos;
};

//************************************************************************************************************
// CMergerParent - Merger of k-mers produced by CMergerChilds
//************************************************************************************************************
template<unsigned SIZE> class CMergerParent
{
public:
	CMergerParent(std::vector<CCircularQueue<SIZE>*>& input_queues, CCircularQueue<SIZE>& output_queue, uint32 n_subthreads);
	void operator()();

private:
	void Process2Inputs();
	void ProcessMoreInputs();

	uint32 LastLowerPlus1(CBundleData<SIZE>& b, uint32 start, uint32 end, CKmer<SIZE> kmer);
	void ProcessWithSubthreads();

	std::vector<CBundleData<SIZE>> input_bundles;
	std::vector<CCircularQueue<SIZE>*>& input_queues;

	CBundleData<SIZE> output_bundle;
	CCircularQueue<SIZE>& output_queue;
	uint32 n_subthreads;
};

//************************************************************************************************************
// CMergerChild - Merger of k-mers from bins
//************************************************************************************************************
template<unsigned SIZE> class CMergerChild
{
	using bin_iter = typename std::vector<CBin<SIZE>>::iterator;
public:
	CMergerChild(bin_iter begin, bin_iter end, CCircularQueue<SIZE>& output_queue);
	void operator()();

private:
	std::vector<std::reference_wrapper<CBin<SIZE>>> bins;
	CCircularQueue<SIZE>& output_queue;
};

//************************************************************************************************************
// CKMC2DbReaderSorted - Produce k-mers in sorted order from KMC2 database
//************************************************************************************************************
template<unsigned SIZE> class CKMC2DbReaderSorted
{
public:
	CKMC2DbReaderSorted(const CKMC_header& header, const CInputDesc& desc);

	void NextBundle(CBundle<SIZE>& bundle, bool& finished);

	void IgnoreRest();

	~CKMC2DbReaderSorted();

private:
	//void get_suf_buf_part(uchar* &buf, uint64 start, uint32 size);

	const CKMC_header& header;
	const CInputDesc& desc;
	uint64* LUTS = nullptr;
	uint32 lut_size = 0;
	uint32 suffix_bytes;
	uint32 record_size;
	FILE* kmc_suf;

	friend class CBin<SIZE>;
	std::vector<CBin<SIZE>> bins;
	CBinBufProvider<SIZE> bin_provider;

	CSufBinReader<SIZE>* suf_bin_reader;
	std::thread suf_bin_reader_th;


	uint32 n_child_threads;
	uint32 n_parent_threads;

	CMergerParent<SIZE>* parent = nullptr;
	std::thread parent_thread;

	CCircularQueue<SIZE>* output_queue;
	std::vector<CCircularQueue<SIZE>*> childs_parent_queues;

	std::vector<CMergerChild<SIZE>*> childs;
	std::vector<std::thread> childs_threads;

	//mutable std::mutex mtx;
};

//************************************************************************************************************
// CKCM2DbReaderSeqCounter_Base - Base class for classes to access k-mers one by one (not sorted) or 
// for counters only from KMC2 database
//************************************************************************************************************
template <unsigned SIZE>
class CKCM2DbReaderSeqCounter_Base
{
protected:
	CKCM2DbReaderSeqCounter_Base(const CKMC_header& header, const CInputDesc& desc);
	~CKCM2DbReaderSeqCounter_Base();

	void open_files();
	bool reload_suf_buff();

	static const uint32 PREFIX_BUFF_BYTES = KMC2_DB_READER_PREFIX_BUFF_BYTES;
	static const uint32 SUFFIX_BUFF_BYTES = KMC2_DB_READER_SUFFIX_BUFF_BYTES;

	const CKMC_header& header;
	const CInputDesc& desc;

	uint32 suffix_bytes;
	uint32 record_size; //of suffix, in bytes
	uint64 suffix_buff_size, suffix_buff_pos, suffix_left_to_read;
	uint64 prefix_buff_size, prefix_buff_pos, prefix_left_to_read;
	uint64 suffix_number;

	uint32 kmer_bytes, prefix_bytes;

	uchar* suffix_buff = nullptr;

	FILE* suffix_file;
	std::string suffix_file_name;
};


//************************************************************************************************************
// CKMC2DbReaderSequential - Produce k-mers sequentialy from KMC2 database (they are not sorted!)
//************************************************************************************************************
template<unsigned SIZE>
class CKMC2DbReaderSequential : public CKCM2DbReaderSeqCounter_Base<SIZE>
{
public:
	CKMC2DbReaderSequential(const CKMC_header& header, const CInputDesc& desc);
	bool NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter);
	~CKMC2DbReaderSequential();

private:
	void allocate_buffers();
	void reload_pref_buff();

	uint32 signle_bin_size, map_size, map_size_bytes, no_of_bins;
	std::string prefix_file_name;
	FILE* prefix_file;
	uint64 current_prefix_index;
	uint64 prefix_mask;

	uint64* prefix_buff = nullptr;
	int a;
};

//************************************************************************************************************
// CKMC2DbReaderCountersOnly - Produce counters of k-mers from KMC2 database
//************************************************************************************************************
template<unsigned SIZE>
class CKMC2DbReaderCountersOnly : CKCM2DbReaderSeqCounter_Base<SIZE>
{
public:
	CKMC2DbReaderCountersOnly(const CKMC_header& header, const CInputDesc& desc);
	bool NextCounter(uint32& counter);

private:
	void allocate_buffers();
};

//************************************************************************************************************
// CKMC2DbReader - reader of KMC2 
//************************************************************************************************************
template<unsigned SIZE>
class CKMC2DbReader : public CInput<SIZE>
{
public:
	CKMC2DbReader(const CKMC_header& header, const CInputDesc& desc, CPercentProgress& percent_progress, KMCDBOpenMode open_mode);

	void NextBundle(CBundle<SIZE>& bundle) override;

	void IgnoreRest() override;

	bool NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter);
	bool NextCounter(uint32& counter);

private:
	CPercentProgress& percent_progress;
	uint32 progress_id;

	std::unique_ptr<CKMC2DbReaderSorted<SIZE>> db_reader_sorted;
	std::unique_ptr<CKMC2DbReaderSequential<SIZE>> db_reader_sequential;
	std::unique_ptr<CKMC2DbReaderCountersOnly<SIZE>> db_reader_counters_only;
};



/*****************************************************************************************************************************/
/**************************************************** CBin IMPLEMENTATION ****************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

template<unsigned SIZE> CBin<SIZE>::CBin(uint32 bin_id, uint64* LUT, CKMC2DbReaderSorted<SIZE>& kmc2_db) :
bin_id(bin_id),
LUT(LUT),
bin_provider(kmc2_db.bin_provider),
kmc2_db(kmc2_db),
suffix_bytes(kmc2_db.suffix_bytes),
counter_size(kmc2_db.header.counter_size),
max_prefix(kmc2_db.lut_size - 1)

{
	kmer_number_start = LUT[0];
	kmer_number_end = LUT[kmc2_db.lut_size];
	kmers_left = kmer_number_end - kmer_number_start;
	kmers_left_for_current_prefix = LUT[1] - LUT[0];
	prefix_bytes = (kmc2_db.header.lut_prefix_len + 3) / 4;
	kmer_bytes = prefix_bytes + suffix_bytes;

	cutoff_min = kmc2_db.desc.cutoff_min;
	cutoff_range = kmc2_db.desc.cutoff_max - kmc2_db.desc.cutoff_min;

	record_size = suffix_bytes + counter_size;

	prefix.clear();
	is_little_endian = CConfig::GetInstance().IsLittleEndian();
	counter_mask = (uint32)((1ull << (counter_size << 3)) - 1);
}


/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CBin<SIZE>::NextKmer(CKmer<SIZE>& kmer, uint32& counter)
{
	while (kmers_left)
	{
		if (pos >= bin_buff.size)
			reload_suf_buf();

		//skip empty
		while (!kmers_left_for_current_prefix)
		{
			++prefix_pos;
			prefix.increment_at(suffix_bytes);
			kmers_left_for_current_prefix = LUT[prefix_pos + 1] - LUT[prefix_pos];
		}

		uchar* record = bin_buff.buf + pos;

		kmer.load_fast(record, suffix_bytes, is_little_endian);
		kmer.set_prefix(prefix, suffix_bytes);

		CCounterBuilder::build_counter(counter, record, counter_size, counter_mask, is_little_endian);
		pos += record_size;

		--kmers_left;
		--kmers_left_for_current_prefix;
		if (counter - cutoff_min <= cutoff_range)
			return true;
	}
	return false;
}
/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CBin<SIZE>::reload_suf_buf()
{
	bin_provider.pop(bin_id, bin_buff);
	pos = 0;
}



/*****************************************************************************************************************************/
/************************************************** CKmerPQ IMPLEMENTATION ***************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

template<unsigned SIZE> CKmerPQ<SIZE>::CKmerPQ(uint32 _no_of_bins)
{
	elems.resize(_no_of_bins + 1);
	descs.resize(_no_of_bins + 1);
	pos = 1;
	desc_pos = 0;
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKmerPQ<SIZE>::reset()
{
	pos = 1;
	desc_pos = 0;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKmerPQ<SIZE>::Process(CCircularQueue<SIZE>& output_queue)
{
	if (pos <= 1)
	{
		output_queue.mark_completed();
		return;
	}
	CBundleData<SIZE> bundle_data;
	uint32 desc_id = 0;
	CBin<SIZE>* bin = descs[elems[1].second].first;
	CKmer<SIZE> kmer;
	uint32 counter;
	uchar* record = nullptr;

	uint32 suffix_bytes = bin->suffix_bytes;
	uint32 counter_size = bin->counter_size;
	uint32 counter_mask = bin->counter_mask;
	uint32 record_size = bin->record_size;
	uint32 cutoff_min = bin->cutoff_min;
	uint32 cutoff_range = bin->cutoff_range;

	bool endian = CConfig::GetInstance().IsLittleEndian();

	while (true)
	{
		if (pos <= 1)
			break;
		bundle_data.Insert(elems[1].first, descs[elems[1].second].second);


		//UPDATE HEAP!
		desc_id = elems[1].second;
		bin = descs[desc_id].first;

		bool exists = false;


		while (bin->kmers_left)
		{

			if (bin->pos >= bin->bin_buff.size)
				bin->reload_suf_buf();

			//skip empty
			while (!bin->kmers_left_for_current_prefix)
			{
				++bin->prefix_pos;
				bin->prefix.increment_at(suffix_bytes);
				bin->kmers_left_for_current_prefix = bin->LUT[bin->prefix_pos + 1] - bin->LUT[bin->prefix_pos];
			}

			record = bin->bin_buff.buf + bin->pos;

			kmer.load_fast(record, suffix_bytes, endian);
			kmer.set_prefix(bin->prefix, suffix_bytes);

			CCounterBuilder::build_counter(counter, record, counter_size, counter_mask, endian);
			bin->pos += record_size;

			--bin->kmers_left;
			--bin->kmers_left_for_current_prefix;
			if (counter - cutoff_min <= cutoff_range)
			{
				exists = true;
				break;
			}
		}

		if (!exists)
		{
			kmer.set(elems[--pos].first);
			desc_id = elems[pos].second;
		}
		else
			descs[desc_id].second = counter;

		uint32 parent, less;
		parent = less = 1;
		while (true)
		{
			if (parent * 2 >= pos)
				break;
			if (parent * 2 + 1 >= pos)
				less = parent * 2;
			else if (elems[parent * 2].first < elems[parent * 2 + 1].first)
				less = parent * 2;
			else
				less = parent * 2 + 1;
			if (elems[less].first < kmer)
			{
				elems[parent] = elems[less];
				parent = less;
			}
			else
				break;
		}
		elems[parent] = std::make_pair(kmer, desc_id);

		if (bundle_data.Full())
		{
			if (!output_queue.push(bundle_data))
				break;
		}
	}
	if (!bundle_data.Empty())
		output_queue.push(bundle_data);
	output_queue.mark_completed();
}

/*****************************************************************************************************************************/
template<unsigned SIZE> inline void CKmerPQ<SIZE>::init_add(CBin<SIZE>* bin)
{
	CKmer<SIZE> kmer;
	uint32 counter;
	if (bin->NextKmer(kmer, counter))
	{
		descs[desc_pos] = std::make_pair(bin, counter);
		elems[pos] = std::make_pair(kmer, desc_pos);
		uint32 child_pos = pos++;

		while (child_pos > 1 && elems[child_pos].first < elems[child_pos / 2].first)
		{
			swap(elems[child_pos], elems[child_pos / 2]);
			child_pos /= 2;
		}

		++desc_pos;
	}
}

/*****************************************************************************************************************************/
/*********************************************** CMergerParent IMPLEMENTATION ************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE>	CMergerParent<SIZE>::CMergerParent(std::vector<CCircularQueue<SIZE>*>& input_queues, CCircularQueue<SIZE>& output_queue, uint32 n_subthreads) :
	input_queues(input_queues),
	output_queue(output_queue),
	n_subthreads(n_subthreads)
{
	input_bundles.resize(input_queues.size());
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CMergerParent<SIZE>::operator()()
{
	if (n_subthreads > 1)
	{
		ProcessWithSubthreads();
	}

	else if (input_queues.size() == 2)
	{
		Process2Inputs();
	}
	else
	{
		ProcessMoreInputs();
	}
}

//************************************************************************************************************
// CParentSubthreadPartDesc - Contains current state of buffers
//************************************************************************************************************
struct CParentSubthreadPartDesc
{
	uint32 start, end, part_end;
	uint32 left()
	{
		return end - part_end;
	};
};

//************************************************************************************************************
// CParentSubthreadSynchronizer - Synchronize subthreads created by CMergerParent
//************************************************************************************************************
class CParentSubthreadSynchronizer
{
	uint32 n_tasks = 0;
	std::mutex mtx;
	std::condition_variable cv;
public:
	void decrement()
	{
		std::lock_guard<std::mutex> lck(mtx);
		--n_tasks;
	}
	void increment()
	{
		std::lock_guard<std::mutex> lck(mtx);
		++n_tasks;
	}

	void wait()
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return !n_tasks; });
	}

	void notify_task_finished()
	{
		std::lock_guard<std::mutex> lck(mtx);
		--n_tasks;
		if (!n_tasks)
			cv.notify_all();
	}
};

//************************************************************************************************************
// CParentSubthreadDesc - Input data of subthreads of CMergerParent
//************************************************************************************************************
template<unsigned SIZE>
struct CParentSubthreadDesc
{
	std::vector<CBundleData<SIZE>>* inputs;
	CBundleData<SIZE>* out;
	std::vector<CParentSubthreadPartDesc> desc;
	uint32 o_start;
};

//************************************************************************************************************
// CParentSubthreadDescQueue - Passes data to subthreads from CMergerParent
//************************************************************************************************************
template<unsigned SIZE>
class CParentSubthreadDescQueue
{
	mutable std::mutex mtx;
	std::condition_variable cv;
	bool empty = true;
	bool completed = false;
public:
	CParentSubthreadDesc<SIZE> desc;
	void start()
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return empty; });
		empty = false;
		cv.notify_all();
	}

	bool pop(CParentSubthreadDesc<SIZE> &_desc)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return !empty || completed; });
		if (completed)
			return false;
		_desc = desc;
		empty = true;
		cv.notify_all();
		return true;
	}
	void mark_completed()
	{
		std::lock_guard<std::mutex> lck(mtx);
		completed = true;
		cv.notify_all();
	}
};

//************************************************************************************************************
// CMergerParentSubthread - Merge data described in CParentSubthreadDescQueue
//************************************************************************************************************
template<unsigned SIZE>
class CMergerParentSubthread
{
	CParentSubthreadDescQueue<SIZE>& task_queue;
	CParentSubthreadSynchronizer& synchronizer;
public:
	CMergerParentSubthread(CParentSubthreadDescQueue<SIZE>& task_queue, CParentSubthreadSynchronizer& synchronizer)
		:
		task_queue(task_queue),
		synchronizer(synchronizer)
	{
	}

	void operator()()
	{
		CParentSubthreadDesc<SIZE> t;
		while (task_queue.pop(t))
		{
			std::vector<CBundleData<SIZE>>& inputs = *t.inputs;
			CBundleData<SIZE>& out = *t.out;
			using heap_elem_t = std::pair<CKmer<SIZE>, uint32>;
			std::vector<heap_elem_t> heap(inputs.size() + 1);

			uint32 pos = 1;

			for (uint32 i = 0; i < inputs.size(); ++i)
			{
				if (t.desc[i].start >= t.desc[i].part_end)
					continue;
				heap[pos] = std::make_pair(inputs[i].kmers_with_counters[t.desc[i].start].kmer, i);
				t.desc[i].start++;
				uint32 child_pos = pos++;
				while (child_pos > 1 && heap[child_pos].first < heap[child_pos / 2].first)
				{
					std::swap(heap[child_pos], heap[child_pos / 2]);
					child_pos /= 2;
				}
			}

			uint32 out_pos = t.o_start;
			while (true)
			{
				if (pos <= 1)
					break;
				uint32 desc_pos = heap[1].second;
				out.kmers_with_counters[out_pos++] = inputs[desc_pos].kmers_with_counters[t.desc[desc_pos].start - 1];

				CKmer<SIZE> kmer;
				uint32 desc_id = heap[1].second;

				if (t.desc[desc_pos].start < t.desc[desc_pos].part_end)
				{
					kmer.set(inputs[desc_pos].kmers_with_counters[t.desc[desc_pos].start].kmer);
					++t.desc[desc_pos].start;
				}
				else
				{
					kmer.set(heap[--pos].first);
					desc_id = heap[pos].second;
				}

				uint32 parent, less;
				parent = less = 1;
				while (true)
				{
					if (parent * 2 >= pos)
						break;
					if (parent * 2 + 1 >= pos)
						less = parent * 2;
					else if (heap[parent * 2].first < heap[parent * 2 + 1].first)
						less = parent * 2;
					else
						less = parent * 2 + 1;
					if (heap[less].first < kmer)
					{
						heap[parent] = heap[less];
						parent = less;
					}
					else
						break;
				}
				heap[parent] = std::make_pair(kmer, desc_id);
			}
			//out.insert_pos = out_pos;
			synchronizer.notify_task_finished();
		}
	}
};

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CMergerParent<SIZE>::Process2Inputs()
{
	CBundleData<SIZE> b1, b2;
	CCircularQueue<SIZE>* q1, *q2;
	q1 = input_queues[0];
	q2 = input_queues[1];
	bool q1_empty = !q1->pop(b1);
	bool q2_empty = !q2->pop(b2);

	if (q1_empty && q2_empty)
	{
		output_queue.mark_completed();
		return;
	}
	if (q1_empty || q2_empty)
	{
		CCircularQueue<SIZE>* q = q1_empty ? q2 : q1;
		CBundleData<SIZE>& b = q1_empty ? b2 : b1;
		while (true)
		{
			if (!output_queue.push(b))
				break;
			if (!q->pop(b))
				break;
		}
		output_queue.mark_completed();
		return;
	}

	uint32 get1 = 0;
	uint32 get2 = 0;

	CKmer<SIZE> kmer2 = b2.kmers_with_counters[get2].kmer;
	uint32 counter2 = b2.kmers_with_counters[get2].counter;
	CKmer<SIZE> kmer1 = b1.kmers_with_counters[get1].kmer;
	uint32 counter1 = b1.kmers_with_counters[get1].counter;

	uint32 left1 = b1.NRecLeft();
	uint32 left2 = b2.NRecLeft();

	uint32 out_insert_pos = 0;
	uint32 out_size = output_bundle.size;

	while (true)
	{
		if (kmer1 < kmer2)
		{
			output_bundle.kmers_with_counters[out_insert_pos].kmer = kmer1;
			output_bundle.kmers_with_counters[out_insert_pos++].counter = counter1;
			if (out_insert_pos == out_size)
			{
				output_bundle.insert_pos = out_insert_pos;

				if (!output_queue.push(output_bundle))
					break;
				out_insert_pos = 0;
				out_size = output_bundle.size;
			}

			++get1;
			if (--left1)
			{
				kmer1 = b1.kmers_with_counters[get1].kmer;
				counter1 = b1.kmers_with_counters[get1].counter;
			}
			else
			{
				b1.get_pos = get1;
				if (q1->pop(b1))
				{
					get1 = 0;
					kmer1 = b1.kmers_with_counters[get1].kmer;
					counter1 = b1.kmers_with_counters[get1].counter;
					left1 = b1.NRecLeft();
				}
				else
					break;

			}
		}
		else
		{
			output_bundle.kmers_with_counters[out_insert_pos].kmer = kmer2;
			output_bundle.kmers_with_counters[out_insert_pos++].counter = counter2;
			if (out_insert_pos == out_size)
			{
				output_bundle.insert_pos = out_insert_pos;
				if (!output_queue.push(output_bundle))
					break;
				out_insert_pos = 0;
				out_size = output_bundle.size;
			}

			++get2;
			if (--left2)
			{
				kmer2 = b2.kmers_with_counters[get2].kmer;
				counter2 = b2.kmers_with_counters[get2].counter;
			}
			else
			{
				b2.get_pos = get2;
				if (q2->pop(b2))
				{
					get2 = 0;
					kmer2 = b2.kmers_with_counters[get2].kmer;
					counter2 = b2.kmers_with_counters[get2].counter;
					left2 = b2.NRecLeft();
				}
				else
					break;
			}
		}
	}

	b1.get_pos = get1;
	b2.get_pos = get2;
	output_bundle.insert_pos = out_insert_pos;

	if (b1.Empty())
		q1->pop(b1);
	if (!b1.Empty())
	{
		while (true)
		{
			if (b1.Empty())
			{
				if (!q1->pop(b1))
					break;
			}
			output_bundle.Insert(b1.TopKmer(), b1.TopCounter());
			b1.Pop();
			if (output_bundle.Full())
			{
				if (!output_queue.push(output_bundle))
					break;
			}
		}
	}

	if (b2.Empty())
		q2->pop(b2);
	if (!b2.Empty())
	{
		while (true)
		{
			if (b2.Empty())
			{
				if (!q2->pop(b2))
					break;
			}
			output_bundle.Insert(b2.TopKmer(), b2.TopCounter());
			b2.Pop();
			if (output_bundle.Full())
			{
				if (!output_queue.push(output_bundle))
					break;
			}
		}
	}

	if (!output_bundle.Empty())
		output_queue.push(output_bundle);
	output_queue.mark_completed();
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CMergerParent<SIZE>::ProcessMoreInputs()
{
	//init
	auto q_iter = input_queues.begin();
	auto b_iter = input_bundles.begin();
	for (; q_iter != input_queues.end();)
	{
		if (!(*q_iter)->pop(*b_iter))
		{
			q_iter = input_queues.erase(q_iter);
			b_iter = input_bundles.erase(b_iter);
		}
		else
			++q_iter, ++b_iter;
	}
	uint32 index_of_min = 0;
	CKmer<SIZE> min_kmer;

	uint32 output_bundle_insert_pos = output_bundle.insert_pos;
	uint32 output_bundle_size = output_bundle.size;
	decltype(output_bundle.kmers_with_counters) kmers_counters = output_bundle.kmers_with_counters;

	while (input_bundles.size())
	{
		index_of_min = 0;
		min_kmer = input_bundles[index_of_min].TopKmer();

		if (input_bundles.size() == 4)
		{
			uint32 tmp_min = 2;
			CKmer<SIZE> tmp_kmer = input_bundles[tmp_min].TopKmer();

			if (input_bundles[1].TopKmer() < min_kmer)
			{
				min_kmer = input_bundles[1].TopKmer();
				index_of_min = 1;
			}
			if (input_bundles[3].TopKmer() < tmp_kmer)
			{
				tmp_kmer = input_bundles[3].TopKmer();
				tmp_min = 3;
			}
			if (tmp_kmer < min_kmer)
			{
				index_of_min = tmp_min;
			}
		}
		else if (input_bundles.size() == 3)
		{
			if (input_bundles[1].TopKmer() < min_kmer)
			{
				min_kmer = input_bundles[1].TopKmer();
				index_of_min = 1;
			}
			if (input_bundles[2].TopKmer() < min_kmer)
			{
				min_kmer = input_bundles[2].TopKmer();
				index_of_min = 2;
			}
		}
		else if (input_bundles.size() == 2)
		{
			if (input_bundles[1].TopKmer() < min_kmer)
			{
				min_kmer = input_bundles[1].TopKmer();
				index_of_min = 1;
			}
		}
		else if (input_bundles.size() == 1)
		{
		}
		else
		{
			for (uint32 i = 1; i < input_bundles.size(); ++i)
			{
				if (input_bundles[i].TopKmer() < min_kmer)
				{
					index_of_min = i;
					min_kmer = input_bundles[index_of_min].TopKmer();
				}
			}
		}

		kmers_counters[output_bundle_insert_pos].kmer = input_bundles[index_of_min].TopKmer();
		kmers_counters[output_bundle_insert_pos++].counter = input_bundles[index_of_min].TopCounter();
		//output_bundle.Insert(input_bundles[index_of_min].TopKmer(), input_bundles[index_of_min].TopCounter());

		input_bundles[index_of_min].Pop();
		if (input_bundles[index_of_min].Empty())
		{

			if (!input_queues[index_of_min]->pop(input_bundles[index_of_min]))
			{
				input_queues.erase(input_queues.begin() + index_of_min);
				input_bundles.erase(input_bundles.begin() + index_of_min);
			}
		}


		if (output_bundle_insert_pos == output_bundle_size)
			//if (output_bundle.Full())
		{
			output_bundle.insert_pos = output_bundle_insert_pos;
			if (!output_queue.push(output_bundle))
			{
				output_bundle_insert_pos = output_bundle.insert_pos; //0
				kmers_counters = output_bundle.kmers_with_counters;
				break;
			}
			output_bundle_insert_pos = output_bundle.insert_pos; //0
			kmers_counters = output_bundle.kmers_with_counters;
		}
	}
	output_bundle.insert_pos = output_bundle_insert_pos;
	if (!output_bundle.Empty())
		output_queue.push(output_bundle);
	output_queue.mark_completed();
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CMergerParent<SIZE>::ProcessWithSubthreads()
{
	std::vector<CBundleData<SIZE>> input_bundles(input_queues.size());
	
	//std::vector<CBundleData<SIZE>> output_bundles;
	//output_bundles.reserve(n_subthreads);
	//for (uint32 i = 0; i < n_subthreads; ++i)
		//output_bundles.emplace_back(KMC2_DB_READER_BUNDLE_CAPACITY);
	CBundleData<SIZE> output_bundle(KMC2_DB_READER_BUNDLE_CAPACITY);
	uint32 curr_end = output_bundle.insert_pos;

	auto q_iter = input_queues.begin();
	auto b_iter = input_bundles.begin();
	for (; q_iter != input_queues.end();)
	{
		if (!(*q_iter)->pop(*b_iter))
		{
			q_iter = input_queues.erase(q_iter);
			b_iter = input_bundles.erase(b_iter);
		}
		else
			++q_iter, ++b_iter;
	}

	std::vector<CParentSubthreadDescQueue<SIZE>> task_descs(n_subthreads);
	std::vector<CParentSubthreadPartDesc> descs(input_queues.size());
	for (uint32 i = 0; i < input_queues.size(); ++i)
	{
		descs[i].start = descs[i].part_end = input_bundles[i].get_pos;
		descs[i].end = input_bundles[i].insert_pos;
	}

	std::vector<CMergerParentSubthread<SIZE>> subtasks;
	std::vector<std::thread> subthreads;

	subtasks.reserve(n_subthreads);
	subthreads.reserve(n_subthreads);

	CParentSubthreadSynchronizer syncer;
	for (uint32 i = 0; i < n_subthreads; ++i)
	{
		subtasks.emplace_back(task_descs[i], syncer);
		subthreads.push_back(std::thread(std::ref(subtasks.back())));
	}

	while (input_bundles.size())
	{
		//prepare threads					
		for (uint32 th = 0; th < n_subthreads; ++th)
		{
			bool any_empty = false;
			for (uint32 i = 0; i < descs.size(); ++i)
			{
				descs[i].part_end = descs[i].start + MIN((output_bundle.size - curr_end) / input_bundles.size() / (n_subthreads - th), descs[i].left());

				//if any is empty it must be refilled or removed
				if (descs[i].part_end == descs[i].start)
				{
					any_empty = true;
					break;
				}
			}
			if (any_empty)
				break;

			//find min kmer
			uint32 min_kmer_i = 0;
			CKmer<SIZE> min_kmer = input_bundles[0].kmers_with_counters[descs[0].part_end - 1].kmer;

			for (uint32 i = min_kmer_i + 1; i < input_bundles.size(); ++i)
			{
				if (input_bundles[i].kmers_with_counters[descs[i].part_end - 1].kmer < min_kmer)
				{
					min_kmer = input_bundles[i].kmers_with_counters[descs[i].part_end - 1].kmer;
					min_kmer_i = i;
				}
			}

			uint32 prev_end = curr_end;

			//correct part_end according to min kmer
			for (uint32 i = 0; i < descs.size(); ++i)
			{
				if (i != min_kmer_i)
				{
					descs[i].part_end = LastLowerPlus1(input_bundles[i], descs[i].start, descs[i].part_end, min_kmer);
				}
				curr_end += descs[i].part_end - descs[i].start;
			}

			task_descs[th].desc.desc = descs;
			task_descs[th].desc.inputs = &input_bundles;
			task_descs[th].desc.out = &output_bundle;
			task_descs[th].desc.o_start = prev_end;
			syncer.increment();
			task_descs[th].start();

			for (uint32 i = 0; i < descs.size(); ++i)
			{
				descs[i].start = descs[i].part_end;
			}
		}

		syncer.wait(); //BARIER

		//send output
		output_bundle.insert_pos = curr_end;
		if (!output_bundle.Empty())
		{
			if (!output_queue.push(output_bundle))
			{
				break;
			}
		}

		curr_end = output_bundle.insert_pos;

		auto q_iter = input_queues.begin();
		auto b_iter = input_bundles.begin();
		auto d_iter = descs.begin();

		for (; b_iter != input_bundles.end();)
		{
			b_iter->get_pos = d_iter->start;
			if ((*b_iter).Empty())
			{
				if ((*q_iter)->pop(*b_iter))
				{
					d_iter->start = d_iter->part_end = b_iter->get_pos;
					d_iter->end = b_iter->insert_pos;
				}
				else
				{
					d_iter = descs.erase(d_iter);
					b_iter = input_bundles.erase(b_iter);
					q_iter = input_queues.erase(q_iter);
					continue;
				}
			}
			++q_iter, ++b_iter, ++d_iter;
		}

	}

	for (auto& t : task_descs)
		t.mark_completed();

	for (auto& t : subthreads)
		t.join();

	output_queue.mark_completed();
}

/*****************************************************************************************************************************/
template<unsigned SIZE> uint32 CMergerParent<SIZE>::LastLowerPlus1(CBundleData<SIZE>& b, uint32 start, uint32 end, CKmer<SIZE> kmer)
{
	auto ub = std::upper_bound(b.kmers_with_counters + start, b.kmers_with_counters + end, kmer,
		[](const CKmer<SIZE>& kmer, typename CBundleData<SIZE>::CKmerWithCounter& kc)
	{
		return kmer < kc.kmer;
	});
	return ub - b.kmers_with_counters;
}

/*****************************************************************************************************************************/
/************************************************ CMergerChild IMPLEMENTATION ************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE> CMergerChild<SIZE>::CMergerChild(bin_iter begin, bin_iter end, CCircularQueue<SIZE>& output_queue) :
bins(begin, end),
output_queue(output_queue)
{

}
/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CMergerChild<SIZE>::operator()()
{
	CKmerPQ<SIZE> kmers_pq(static_cast<uint32>(bins.size()));
	for (uint32 i = 0; i < bins.size(); ++i)
		kmers_pq.init_add(&bins[i].get());

	kmers_pq.Process(output_queue);
}



/*****************************************************************************************************************************/
/********************************************* CKMC2DbReaderSorted IMPLEMENTATION ********************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC2DbReaderSorted<SIZE>::CKMC2DbReaderSorted(const CKMC_header& header, const CInputDesc& desc) :
header(header),
desc(desc)
{
	LUTS = nullptr;
	lut_size = 1 << 2 * header.lut_prefix_len;
	uint32 lut_recs = (1 << 2 * header.lut_prefix_len) * header.no_of_bins + 1;
	LUTS = new uint64[lut_recs];
	suffix_bytes = (header.kmer_len - header.lut_prefix_len) / 4;
	record_size = suffix_bytes + header.counter_size;
	if (!LUTS)
	{
		std::cerr << "Error: cannot allocate memory for LUTS of KMC2 database\n";
		exit(1);
	}

	std::string kmc_pre_file_name = desc.file_src + ".kmc_pre";
	FILE* kmc_pre = fopen(kmc_pre_file_name.c_str(), "rb");
	if (!kmc_pre)
	{
		std::cerr << "Error: cannot open kmc2 prefix file to read LUTS";
		exit(1);
	}

	my_fseek(kmc_pre, 4, SEEK_SET);
	if (fread(LUTS, sizeof(uint64), lut_recs, kmc_pre) != lut_recs)
	{
		std::cerr << "Some error occured while reading LUTS from kmc2 prefix file \n";
		exit(1);
	}
	fclose(kmc_pre);

	std::string kmc_suf_file_name = desc.file_src + ".kmc_suf";
	kmc_suf = fopen(kmc_suf_file_name.c_str(), "rb");

	if (!kmc_suf)
	{
		std::cerr << "Error: cannot open kmc2 suffix file\n";
		exit(1);
	}
	setvbuf(kmc_suf, NULL, _IONBF, 0);

	bins.reserve(header.no_of_bins);
	for (uint32 i = 0; i < header.no_of_bins; ++i)
		bins.emplace_back(i, LUTS + i * lut_size, *this);

	//starting threads

	bin_provider.init(bins);

	suf_bin_reader = new CSufBinReader<SIZE>(bin_provider, kmc_suf);
	suf_bin_reader_th = std::thread(std::ref(*suf_bin_reader));

	uint32 n_threads = desc.threads;

	if (n_threads < 3)
	{
		n_threads = n_child_threads = 1;		
		output_queue = new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY);
		childs.push_back(new CMergerChild<SIZE>(bins.begin(), bins.end(), *output_queue));
		childs_threads.push_back(std::thread(std::ref(*childs.front())));
		return;
	}

	else if (n_threads == 3)
	{
		n_child_threads = 2;
		n_parent_threads = 1; 
	}
	//based on experiment on 24 core machine
	else if (n_threads < 6)
	{
		n_child_threads = 3;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 9)
	{
		n_child_threads = 4;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 11)
	{
		n_child_threads = 5;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 14)
	{
		n_child_threads = 6;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 17)
	{
		n_child_threads = 7;
		n_parent_threads = n_threads - n_child_threads;
	}
	else
	{
		n_child_threads = (n_threads - 17) / 5 + 8;
		n_parent_threads = n_threads - n_child_threads;
	}
	childs_parent_queues.reserve(n_child_threads);
	childs.reserve(n_child_threads);

	uint32 bundle_size = KMC2_DB_READER_BUNDLE_CAPACITY;
	if (n_parent_threads < 2)
	{
		bundle_size = BUNDLE_CAPACITY;
	}

	for (uint32 i = 0; i < n_child_threads; ++i)
		childs_parent_queues.push_back(new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY, bundle_size));

	uint32 bins_per_thread = header.no_of_bins / n_child_threads;

	for (uint32 i = 0; i < n_child_threads - 1; ++i)
	{
		childs.push_back(new CMergerChild<SIZE>(bins.begin() + i * bins_per_thread, bins.begin() + (i + 1) * bins_per_thread, *childs_parent_queues[i]));
		childs_threads.push_back(std::thread(std::ref(*childs.back())));
	}

	//last one
	childs.push_back(new CMergerChild<SIZE>(bins.begin() + (n_child_threads - 1) * bins_per_thread, bins.end(), *childs_parent_queues.back()));
	childs_threads.push_back(std::thread(std::ref(*childs.back())));

	output_queue = new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY, bundle_size);

	parent = new CMergerParent<SIZE>(childs_parent_queues, *output_queue, n_parent_threads);
	parent_thread = std::thread(std::ref(*parent));
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderSorted<SIZE>::NextBundle(CBundle<SIZE>& bundle, bool& finished)
{
	if (output_queue->pop(bundle.Data()))
	{
		return;
	}

	for (auto& child_thread : childs_threads)
		child_thread.join();
	for (auto& child : childs)
		delete child;

	if(parent_thread.joinable()) //for small number of threads there is only one child thread and parent threads is not needed
		parent_thread.join();
	if (parent) //as above
		delete parent;

	delete output_queue;

	for (auto& q : childs_parent_queues)
		delete q;

	suf_bin_reader_th.join();
	delete suf_bin_reader;
	finished = true;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderSorted<SIZE>::IgnoreRest()
{
	output_queue->force_finish();

	for (auto& q : childs_parent_queues)
		q->force_finish();

	for (auto& child_thread : childs_threads)
		child_thread.join();
	for (auto& child : childs)
		delete child;

	if (parent_thread.joinable()) //for small number of threads there is only one child thread and parent threads is not needed
		parent_thread.join();

	if(parent) //as above
		delete parent;

	delete output_queue;

	for (auto& q : childs_parent_queues)
		delete q;

	bin_provider.force_to_finish();

	suf_bin_reader_th.join();
	delete suf_bin_reader;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC2DbReaderSorted<SIZE>::~CKMC2DbReaderSorted()
{
	delete[] LUTS;
	fclose(kmc_suf);
}

/*****************************************************************************************************************************/
/******************************************* CKCM2DbReaderSeqCounter_Base IMPLEMENTATION *************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE> CKCM2DbReaderSeqCounter_Base<SIZE>::CKCM2DbReaderSeqCounter_Base(const CKMC_header& header, const CInputDesc& desc) :
header(header),
desc(desc)
{
	suffix_bytes = (header.kmer_len - header.lut_prefix_len) / 4;
	record_size = suffix_bytes + header.counter_size;
	suffix_buff_size = SUFFIX_BUFF_BYTES / record_size * record_size;

	suffix_left_to_read = header.total_kmers * record_size;

	if (suffix_left_to_read < suffix_buff_size)
		suffix_buff_size = suffix_left_to_read;

	prefix_bytes = (header.lut_prefix_len + 3) / 4;

	kmer_bytes = prefix_bytes + suffix_bytes;

}

/*****************************************************************************************************************************/
/********************************************************* PROTECTED**********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> CKCM2DbReaderSeqCounter_Base<SIZE>::~CKCM2DbReaderSeqCounter_Base()
{
	if (suffix_file)
		fclose(suffix_file);
	delete[] suffix_buff;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKCM2DbReaderSeqCounter_Base<SIZE>::open_files()
{
	suffix_file_name = desc.file_src + ".kmc_suf";

	suffix_file = fopen(suffix_file_name.c_str(), "rb");


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
}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKCM2DbReaderSeqCounter_Base<SIZE>::reload_suf_buff()
{
	uint64 to_read = MIN(suffix_left_to_read, suffix_buff_size);
	if (to_read == 0)
		return false;
	uint64 readed = fread(suffix_buff, 1, to_read, suffix_file);
	if (readed != to_read)
	{
		std::cerr << "Error: some error while reading " << suffix_file_name << "\n";
		exit(1);
	}
	suffix_buff_pos = 0;
	suffix_left_to_read -= to_read;
	return true;
}


/*****************************************************************************************************************************/
/********************************************* CKMC2DbReaderSequential IMPLEMENTATION ****************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC2DbReaderSequential<SIZE>::CKMC2DbReaderSequential(const CKMC_header& header, const CInputDesc& desc) :
CKCM2DbReaderSeqCounter_Base<SIZE>(header, desc)
{
	this->open_files();

	prefix_file_name = desc.file_src + ".kmc_pre";
	prefix_file = fopen(prefix_file_name.c_str(), "rb");

	if (!prefix_file)
	{
		std::cerr << "Error: cannot open file: " << prefix_file_name << "\n";
		exit(1);
	}
	setvbuf(prefix_file, NULL, _IONBF, 0);
	my_fseek(prefix_file, 4 + sizeof(uint64), SEEK_SET);//skip KMCP and first value as it must be 0

	signle_bin_size = 1 << 2 * header.lut_prefix_len;
	map_size = (1 << 2 * header.signature_len) + 1;

	map_size_bytes = map_size * sizeof(uint32);

	no_of_bins = header.no_of_bins;

	this->prefix_buff_size = this->PREFIX_BUFF_BYTES / sizeof(uint64);

	this->suffix_left_to_read = this->header.total_kmers * this->record_size;

	if (this->suffix_left_to_read < this->suffix_buff_size)
		this->suffix_buff_size = this->suffix_left_to_read;

	this->prefix_left_to_read = (1 << this->header.lut_prefix_len * 2) * this->no_of_bins;

	if (this->prefix_left_to_read < this->prefix_buff_size)
		this->prefix_buff_size = this->prefix_left_to_read;

	prefix_mask = (1 << 2 * this->header.lut_prefix_len) - 1;

	allocate_buffers();

	my_fseek(prefix_file, 4 + sizeof(uint64), SEEK_SET);
	reload_pref_buff();

	this->reload_suf_buff();
	current_prefix_index = 0;
	this->suffix_number = 0;
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC2DbReaderSequential<SIZE>::NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter)
{
	while (true)
	{
		if (this->suffix_number >= this->header.total_kmers)
			return false;

		while (this->prefix_buff[this->prefix_buff_pos] <= this->suffix_number)
		{
			++current_prefix_index;
			++this->prefix_buff_pos;
			if (this->prefix_buff_pos >= this->prefix_buff_size)
				this->reload_pref_buff();
		}

		uchar* record = this->suffix_buff + this->suffix_buff_pos;
		uint32 pos = this->kmer_bytes - 1;

		uint32 current_prefix = static_cast<uint32>(current_prefix_index & prefix_mask);

		kmer.load(record, this->suffix_bytes);
		for (int32 i = this->prefix_bytes - 1; i >= 0; --i)
			kmer.set_byte(pos--, current_prefix >> (i << 3));

		counter = 0;
		for (int32 i = this->header.counter_size - 1; i >= 0; --i)
		{
			counter <<= 8;
			counter += record[i];
		}

		++this->suffix_number;
		this->suffix_buff_pos += this->record_size;

		if (this->suffix_buff_pos >= this->suffix_buff_size)
			this->reload_suf_buff();

		if (counter >= this->desc.cutoff_min && counter <= this->desc.cutoff_max)
			return true;
	}
}

/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC2DbReaderSequential<SIZE>::~CKMC2DbReaderSequential()
{
	if (prefix_file)
		fclose(prefix_file);
	delete[] prefix_buff;
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderSequential<SIZE>::allocate_buffers()
{
	this->suffix_buff = new uchar[this->suffix_buff_size + sizeof(uint64)]; //+ sizeof(uint64) because CKmer::load_fast may look out of buffer
	this->prefix_buff = new uint64[this->prefix_buff_size];
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderSequential<SIZE>::reload_pref_buff()
{
	uint64 to_read = MIN(this->prefix_left_to_read, this->prefix_buff_size);
	this->prefix_buff_pos = 0;
	if (fread(prefix_buff, sizeof(uint64), to_read, prefix_file) != to_read)
	{
		std::cerr << "Error: some error while reading " << prefix_file_name << "\n";
		exit(1);
	}
	this->prefix_left_to_read -= to_read;
	if (to_read < this->prefix_buff_size)
	{
		this->prefix_buff[to_read] = this->header.total_kmers;//guard
	}
}


/*****************************************************************************************************************************/
/******************************************** CKMC2DbReaderCountersOnly IMPLEMENTATION ***************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC2DbReaderCountersOnly<SIZE>::CKMC2DbReaderCountersOnly(const CKMC_header& header, const CInputDesc& desc) :
CKCM2DbReaderSeqCounter_Base<SIZE>(header, desc)
{
	this->open_files();
	allocate_buffers();
	this->reload_suf_buff();
	this->suffix_number = 0;
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC2DbReaderCountersOnly<SIZE>::NextCounter(uint32& counter)
{
	while (true)
	{
		if (this->suffix_number >= this->header.total_kmers)
			return false;

		uchar* record = this->suffix_buff + this->suffix_buff_pos + this->suffix_bytes;

		counter = 0;
		for (int32 i = this->header.counter_size - 1; i >= 0; --i)
		{
			counter <<= 8;
			counter += record[i];
		}

		++this->suffix_number;
		this->suffix_buff_pos += this->record_size;

		if (this->suffix_buff_pos >= this->suffix_buff_size)
			this->reload_suf_buff();

		if (counter >= this->desc.cutoff_min && counter <= this->desc.cutoff_max)
			return true;
	}
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderCountersOnly<SIZE>::allocate_buffers()
{
	this->suffix_buff = new uchar[this->suffix_buff_size];
}



/*****************************************************************************************************************************/
/************************************************* CKMC2DbReader IMPLEMENTATION **********************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC2DbReader<SIZE>::CKMC2DbReader(const CKMC_header& header, const CInputDesc& desc, CPercentProgress& percent_progress, KMCDBOpenMode open_mode) :
	percent_progress(percent_progress)
{
	progress_id = percent_progress.RegisterItem(header.total_kmers);
	switch (open_mode)
	{
	case KMCDBOpenMode::sorted:
		db_reader_sorted = std::make_unique <CKMC2DbReaderSorted<SIZE>>(header, desc);
		break;
	case KMCDBOpenMode::sequential:
		db_reader_sequential = std::make_unique<CKMC2DbReaderSequential<SIZE>>(header, desc);
		break;
	case KMCDBOpenMode::counters_only:
		db_reader_counters_only = std::make_unique<CKMC2DbReaderCountersOnly<SIZE>>(header, desc);
		break;
	default: //should never be here
		std::cerr << "Error: unknow open mode \n";
		exit(1);
	}
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReader<SIZE>::NextBundle(CBundle<SIZE>& bundle)
{
	db_reader_sorted->NextBundle(bundle, this->finished);
	percent_progress.UpdateItem(progress_id, bundle.Size());
	if (this->finished)
	{
		percent_progress.Complete(progress_id);
	}
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReader<SIZE>::IgnoreRest()
{
	if (this->finished)
		return;
	db_reader_sorted->IgnoreRest();
	this->finished = true;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC2DbReader<SIZE>::NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter)
{
	if (db_reader_sequential->NextKmerSequential(kmer, counter))
	{
		percent_progress.UpdateItem(progress_id);
		return true;
	}
	percent_progress.Complete(progress_id);
	return false;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC2DbReader<SIZE>::NextCounter(uint32& counter)
{
	if (db_reader_counters_only->NextCounter(counter))
	{
		percent_progress.UpdateItem(progress_id);
		return true;
	}
	percent_progress.Complete(progress_id);
	return false;
}



#endif

