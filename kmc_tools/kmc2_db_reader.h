/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
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

#include <condition_variable>


//Forward declaration
template<unsigned SIZE> class CKMC2DbReaderSorted;

template<unsigned SIZE> class CBin;


struct CBinBuff //must be moveable
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
	using to_read_t = std::tuple<uint32, uint64, uchar*, uint32>;//bin_id, file_pos, bufer to read, size to read

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
#ifdef ENABLE_LOGGER
		CTimer timer;
#endif

		while (bin_provider.get_next_to_read(bin_id, file_pos, buf, size))
		{		
			my_fseek(suf_file, file_pos, SEEK_SET);
#ifdef ENABLE_LOGGER		
			timer.start();
#endif						
			if (fread(buf, 1, size, suf_file) != size)			
			{
				std::cout << "Error while reading sufix file\n";
				exit(1);
			}
#ifdef ENABLE_LOGGER
			CLoger::GetLogger().log_operation("fread", this, timer.get_time());
			timer.start();
#endif
			bin_provider.notify_bin_filled(bin_id);
		}
	}

};


template<unsigned SIZE> class CBin
{
public:
	CBin(uint32 bin_id, uint64* LUT, CKMC2DbReaderSorted<SIZE>& kmc2_db);
	bool NextKmer(CKmer<SIZE>& kmer, uint32& counter);

	uint64 get_kmer_number()
	{
		return kmer_number;
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
	//Because VS2013 does generate default move ctro here
	CBin(CBin&& o) throw():
		bin_id(o.bin_id),
		bin_buff(std::move(o.bin_buff)),
		LUT(o.LUT),
		pos(o.pos),
		bin_provider(o.bin_provider),
		kmc2_db(o.kmc2_db),
	
		kmer_number(o.kmer_number), kmer_number_end(o.kmer_number_end),
		kmer_bytes(o.kmer_bytes), prefix_bytes(o.prefix_bytes), suffix_bytes(o.suffix_bytes), counter_size(o.counter_size), record_size(o.record_size),
		prefix(o.prefix),
		max_prefix(o.max_prefix)
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
	uint64 kmer_number, kmer_number_end;
	uint32 kmer_bytes, prefix_bytes, suffix_bytes, counter_size, record_size;
	uint64 prefix = 0;
	uint64 max_prefix;
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
		start = b.get_kmer_number();
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
	inline bool get_min(CBundleData<SIZE>& bundle_data);

	inline void reset();

private:
	inline void update_heap();

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
	CMergerParent(std::vector<CCircularQueue<SIZE>*>& input_queues, CCircularQueue<SIZE>& output_queue);

	void operator()();

private:
	std::vector<CBundleData<SIZE>> input_bundles;
	std::vector<CCircularQueue<SIZE>*>& input_queues;
	CBundleData<SIZE> output_bundle;
	CCircularQueue<SIZE>& output_queue;
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

	CMergerParent<SIZE>* parent;
	std::thread parent_thread;

	CCircularQueue<SIZE> output_queue;
	std::vector<CCircularQueue<SIZE>*> childs_parent_queues;

	std::vector<CMergerChild<SIZE>*> childs;
	std::vector<std::thread> childs_threads;

	//mutable std::mutex mtx;
};

//************************************************************************************************************
// CKCM2DbReaderSeqCounter_Base - Base class for classes to access k-mers one by one (not sorted) or 
// for counters only from KMC2 database
//************************************************************************************************************
template <unsigned SIZE> class CKCM2DbReaderSeqCounter_Base
{
protected:
	CKCM2DbReaderSeqCounter_Base(const CKMC_header& header, const CInputDesc& desc);
	~CKCM2DbReaderSeqCounter_Base();

	void open_files();
	bool reload_suf_buff();

	static const uint32 PREFIX_BUFF_BYTES = KMC2_DB_READER_PREFIX_BUFF_BYTES;
	static const uint32 SUFIX_BUFF_BYTES = KMC2_DB_READER_SUFIX_BUFF_BYTES;

	const CKMC_header& header;
	const CInputDesc& desc;

	uint32 sufix_bytes;
	uint32 record_size; //of sufix, in bytes
	uint64 sufix_buff_size, sufix_buff_pos, sufix_left_to_read;
	uint64 prefix_buff_size, prefix_buff_pos, prefix_left_to_read;
	uint64 sufix_number;

	uint32 kmer_bytes, prefix_bytes;

	uchar* sufix_buff = nullptr;

	FILE* sufix_file;
	std::string sufix_file_name;
};

//************************************************************************************************************
// CKMC2DbReaderSequential - Produce k-mers sequentialy from KMC2 database (they are not sorted!)
//************************************************************************************************************
template<unsigned SIZE> class CKMC2DbReaderSequential : public CKCM2DbReaderSeqCounter_Base<SIZE>
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
};

//************************************************************************************************************
// CKMC2DbReaderCountersOnly - Produce counters of k-mers from KMC2 database
//************************************************************************************************************
template<unsigned SIZE> class CKMC2DbReaderCountersOnly : CKCM2DbReaderSeqCounter_Base<SIZE>
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
template<unsigned SIZE> class CKMC2DbReader : public CInput<SIZE>
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
	kmer_number = LUT[0];
	kmer_number_end = LUT[kmc2_db.lut_size];
	prefix_bytes = (kmc2_db.header.lut_prefix_len + 3) / 4;
	kmer_bytes = prefix_bytes + suffix_bytes;

	record_size = suffix_bytes + counter_size;
}


/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/


/*****************************************************************************************************************************/
template<unsigned SIZE> bool CBin<SIZE>::NextKmer(CKmer<SIZE>& kmer, uint32& counter)
{
	while (true)
	{
		if (kmer_number >= kmer_number_end)
			return false;

		if (pos >= bin_buff.size)
			reload_suf_buf();

		//skip empty
		while (LUT[prefix + 1] <= kmer_number)
		{
			++prefix;
		}

		uint32 in_kmer_pos = kmer_bytes - 1;
		uchar* record = bin_buff.buf + pos;
		kmer.load(record, suffix_bytes);
		for (int32 i = prefix_bytes - 1; i >= 0; --i)
			kmer.set_byte(in_kmer_pos--, uchar(prefix >> (i << 3)));

		counter = 0;
		for (int32 i = counter_size - 1; i >= 0; --i)
		{
			counter <<= 8;
			counter += record[i];
		}

		++kmer_number;
		pos += record_size;

		if (counter >= kmc2_db.desc.cutoff_min && counter <= kmc2_db.desc.cutoff_max)
			return true;
	}
	return true;
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
template<unsigned SIZE> inline bool CKmerPQ<SIZE>::get_min(CBundleData<SIZE>& bundle_data)
{
	if (pos <= 1)
		return false;
	bundle_data.Insert(elems[1].first, descs[elems[1].second].second);

	update_heap();
	return true;
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
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> inline void CKmerPQ<SIZE>::update_heap()
{
	uint32 desc_id = elems[1].second;
	CBin<SIZE>* bin = descs[desc_id].first;
	CKmer<SIZE> kmer;
	uint32 counter;
	if (!bin->NextKmer(kmer, counter))
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
}



/*****************************************************************************************************************************/
/*********************************************** CMergerParent IMPLEMENTATION ************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE>	CMergerParent<SIZE>::CMergerParent(std::vector<CCircularQueue<SIZE>*>& input_queues, CCircularQueue<SIZE>& output_queue) :
	input_queues(input_queues),
	output_queue(output_queue)
{
	input_bundles.resize(input_queues.size());
}


/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CMergerParent<SIZE>::operator()()
{
	//init
	//for (uint32 i = 0; i < input_queues.size(); ++i)
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

	//run
	uint32 index_of_min = 0;
	while (input_bundles.size())
	{
		index_of_min = 0;
		for (uint32 i = 1; i < input_bundles.size(); ++i)
		{
			if (input_bundles[i].TopKmer() < input_bundles[index_of_min].TopKmer())
				index_of_min = i;
		}

		output_bundle.Insert(input_bundles[index_of_min].TopKmer(), input_bundles[index_of_min].TopCounter());
		input_bundles[index_of_min].Pop();
		if (input_bundles[index_of_min].Empty())
		{
			if (!input_queues[index_of_min]->pop(input_bundles[index_of_min]))
			{
				input_queues.erase(input_queues.begin() + index_of_min);
				input_bundles.erase(input_bundles.begin() + index_of_min);
			}
		}


		if (output_bundle.Full())
		{
			if (!output_queue.push(output_bundle))
				break;
		}
	}
	if (!output_bundle.Empty())
		output_queue.push(output_bundle);
	output_queue.mark_completed();
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

	CBundleData<SIZE> bundle_data;
	while (kmers_pq.get_min(bundle_data))
	{
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
/********************************************* CKMC2DbReaderSorted IMPLEMENTATION ********************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
template<unsigned SIZE> CKMC2DbReaderSorted<SIZE>::CKMC2DbReaderSorted(const CKMC_header& header, const CInputDesc& desc) :
	header(header),
	desc(desc),
	output_queue(DEFAULT_CIRCULAL_QUEUE_CAPACITY)
{
	LUTS = nullptr;
	lut_size = 1 << 2 * header.lut_prefix_len;
	uint32 lut_recs = (1 << 2 * header.lut_prefix_len) * header.no_of_bins + 1;
	LUTS = new uint64[lut_recs];
	suffix_bytes = (header.kmer_len - header.lut_prefix_len) / 4;
	record_size = suffix_bytes + header.counter_size;
	if (!LUTS)
	{
		std::cout << "cannot allocate memory for LUTS of KMC2 database\n";
		exit(1);
	}

	std::string kmc_pre_file_name = desc.file_src + ".kmc_pre";
	FILE* kmc_pre = fopen(kmc_pre_file_name.c_str(), "rb");
	if (!kmc_pre)
	{
		std::cout << "Cannot open kmc2 prefix file to read LUTS";
		exit(1);
	}

	my_fseek(kmc_pre, 4, SEEK_SET);
	if (fread(LUTS, sizeof(uint64), lut_recs, kmc_pre) != lut_recs)
	{
		std::cout << "Some error occured while reading LUTS from kmc2 prefix file \n";
		exit(1);
	}
	fclose(kmc_pre);

	std::string kmc_suf_file_name = desc.file_src + ".kmc_suf";
	kmc_suf = fopen(kmc_suf_file_name.c_str(), "rb");

	if (!kmc_suf)
	{
		std::cout << "Cannot open kmc2 suffix file\n";
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

	n_child_threads = desc.threads;

	childs_parent_queues.reserve(n_child_threads);

	for (uint32 i = 0; i < n_child_threads; ++i)
		childs_parent_queues.push_back(new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY));


	uint32 bins_per_thread = header.no_of_bins / n_child_threads;

	for (uint32 i = 0; i < n_child_threads - 1; ++i)
	{
		childs.push_back(new CMergerChild<SIZE>(bins.begin() + i * bins_per_thread, bins.begin() + (i + 1) * bins_per_thread, *childs_parent_queues[i]));
		childs_threads.push_back(std::thread(std::ref(*childs.back())));
	}

	//last one
	childs.push_back(new CMergerChild<SIZE>(bins.begin() + (n_child_threads - 1) * bins_per_thread, bins.end(), *childs_parent_queues.back()));
	childs_threads.push_back(std::thread(std::ref(*childs.back())));

	parent = new CMergerParent<SIZE>(childs_parent_queues, output_queue);
	parent_thread = std::thread(std::ref(*parent));




}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderSorted<SIZE>::NextBundle(CBundle<SIZE>& bundle, bool& finished)
{
	if (output_queue.pop(bundle.Data()))
	{
		return;
	}

	for (auto& child_thread : childs_threads)
		child_thread.join();
	for (auto& child : childs)
		delete child;

	parent_thread.join();
	delete parent;

	for (auto& q : childs_parent_queues)
		delete q;

	suf_bin_reader_th.join();
	delete suf_bin_reader;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderSorted<SIZE>::IgnoreRest()
{
	output_queue.force_finish();

	for (auto& q : childs_parent_queues)
		q->force_finish();

	for (auto& child_thread : childs_threads)
		child_thread.join();
	for (auto& child : childs)
		delete child;

	parent_thread.join();
	delete parent;

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
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
//template<unsigned SIZE> void CKMC2DbReaderSorted<SIZE>::get_suf_buf_part(uchar* &buf, uint64 start, uint32 size)
//{
//#ifdef ENABLE_LOGGER
//	CTimer timer;
//	timer.start();
//#endif
//	std::unique_lock<std::mutex> lck(mtx);
//#ifdef ENABLE_LOGGER
//	CLoger::GetLogger().log_operation("waiting for lock", this, timer.get_time());
//	timer.start();
//#endif
//	start = 4 + start * record_size;
//	size *= record_size;
//
//
//	my_fseek(kmc_suf, start, SEEK_SET);
//	if (fread(buf, 1, size, kmc_suf) != size)
//	{
//		std::cout << "Error: some error occured while reading " << desc.file_src << ".kmc_suf file\n";
//		exit(1);
//	}
//#ifdef ENABLE_LOGGER
//	CLoger::GetLogger().log_operation("fread time", this, timer.get_time());
//#endif
//}


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
	sufix_bytes = (header.kmer_len - header.lut_prefix_len) / 4;
	record_size = sufix_bytes + header.counter_size;
	sufix_buff_size = SUFIX_BUFF_BYTES / record_size * record_size;

	sufix_left_to_read = header.total_kmers * record_size;

	if (sufix_left_to_read < sufix_buff_size)
		sufix_buff_size = sufix_left_to_read;

	prefix_bytes = (header.lut_prefix_len + 3) / 4;

	kmer_bytes = prefix_bytes + sufix_bytes;

}

/*****************************************************************************************************************************/
/********************************************************* PROTECTED**********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> CKCM2DbReaderSeqCounter_Base<SIZE>::~CKCM2DbReaderSeqCounter_Base()
{
	if (sufix_file)
		fclose(sufix_file);
	delete[] sufix_buff;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKCM2DbReaderSeqCounter_Base<SIZE>::open_files()
{
	sufix_file_name = desc.file_src + ".kmc_suf";

	sufix_file = fopen(sufix_file_name.c_str(), "rb");


	if (!sufix_file)
	{
		std::cout << "Error: cannot open file: " << sufix_file_name << "\n";
		exit(1);
	}
	setvbuf(sufix_file, NULL, _IONBF, 0);

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
}

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKCM2DbReaderSeqCounter_Base<SIZE>::reload_suf_buff()
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
			std::cout << "Error: cannot open file: " << prefix_file_name << "\n";
			exit(1);
		}
		setvbuf(prefix_file, NULL, _IONBF, 0);
		my_fseek(prefix_file, 4 + sizeof(uint64), SEEK_SET);//skip KMCP and first value as it must be 0

		signle_bin_size = 1 << 2 * header.lut_prefix_len;
		map_size = (1 << 2 * header.signature_len) + 1;

		map_size_bytes = map_size * sizeof(uint32);

		no_of_bins = header.no_of_bins;

		this->prefix_buff_size = this->PREFIX_BUFF_BYTES / sizeof(uint64);

		this->sufix_left_to_read = this->header.total_kmers * this->record_size;

		if (this->sufix_left_to_read < this->sufix_buff_size)
			this->sufix_buff_size = this->sufix_left_to_read;

		this->prefix_left_to_read = (1 << this->header.lut_prefix_len * 2) * this->no_of_bins;

		if (this->prefix_left_to_read < this->prefix_buff_size)
			this->prefix_buff_size = this->prefix_left_to_read;

		prefix_mask = (1 << 2 * this->header.lut_prefix_len) - 1;

		allocate_buffers();

		my_fseek(prefix_file, 4 + sizeof(uint64), SEEK_SET);
		reload_pref_buff();

		this->reload_suf_buff();
		current_prefix_index = 0;
		this->sufix_number = 0;
	}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC2DbReaderSequential<SIZE>::NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter)
{
	while (true)
	{
		if (this->sufix_number >= this->header.total_kmers)
			return false;

		while (this->prefix_buff[this->prefix_buff_pos] <= this->sufix_number)
		{
			++current_prefix_index;
			++this->prefix_buff_pos;
			if (this->prefix_buff_pos >= this->prefix_buff_size)
				this->reload_pref_buff();
		}

		uchar* record = this->sufix_buff + this->sufix_buff_pos;
		uint32 pos = this->kmer_bytes - 1;

		uint32 current_prefix = static_cast<uint32>(current_prefix_index & prefix_mask);

		kmer.load(record, this->sufix_bytes);
		for (int32 i = this->prefix_bytes - 1; i >= 0; --i)
			kmer.set_byte(pos--, current_prefix >> (i << 3));

		counter = 0;
		for (int32 i = this->header.counter_size - 1; i >= 0; --i)
		{
			counter <<= 8;
			counter += record[i];
		}

		++this->sufix_number;
		this->sufix_buff_pos += this->record_size;

		if (this->sufix_buff_pos >= this->sufix_buff_size)
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
	this->sufix_buff = new uchar[this->sufix_buff_size];
	this->prefix_buff = new uint64[this->prefix_buff_size];
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReaderSequential<SIZE>::reload_pref_buff()
{
	uint64 to_read = MIN(this->prefix_left_to_read, this->prefix_buff_size);
	this->prefix_buff_pos = 0;
	if (fread(prefix_buff, sizeof(uint64), to_read, prefix_file) != to_read)
	{
		std::cout << "Error: some error while reading " << prefix_file_name << "\n";
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
		this->sufix_number = 0;
	}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> bool CKMC2DbReaderCountersOnly<SIZE>::NextCounter(uint32& counter)
{
	while (true)
	{
		if (this->sufix_number >= this->header.total_kmers)
			return false;

		uchar* record = this->sufix_buff + this->sufix_buff_pos + this->sufix_bytes;

		counter = 0;
		for (int32 i = this->header.counter_size - 1; i >= 0; --i)
		{
			counter <<= 8;
			counter += record[i];
		}

		++this->sufix_number;
		this->sufix_buff_pos += this->record_size;

		if (this->sufix_buff_pos >= this->sufix_buff_size)
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
	this->sufix_buff = new uchar[this->sufix_buff_size];
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
		std::cout << "Error: unknow open mode \n";
		exit(1);
	}
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReader<SIZE>::NextBundle(CBundle<SIZE>& bundle)
{
#ifdef ENABLE_LOGGER
	CTimer timer;
	timer.start();
#endif
	db_reader_sorted->NextBundle(bundle, this->finished);
	percent_progress.UpdateItem(progress_id, bundle.Size());
	if(this->finished)
	{
		percent_progress.Complete(progress_id);
	}
#ifdef ENABLE_LOGGER
	CLoger::GetLogger().log_operation("pobranie bundla z wejscia", this, timer.get_time());
#endif
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void CKMC2DbReader<SIZE>::IgnoreRest()
{
	db_reader_sorted->IgnoreRest();
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