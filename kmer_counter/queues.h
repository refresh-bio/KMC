/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _QUEUES_H
#define _QUEUES_H

#include "defs.h"
#include <stdio.h>
#include <iostream>
#include <tuple>
#include <queue>
#include <list>
#include <set>
#include <map>
#include <string>
#include "mem_disk_file.h"

using namespace std;

#include <thread>
#include <mutex>
#include <condition_variable>

using std::thread;

//************************************************************************************************************

enum class FilePart { Begin, Middle, End };
enum class CompressionType { plain, gzip, bzip2 };


// Reader will clasify reads as normal or long. Distinction is not strict, it depends on current configuration (mem limit and no. of threads). In case of fastq reader will assure, that
// only header and read are sent to splitter (qual and its header will be ignored) 
enum class ReadType { normal_read, long_read, na }; //there is no strict distincion between normal and long reads, read will be clasified as long if reader is not able to determine whole fastq/fasta record (header, read, header, qual/header, read)


class CBinaryPackQueue
{
	std::queue<tuple<uchar*, uint64, FilePart, CompressionType>> q;
	std::mutex mtx;
	std::condition_variable cv_pop, cv_push;
	bool completed = false;
	bool stop = false;
public:
	bool push(uchar* data, uint64 size, FilePart file_part, CompressionType mode)
	{
		std::lock_guard<std::mutex> lck(mtx);
		if (stop)
			return false;
		q.emplace(data, size, file_part, mode);
		if (q.size() == 1) //was empty
			cv_pop.notify_all();
		return true;
	}

	void mark_completed()
	{
		std::lock_guard<std::mutex> lck(mtx);
		completed = true;
		cv_pop.notify_all();
	}
	bool pop(uchar* &data, uint64 &size, FilePart &file_part, CompressionType &mode)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this]{return !q.empty() || completed; });
		if (q.empty())
			return false;

		data = get<0>(q.front());
		size = get<1>(q.front());
		file_part = get<2>(q.front());
		mode = get<3>(q.front());

		q.pop();
		return true;
	}

	bool peek_next_pack(uchar* &data, uint64 &size)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] {return !q.empty() || completed; });
		if (q.empty())
			return false;
		data = get<0>(q.front());
		size = get<1>(q.front());
		return get<2>(q.front()) != FilePart::End;
	}

	bool is_next_last()
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this]{return !q.empty() || completed; });
		if (q.empty())
			return true;
		return get<2>(q.front()) == FilePart::End;
	}

	void ignore_rest()
	{
		lock_guard<mutex> lck(mtx);
		stop = true;
		cv_push.notify_all();
	}
};


class CMissingEOL_at_EOF_counter
{
	uint32 n_missing;
	std::mutex mtx;
public:
	CMissingEOL_at_EOF_counter()
	{
		Reset();
	}
	void Reset()
	{
		n_missing = 0;
	}
	void RegisterMissingEOL()
	{
		std::lock_guard<std::mutex> lck(mtx);
		++n_missing;
	}
	uint32 Get()
	{
		return n_missing;
	}
	
};

//************************************************************************************************************
class CInputFilesQueue {
	typedef string elem_t;
	typedef queue<elem_t, list<elem_t>> queue_t;

	queue_t q;
	bool is_completed;

	mutable mutex mtx;								// The mutex to synchronise on

public:
	queue_t GetCopy()
	{
		return q;
	}
	CInputFilesQueue(const vector<string> &file_names) {
		unique_lock<mutex> lck(mtx);

		for(vector<string>::const_iterator p = file_names.begin(); p != file_names.end(); ++p)
			q.push(*p);

		is_completed = false;
	};
	~CInputFilesQueue() {};

	bool empty() {
		lock_guard<mutex> lck(mtx);
		return q.empty();
	}
	bool completed() {
		lock_guard<mutex> lck(mtx);
		return q.empty() && is_completed;
	}
	void mark_completed() {
		lock_guard<mutex> lck(mtx);
		is_completed = true;
	}
	bool pop(string &file_name) {
		lock_guard<mutex> lck(mtx);

		if(q.empty())
			return false;

		file_name = q.front();
		q.pop();

		return true;
	}
};

//************************************************************************************************************
class CPartQueue {
	typedef tuple<uchar *, uint64, ReadType> elem_t;
	typedef queue<elem_t, list<elem_t>> queue_t;

	queue_t q;
	bool is_completed;
	int n_readers;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	CPartQueue(int _n_readers) {
		unique_lock<mutex> lck(mtx);
		is_completed    = false;
		n_readers       = _n_readers;
	};
	~CPartQueue() {};

	bool empty() {
		lock_guard<mutex> lck(mtx);
		return q.empty();
	}

	bool completed() {
		lock_guard<mutex> lck(mtx);
		return q.empty() && !n_readers;
	}

	void mark_completed() {
		lock_guard<mutex> lck(mtx);
		n_readers--;
		if(!n_readers)
			cv_queue_empty.notify_all();
	}

	void push(uchar *part, uint64 size, ReadType read_type) {
		unique_lock<mutex> lck(mtx);
		
		bool was_empty = q.empty();
		q.push(make_tuple(part, size, read_type));

		if(was_empty)
			cv_queue_empty.notify_all();
	}

	bool pop(uchar *&part, uint64 &size, ReadType& read_type) {
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this]{return !this->q.empty() || !this->n_readers;}); 

		if(q.empty())
			return false;
		std::tie(part, size, read_type) = q.front();
		q.pop();

		return true;
	}
};

//************************************************************************************************************
class CStatsPartQueue
{
	typedef tuple<uchar *, uint64, ReadType> elem_t;
	typedef queue<elem_t, list<elem_t>> queue_t;

	queue_t q;

	mutable mutex mtx;
	condition_variable cv_queue_empty;
	int n_readers;
	int64 bytes_to_read;
public:
	CStatsPartQueue(int _n_readers, int64 _bytes_to_read)
	{
		unique_lock<mutex> lck(mtx);
		n_readers = _n_readers;
		bytes_to_read = _bytes_to_read;
	}

	~CStatsPartQueue() {};

	void mark_completed() {
		lock_guard<mutex> lck(mtx);
		n_readers--;
		if (!n_readers)
			cv_queue_empty.notify_all();
	}

	bool completed() {
		lock_guard<mutex> lck(mtx);
		return q.empty() && !n_readers;
	}

	bool push(uchar *part, uint64 size, ReadType read_type) {
		unique_lock<mutex> lck(mtx);

		if (bytes_to_read <= 0)
			return false;

		bool was_empty = q.empty();
		q.push(make_tuple(part, size, read_type));
		bytes_to_read -= size;
		if (was_empty)
			cv_queue_empty.notify_one();

		return true;
	}

	bool pop(uchar *&part, uint64 &size, ReadType& read_type) {
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this]{return !this->q.empty() || !this->n_readers; });

		if (q.empty())
			return false;

		std::tie(part, size, read_type) = q.front();
		q.pop();

		return true;
	}


};

//************************************************************************************************************
class CBinPartQueue {
	typedef tuple<int32, uchar *, uint32, uint32, list<pair<uint64, uint64>>> elem_t;
	typedef queue<elem_t, list<elem_t>> queue_t;
	queue_t q;

	int n_writers;
	bool is_completed;

	mutable mutex mtx;						// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	CBinPartQueue(int _n_writers) {
		lock_guard<mutex> lck(mtx);

		n_writers       = _n_writers;
		is_completed    = false;
	}
	~CBinPartQueue() {}

	bool empty() {
		lock_guard<mutex> lck(mtx);
		return q.empty();
	}
	bool completed() {
		lock_guard<mutex> lck(mtx);
		return q.empty() && !n_writers;
	}
	void mark_completed() {
		lock_guard<mutex> lck(mtx);
		n_writers--;
		if(!n_writers)
			cv_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *part, uint32 true_size, uint32 alloc_size, list<pair<uint64, uint64>>& expander_parts) {
		unique_lock<mutex> lck(mtx);

		bool was_empty = q.empty();
		q.push(std::make_tuple(bin_id, part, true_size, alloc_size, std::move(expander_parts)));
		
		if(was_empty)
			cv_queue_empty.notify_all();
	}
	bool pop(int32 &bin_id, uchar *&part, uint32 &true_size, uint32 &alloc_size, list<pair<uint64, uint64>>& expander_parts) {
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this]{return !q.empty() || !n_writers;}); 

		if(q.empty())
			return false;

		bin_id     = get<0>(q.front());
		part       = get<1>(q.front());
		true_size  = get<2>(q.front());
		alloc_size = get<3>(q.front());
		expander_parts = std::move(get<4>(q.front()));
		q.pop();

		return true;
	}
};

//************************************************************************************************************
class CExpanderPackDesc
{
	vector<list<pair<uint64, uint64>>> data; //end_pos, n_plus_x_recs
	std::mutex mtx;
public:
	CExpanderPackDesc(uint32 n_bins)
	{
		data.resize(n_bins);
	}
	void push(uint32 bin_id, list<pair<uint64, uint64>>& l)
	{
		lock_guard<mutex> lck(mtx);
		data[bin_id].splice(data[bin_id].end(), l);
	}
	void pop(uint32 bin_id, list<pair<uint64, uint64>>& l)
	{
		lock_guard<mutex> lck(mtx);
		l = std::move(data[bin_id]);
		data[bin_id].clear();
	}
};
class CBinDesc {
	typedef tuple<string, int64, uint64, uint32, uint32, CMemDiskFile*, uint64, uint64> desc_t;
	typedef map<int32, desc_t> map_t;

	map_t m;
	int32 bin_id;

	vector<int32> random_bins;

	vector<int32> sorted_bins;

	mutable mutex mtx;

public:
	CBinDesc() {
		lock_guard<mutex> lck(mtx);
		bin_id = -1;
	}
	~CBinDesc() {}

	void reset_reading() {
		lock_guard<mutex> lck(mtx);
		bin_id = -1;
	}

	bool empty() {
		lock_guard<mutex> lck(mtx);
		return m.empty();
	}

	void init_sort(vector<pair<int32, int64>>& bin_sizes)
	{
		//for (auto& p : m)
		//	bin_sizes.push_back(make_pair(p.first, get<2>(p.second)));

		//sort(bin_sizes.begin(), bin_sizes.end(), [](const pair<int32, int64>& l, const pair<int32, int64>& r){
		//	return l.second > r.second;
		//});

		sorted_bins.reserve(bin_sizes.size());

		for (auto& e : bin_sizes)
			sorted_bins.push_back(e.first);

	}

	uint64 get_req_size(uint32 sorting_phases, int64 file_size, int64 kxmers_size, int64 out_buffer_size, int64 kxmer_counter_size, int64 lut_size)
	{
		int64 part1_size;
		int64 part2_size;

		if (sorting_phases % 2 == 0)
		{
			part1_size = kxmers_size + kxmer_counter_size;
			part2_size = max(max(file_size, kxmers_size), out_buffer_size + lut_size);
		}
		else
		{
			part1_size = max(kxmers_size + kxmer_counter_size, file_size);
			part2_size = max(kxmers_size, out_buffer_size + lut_size);
		}
		return part1_size + part2_size;
	}

	int64 round_up_to_alignment(int64 x)
	{
		return (x + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;
	}

	uint64 get_n_rec_sum()
	{
		uint64 res = 0;
		for (auto& e : m)
			res += get<2>(e.second);		
		return res;
	}

	vector<pair<int32, int64>> get_sorted_req_sizes(uint32 max_x, const uint64 size_of_kmer_t, uint32 cutoff_min, int64 cutoff_max, int64 counter_max, uint32 lut_prefix_len)
	{
		lock_guard<mutex> lck(mtx);
		vector<pair<int32, int64>> bin_sizes;

		uint64 size;
		uint64 n_rec;
		uint64 n_plus_x_recs;		
		uint32 kmer_len;

		for (auto& p : m)
		{
			int32 bin_id = p.first;

			size = (uint64)get<1>(m[bin_id]);
			n_rec = get<2>(m[bin_id]);			
			kmer_len = get<4>(m[bin_id]);
			n_plus_x_recs = get<6>(m[bin_id]);

			uint64 input_kmer_size;
			uint64 kxmer_counter_size;
			uint32 kxmer_symbols;
			if (max_x)
			{
				input_kmer_size = n_plus_x_recs * size_of_kmer_t;
				kxmer_counter_size = n_plus_x_recs * sizeof(uint32);
				kxmer_symbols = kmer_len + max_x + 1;
			}
			else
			{
				input_kmer_size = n_rec * size_of_kmer_t;
				kxmer_counter_size = 0;
				kxmer_symbols = kmer_len;
			}
			
			uint64 max_out_recs = (n_rec + 1) / max(cutoff_min, 1u);

			uint64 counter_size = MIN(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));

			uint32 kmer_symbols = kmer_len - lut_prefix_len;
			uint64 kmer_bytes = kmer_symbols / 4;
			uint64 out_buffer_size = max_out_recs * (kmer_bytes + counter_size);

			uint32 rec_len = (kxmer_symbols + 3) / 4;

			uint64 lut_recs = 1ull << (2 * lut_prefix_len);
			uint64 lut_size = lut_recs * sizeof(uint64);

			auto req_size = get_req_size(rec_len, round_up_to_alignment(size), round_up_to_alignment(input_kmer_size), round_up_to_alignment(out_buffer_size), round_up_to_alignment(kxmer_counter_size), round_up_to_alignment(lut_size));
			bin_sizes.push_back(make_pair(p.first, req_size));
		}

		sort(bin_sizes.begin(), bin_sizes.end(), [](const pair<int32, int64>& l, const pair<int32, int64>& r){
			return l.second > r.second;
		});
		return bin_sizes;
	}


	int32 get_next_sort_bin()
	{
		lock_guard<mutex> lck(mtx);
		if (bin_id == -1)
			bin_id = 0;
		else
			++bin_id;
		if (bin_id >= (int32)m.size())
			return -1000;
		return sorted_bins[bin_id];
	}
	void init_random()
	{
		lock_guard<mutex> lck(mtx);
		vector<pair<int32, int64>> bin_sizes;

		for (auto& p : m)
			bin_sizes.push_back(make_pair(p.first, get<2>(p.second)));

		sort(bin_sizes.begin(), bin_sizes.end(), [](const pair<int32, int64>& l, const pair<int32, int64>& r){
			return l.second > r.second;
		});

		uint32 no_sort_start = uint32(0.6 * bin_sizes.size());
		uint32 no_sort_end = uint32(0.8 * bin_sizes.size());

		for (uint32 i = 0; i < no_sort_start; ++i)
			random_bins.push_back(bin_sizes[i].first);

		for (uint32 i = no_sort_end; i < bin_sizes.size(); ++i)
			random_bins.push_back(bin_sizes[i].first);

		random_shuffle(random_bins.begin(), random_bins.end());

		for (uint32 i = no_sort_start; i < no_sort_end; ++i)
			random_bins.push_back(bin_sizes[i].first);
	}

	int32 get_next_random_bin()
	{
		lock_guard<mutex> lck(mtx);
		if (bin_id == -1)
			bin_id = 0;
		else
			++bin_id;

		if (bin_id >= (int32)m.size())
			return -1000;
		return random_bins[bin_id];
	}

	int32 get_next_bin()
	{
		lock_guard<mutex> lck(mtx);
		map_t::iterator p;
		if(bin_id == -1)
			p = m.begin();
		else
		{
			p = m.find(bin_id);
			if(p != m.end())
				++p;
		}

		if(p == m.end())
			bin_id = -1000;
		else
			bin_id = p->first;		

		return bin_id;
	}
	void insert(int32 bin_id, CMemDiskFile *file, string desc, int64 size, uint64 n_rec, uint64 n_plus_x_recs, uint64 n_super_kmers, uint32 buffer_size = 0, uint32 kmer_len = 0) {
		lock_guard<mutex> lck(mtx);

		map_t::iterator p = m.find(bin_id);
		if(p != m.end())
		{
			if(desc != "")
			{
				get<0>(m[bin_id]) = desc;
				get<5>(m[bin_id]) = file;
			}
			get<1>(m[bin_id]) += size;
			get<2>(m[bin_id]) += n_rec;
			get<6>(m[bin_id]) += n_plus_x_recs;
			get<7>(m[bin_id]) += n_super_kmers;
			if(buffer_size)
			{
				get<3>(m[bin_id]) = buffer_size;
				get<4>(m[bin_id]) = kmer_len;
			}
		}
		else
			m[bin_id] = std::make_tuple(desc, size, n_rec, buffer_size, kmer_len, file, n_plus_x_recs, n_super_kmers);
	}
	void read(int32 bin_id, CMemDiskFile *&file, string &desc, uint64 &size, uint64 &n_rec, uint64 &n_plus_x_recs, uint32 &buffer_size, uint32 &kmer_len) {
		lock_guard<mutex> lck(mtx);

		desc			= get<0>(m[bin_id]);
		file			= get<5>(m[bin_id]);
		size			= (uint64) get<1>(m[bin_id]);
		n_rec			= get<2>(m[bin_id]);
		buffer_size		= get<3>(m[bin_id]);
		kmer_len		= get<4>(m[bin_id]);
		n_plus_x_recs	= get<6>(m[bin_id]);
	}
	void read(int32 bin_id, CMemDiskFile *&file, string &desc, uint64 &size, uint64 &n_rec, uint64 &n_plus_x_recs, uint64 &n_super_kmers) {
		lock_guard<mutex> lck(mtx);

		desc			= get<0>(m[bin_id]);
		file			= get<5>(m[bin_id]);
		size			= (uint64) get<1>(m[bin_id]);
		n_rec			= get<2>(m[bin_id]);
		n_plus_x_recs	= get<6>(m[bin_id]);
		n_super_kmers		= get<7>(m[bin_id]);
	}
};

//************************************************************************************************************
class CBinQueue {
	typedef tuple<int32, uchar *, uint64, uint64> elem_t;
	typedef queue<elem_t, list<elem_t>> queue_t;
	queue_t q;

	int n_writers;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	CBinQueue(int _n_writers) {
		lock_guard<mutex> lck(mtx);
		n_writers = _n_writers;
	}
	~CBinQueue() {}

	bool empty() {
		lock_guard<mutex> lck(mtx);
		return q.empty();
	}
	bool completed() {
		lock_guard<mutex> lck(mtx);
		return q.empty() && !n_writers;
	}
	void mark_completed() {
		lock_guard<mutex> lck(mtx);
		n_writers--;
		if(n_writers == 0)
			cv_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *part, uint64 size, uint64 n_rec) {
		lock_guard<mutex> lck(mtx);
		bool was_empty = q.empty();
		q.push(std::make_tuple(bin_id, part, size, n_rec));
		if(was_empty)
			cv_queue_empty.notify_all();
	}
	int32 top_bin_id()
	{
		unique_lock<mutex> lck(mtx);
		return get<0>(q.front());
	}


	bool pop(int32 &bin_id, uchar *&part, uint64 &size, uint64 &n_rec, bool& is_allowed, const set<int>& allowed) {
		unique_lock<mutex> lck(mtx);

		cv_queue_empty.wait(lck, [this]{return !q.empty() || !n_writers; });

		if (q.empty())
			return false;

		bin_id = get<0>(q.front());

		is_allowed = true;
		if (allowed.find(bin_id) == allowed.end())
		{
			is_allowed = false;
			return true;
		}

		part = get<1>(q.front());
		size = get<2>(q.front());
		n_rec = get<3>(q.front());
		q.pop();

		return true;
	}

	bool pop_if_any(int32 &bin_id, uchar *&part, uint64 &size, uint64 &n_rec)
	{
		lock_guard<mutex> lck(mtx);
		
		if (q.empty())
			return false;

		bin_id = get<0>(q.front());
		part = get<1>(q.front());
		size = get<2>(q.front());
		n_rec = get<3>(q.front());
		q.pop();

		return true;
	}

	bool pop(int32 &bin_id, uchar *&part, uint64 &size, uint64 &n_rec) {
		unique_lock<mutex> lck(mtx);

		cv_queue_empty.wait(lck, [this]{return !q.empty() || !n_writers;}); 

		if(q.empty())
			return false;

		bin_id = get<0>(q.front());
		part   = get<1>(q.front());
		size   = get<2>(q.front());
		n_rec  = get<3>(q.front());
		q.pop();

		return true;
	}
};

//************************************************************************************************************
class CKmerQueue {
	typedef tuple<int32, uchar*, list<pair<uint64, uint64>>, uchar*, uint64, uint64, uint64, uint64, uint64> data_t;
	typedef list<data_t> list_t;
	int n_writers;
private:
	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

	list_t l;
	int32 n_bins;
public:
	CKmerQueue(int32 _n_bins, int _n_writers) {
		lock_guard<mutex> lck(mtx);
		n_bins = _n_bins;
		n_writers = _n_writers;
	}
	~CKmerQueue() {
	}

	bool empty() {
		lock_guard<mutex> lck(mtx);
		return l.empty() && !n_writers;
	}
	void mark_completed() {
		lock_guard<mutex> lck(mtx);
		n_writers--;
		if (!n_writers)
			cv_queue_empty.notify_all();
	}

	void add_n_writers(int32 n)
	{
		lock_guard<mutex> lck(mtx);
		n_writers += n;
	}

	void push(int32 bin_id, uchar *data, list<pair<uint64, uint64>> data_packs, uchar *lut, uint64 lut_size, uint64 n_unique, uint64 n_cutoff_min, uint64 n_cutoff_max, uint64 n_total) {
		lock_guard<mutex> lck(mtx);
		l.push_back(std::make_tuple(bin_id, data, std::move(data_packs), lut, lut_size, n_unique, n_cutoff_min, n_cutoff_max, n_total));
		cv_queue_empty.notify_all();
	}
	bool pop(int32 &bin_id, uchar *&data, list<pair<uint64, uint64>>& data_packs, uchar *&lut, uint64 &lut_size, uint64 &n_unique, uint64 &n_cutoff_min, uint64 &n_cutoff_max, uint64 &n_total) {
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this]{return !l.empty() || !n_writers; });
		if (l.empty())
			return false;

		bin_id = get<0>(l.front());
		data = get<1>(l.front());
		data_packs = std::move(get<2>(l.front()));
		lut = get<3>(l.front());
		lut_size = get<4>(l.front());
		n_unique = get<5>(l.front());
		n_cutoff_min = get<6>(l.front());
		n_cutoff_max = get<7>(l.front());
		n_total = get<8>(l.front());

		l.pop_front();

		if (l.empty())
			cv_queue_empty.notify_all();

		return true;
	}
};




//************************************************************************************************************
class CMemoryMonitor {
	uint64 max_memory;
	uint64 memory_in_use;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_memory_full;				// The condition to wait for

public:
	CMemoryMonitor(uint64 _max_memory) {
		lock_guard<mutex> lck(mtx);
		max_memory    = _max_memory;
		memory_in_use = 0;
	}
	~CMemoryMonitor() {
	}

	void increase(uint64 n) {
		unique_lock<mutex> lck(mtx);
		cv_memory_full.wait(lck, [this, n]{return memory_in_use + n <= max_memory;});
		memory_in_use += n;
	}
	void force_increase(uint64 n) {
		unique_lock<mutex> lck(mtx);
		cv_memory_full.wait(lck, [this, n]{return memory_in_use + n <= max_memory || memory_in_use == 0;});
		memory_in_use += n;
	}
	void decrease(uint64 n) {
		lock_guard<mutex> lck(mtx);
		memory_in_use -= n;
		cv_memory_full.notify_all();
	}
	void info(uint64 &_max_memory, uint64 &_memory_in_use)
	{
		lock_guard<mutex> lck(mtx);
		_max_memory    = max_memory;
		_memory_in_use = memory_in_use;
	}
};

//************************************************************************************************************
class CMemoryPool {
protected:
	int64 total_size;
	int64 part_size;
	int64 n_parts_total;
	int64 n_parts_free;

	uchar *buffer, *raw_buffer;
	uint32 *stack;

	mutable mutex mtx;							// The mutex to synchronise on
	condition_variable cv;						// The condition to wait for

public:
	CMemoryPool(int64 _total_size, int64 _part_size) {
		raw_buffer = nullptr;
		buffer = nullptr;
		stack  = nullptr;
		prepare(_total_size, _part_size);
	}
	~CMemoryPool() {
		release();
	}

	void prepare(int64 _total_size, int64 _part_size) {
		release();

		n_parts_total = _total_size / _part_size;
		part_size     = (_part_size + 15) / 16 * 16;			// to allow mapping pointer to int*
		n_parts_free  = n_parts_total;

		total_size = n_parts_total * part_size;

		raw_buffer = new uchar[total_size+64];
		buffer     = raw_buffer;
		while(((uint64) buffer) % 64)
			buffer++;

		stack = new uint32[n_parts_total];
		for(uint32 i = 0; i < n_parts_total; ++i)
			stack[i] = i;
	}

	void release(void) {
		if(raw_buffer)
			delete[] raw_buffer;
		raw_buffer = nullptr;
		buffer     = nullptr;

		if(stack)
			delete[] stack;
		stack = nullptr;
	}

	// Allocate memory buffer - uchar*
	void reserve(uchar* &part)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0;});

		part = buffer + stack[--n_parts_free]*part_size;
	}
	// Allocate memory buffer - char*
	void reserve(char* &part)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0;});

		part = (char*) (buffer + stack[--n_parts_free]*part_size);
	}
	// Allocate memory buffer - uint32*
	void reserve(uint32* &part)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0;});

		part = (uint32*) (buffer + stack[--n_parts_free]*part_size);
	}
	// Allocate memory buffer - uint64*
	void reserve(uint64* &part)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0;});

		part = (uint64*) (buffer + stack[--n_parts_free]*part_size);
	}
	// Allocate memory buffer - double*
	void reserve(double* &part)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0;});

		part = (double*) (buffer + stack[--n_parts_free]*part_size);
	}
	// Allocate memory buffer - float*
	void reserve(float* &part)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0; });

		part = (float*)(buffer + stack[--n_parts_free] * part_size);
	}

	// Deallocate memory buffer - uchar*
	void free(uchar* part)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32) ((part - buffer) / part_size);
		
		cv.notify_all();
	}
	// Deallocate memory buffer - char*
	void free(char* part)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32) (((uchar*) part - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - uint32*
	void free(uint32* part)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32) ((((uchar *) part) - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - uint64*
	void free(uint64* part)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32) ((((uchar *) part) - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - double*
	void free(double* part)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32) ((((uchar *) part) - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - float*
	void free(float* part)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32)((((uchar *)part) - buffer) / part_size);
		cv.notify_all();
	}
};


class CMemoryPoolWithBamSupport : public CMemoryPool
{
private:
	uint32 bam_current_id = 0;
	std::set<uint32_t> finished_ids;
public:
	CMemoryPoolWithBamSupport(int64 _total_size, int64 _part_size) : CMemoryPool(_total_size, _part_size) 
	{
		
	}

	void bam_notify_id_finished(uint32_t id)
	{
		lock_guard<mutex> lck(mtx);
		finished_ids.insert(id);
		std::set<uint32_t>::iterator it;
		while ((it = finished_ids.find(bam_current_id)) != finished_ids.end())
		{
			finished_ids.erase(it);
			++bam_current_id;
		}
		cv.notify_all();
	}

	// Allock memory buffer for bam input, assure at least one part available for PrepareForSplitter task of CBamTaskManager and for current id
	void bam_reserve_gunzip(uchar* &part, uint32_t id)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this, id] {
			if (id == bam_current_id)
				return n_parts_free > 1;
			else
				return n_parts_free > 2;
		});

		part = buffer + stack[--n_parts_free] * part_size;
	}
};

class CMemoryBins {
	int64 total_size;
	int64 free_size;
	uint32 n_bins;

	bool use_strict_mem;

	typedef std::tuple<uchar*, uchar*, uchar*, uchar*, uchar*, uchar*, uchar*, int64> bin_ptrs_t;

public:
	typedef enum{ mba_input_file, mba_input_array, mba_tmp_array, mba_suffix, mba_kxmer_counters, mba_lut } mba_t;

private:
	uchar *buffer, *raw_buffer;
	bin_ptrs_t *bin_ptrs;
	uint32 n_threads;

	map<uint64, uint64> map_reserved;

	mutable mutex mtx;							// The mutex to synchronise on
	condition_variable cv;						// The condition to wait for

public:
	CMemoryBins(int64 _total_size, uint32 _n_bins, bool _use_strict_mem, uint32 _n_threads) {
		raw_buffer = nullptr;
		buffer = nullptr;
		bin_ptrs = nullptr;
		use_strict_mem = _use_strict_mem;
		n_threads = _n_threads;
		prepare(_total_size, _n_bins);
	}
	~CMemoryBins() {
		release();
	}

	uint64 GetTotalSize()
	{
		return total_size;
	}

	int64 round_up_to_alignment(int64 x)
	{
		return (x + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;
	}

	void parallel_memory_init()
	{
		// Initialize the memory in parallel
		vector<thread> v_thr;
		for (uint32 i = 0; i < n_threads; ++i)
			v_thr.push_back(thread([&, i]{
			uint64 part_size = total_size / n_threads;
			for (uint64 pos = part_size * i; pos < part_size*(i + 1); pos += 1 << 20)
				buffer[pos] = 0;
		}));

		for (auto &x : v_thr)
			x.join();
	}

	void prepare(int64 _total_size, uint32 _n_bins) {
		release();

		n_bins = _n_bins;
		bin_ptrs = new bin_ptrs_t[n_bins];

		total_size = round_up_to_alignment(_total_size - n_bins * sizeof(bin_ptrs_t));
		free_size = total_size;

		raw_buffer = (uchar*)malloc(total_size + ALIGNMENT);
		buffer = raw_buffer;
		while (((uint64)buffer) % ALIGNMENT)
			buffer++;

		parallel_memory_init();

		map_reserved.clear();
		map_reserved[total_size] = 0;							// guard
	}

	void release(void) {
		if (raw_buffer)
			::free(raw_buffer);
		raw_buffer = nullptr;
		buffer = nullptr;

		if (bin_ptrs)
			delete[] bin_ptrs;
		bin_ptrs = nullptr;
	}

	void log(string info, uint64 size = 0)
	{
		return;			// !!! Log will be removed 

		cerr << info;
		if (size)
			cerr << " [" << size << "]: ";
		else
			cerr << ": ";
		for (auto &p : map_reserved)
			cerr << "(" << p.first << ", " << p.second << ")   ";
		cerr << endl;

	}

	// Prepare memory buffer for bin of given id
/*	bool init_old(uint32 bin_id, uint32 sorting_phases, int64 file_size, int64 kxmers_size, int64 out_buffer_size, int64 kxmer_counter_size, int64 lut_size)
	{
		unique_lock<mutex> lck(mtx);
		int64 part1_size;
		int64 part2_size;

		if (sorting_phases % 2 == 0)
		{
			part1_size = kxmers_size + kxmer_counter_size;
			part2_size = max(max(file_size, kxmers_size), out_buffer_size + lut_size);
		}
		else
		{
			part1_size = max(kxmers_size + kxmer_counter_size, file_size);
			part2_size = max(kxmers_size, out_buffer_size + lut_size);
		}
		int64 req_size = part1_size + part2_size;
		if (use_strict_mem && req_size > total_size) 				
		{
			return false;
		}
		uint64 found_pos;
		uint64 last_found_pos;

		// Look for space to insert
		cv.wait(lck, [&]() -> bool{
			found_pos = total_size;
			uint64 prev_end_pos = 0;
			for (auto &p : map_reserved)
				if (prev_end_pos + req_size < p.first)
				{
					found_pos = prev_end_pos;
					return true;
				}
				else
					prev_end_pos = p.first + p.second;

			// Reallocate memory for buffer if necessary
			if (map_reserved.size() == 1 && req_size > (int64) total_size)
			{
				::free(raw_buffer);
				total_size = round_up_to_alignment(req_size);
				free_size = total_size;

				raw_buffer = (uchar*)malloc(total_size + ALIGNMENT);
				buffer = raw_buffer;
				while (((uint64)buffer) % ALIGNMENT)
					buffer++;

				parallel_memory_init(buffer, total_size, n_threads);

				map_reserved.clear();
				map_reserved[total_size] = 0;
				found_pos = 0;
				return true;
			}

			return false;
		});

		// Reserve found free space
		map_reserved[found_pos] = req_size;

		uchar *base_ptr = get<0>(bin_ptrs[bin_id]) = buffer + found_pos;

		if (sorting_phases % 2 == 0)				// the result of sorting is in the same place as input
		{
			get<1>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			get<2>(bin_ptrs[bin_id]) = base_ptr;
			get<3>(bin_ptrs[bin_id]) = base_ptr + part1_size;
		}
		else
		{
			get<1>(bin_ptrs[bin_id]) = base_ptr;
			get<2>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			get<3>(bin_ptrs[bin_id]) = base_ptr;
		}
		get<4>(bin_ptrs[bin_id]) = base_ptr + part1_size;									// data
		get<5>(bin_ptrs[bin_id]) = get<4>(bin_ptrs[bin_id]) + out_buffer_size;
		if (kxmer_counter_size)
			get<6>(bin_ptrs[bin_id]) = base_ptr + kxmers_size;								//kxmers counter
		else
			get<6>(bin_ptrs[bin_id]) = nullptr;
		free_size -= req_size;
		get<7>(bin_ptrs[bin_id]) = req_size;

		return true;
	}
	*/
	// Prepare memory buffer for bin of given id - in fact alllocate only for the bin file size
	bool init(uint32 bin_id, uint32 sorting_phases, int64 file_size, int64 kxmers_size, int64 out_buffer_size, int64 kxmer_counter_size, int64 lut_size)
	{
		unique_lock<mutex> lck(mtx);
		int64 part1_size;
		int64 part2_size;

		if (sorting_phases % 2 == 0)
		{
			part1_size = kxmers_size + kxmer_counter_size;
			part2_size = max(max(file_size, kxmers_size), out_buffer_size + lut_size);
		}
		else
		{
			part1_size = max(kxmers_size + kxmer_counter_size, file_size);
			part2_size = max(kxmers_size, out_buffer_size + lut_size);
		}
		int64 req_size = part1_size + part2_size;
		if (use_strict_mem && req_size > total_size)
		{
			return false;
		}

		log("Init begin", file_size);

		uint64 found_pos;

		// Look for space to insert
		cv.wait(lck, [&]() -> bool{
			found_pos = total_size;
			uint64 prev_end_pos = 0;

			// Look for space for complete bin
			for (auto &p : map_reserved)
				if (prev_end_pos + req_size < p.first)
				{
					found_pos = prev_end_pos;

					// Reserve found free space
					map_reserved[found_pos] = req_size;
					return true;
				}
				else
					prev_end_pos = p.first + p.second;

			// Free space is small, so look for the space just for bin file
			for (auto &p : map_reserved)
				if (prev_end_pos + file_size < p.first)
				{
					found_pos = prev_end_pos;

					// Reserve found free space
					map_reserved[found_pos] = file_size;
					return true;
				}
				else
					prev_end_pos = p.first + p.second;

			// Reallocate memory for buffer if necessary
			if (map_reserved.size() == 1 && req_size >(int64) total_size)
			{
				::free(raw_buffer);
				total_size = round_up_to_alignment(req_size);
				free_size = total_size;

				raw_buffer = (uchar*)malloc(total_size + ALIGNMENT);
				buffer = raw_buffer;
				while (((uint64)buffer) % ALIGNMENT)
					buffer++;

				parallel_memory_init();

				map_reserved.clear();
				map_reserved[total_size] = 0;
				found_pos = 0;

				// Reserve found free space
				map_reserved[found_pos] = req_size;

				return true;
			}

			return false;
		});

		uchar *base_ptr = get<0>(bin_ptrs[bin_id]) = buffer + found_pos;

		if (map_reserved[found_pos] == (uint64)file_size)		// Only memory for bin file is allocated
		{
			get<1>(bin_ptrs[bin_id]) = base_ptr;

			free_size -= file_size;
		}
		else											// Memory for whole bin is allocated
		{
			if (sorting_phases % 2 == 0)				// the result of sorting is in the same place as input
			{
				get<1>(bin_ptrs[bin_id]) = base_ptr + part1_size;
				get<2>(bin_ptrs[bin_id]) = base_ptr;
				get<3>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			}
			else
			{
				get<1>(bin_ptrs[bin_id]) = base_ptr;
				get<2>(bin_ptrs[bin_id]) = base_ptr + part1_size;
				get<3>(bin_ptrs[bin_id]) = base_ptr;
			}
			get<4>(bin_ptrs[bin_id]) = base_ptr + part1_size;									// data
			get<5>(bin_ptrs[bin_id]) = get<4>(bin_ptrs[bin_id]) + out_buffer_size;
			if (kxmer_counter_size)
				get<6>(bin_ptrs[bin_id]) = base_ptr + kxmers_size;								// kxmers counter
			else
				get<6>(bin_ptrs[bin_id]) = nullptr;
			get<7>(bin_ptrs[bin_id]) = req_size;

			free_size -= req_size;
		}

		log("Init end");

		return true;
	}

	// Prepare memory buffer for bin of given id - in fact alllocate only for the bin file size
	bool extend(uint32 bin_id, uint32 sorting_phases, int64 file_size, int64 kxmers_size, int64 out_buffer_size, int64 kxmer_counter_size, int64 lut_size)
	{
		unique_lock<mutex> lck(mtx);
		int64 part1_size;
		int64 part2_size;

		if (sorting_phases % 2 == 0)
		{
			part1_size = kxmers_size + kxmer_counter_size;
			part2_size = max(max(file_size, kxmers_size), out_buffer_size + lut_size);
		}
		else
		{
			part1_size = max(kxmers_size + kxmer_counter_size, file_size);
			part2_size = max(kxmers_size, out_buffer_size + lut_size);
		}
		int64 req_size = part1_size + part2_size;

		log("Ext. begin", req_size);

		uint64 found_pos;
		bool must_reallocate = false;
		uint64 present_pos = get<0>(bin_ptrs[bin_id]) - buffer;

		// Case 1 - whole buffer was allocated during init()
		// Check whether we already have whole buffer allocated or not
		if (map_reserved[present_pos] == (uint64)req_size)
		{
			log("Ext-r end");
			return true;
		}

		// We have to extend the buffer
		
		// Look for space to insert
		cv.wait(lck, [&]() -> bool{
			auto p = map_reserved.find(present_pos);
			auto q = p;
			++q;

			// Check whether we can just extend the memory buffer to the required size
			if (p->first + req_size <= q->first)
			{
				must_reallocate = false;
				found_pos = present_pos;
				return true;
			}
			
			// Look for a free space of required size
			found_pos = total_size;
			uint64 prev_end_pos = 0;
			for (auto &p : map_reserved)
				if (prev_end_pos + req_size < p.first)
				{
					found_pos = prev_end_pos;
					must_reallocate = true;
					return true;
				}
				else
					prev_end_pos = p.first + p.second;

			return false;
		});

		// Case 2 - buffer must be extended but without reallocation
		if (!must_reallocate)
		{
			map_reserved[present_pos] = req_size;

			uchar *base_ptr = get<0>(bin_ptrs[bin_id]) = buffer + found_pos;

			if (sorting_phases % 2 == 0)				// the result of sorting is in the same place as input
			{
				get<1>(bin_ptrs[bin_id]) = base_ptr + part1_size;
				get<2>(bin_ptrs[bin_id]) = base_ptr;
				get<3>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			}
			else
			{
				get<1>(bin_ptrs[bin_id]) = base_ptr;
				get<2>(bin_ptrs[bin_id]) = base_ptr + part1_size;
				get<3>(bin_ptrs[bin_id]) = base_ptr;
			}
			get<4>(bin_ptrs[bin_id]) = base_ptr + part1_size;									// data
			get<5>(bin_ptrs[bin_id]) = get<4>(bin_ptrs[bin_id]) + out_buffer_size;
			if (kxmer_counter_size)
				get<6>(bin_ptrs[bin_id]) = base_ptr + kxmers_size;								//kxmers counter
			else
				get<6>(bin_ptrs[bin_id]) = nullptr;
			
			get<7>(bin_ptrs[bin_id]) = req_size;

			free_size += file_size;
			free_size -= req_size;
			
			// Sub-Case 2a - input data must be moved to correct position
			if (sorting_phases % 2 == 0)				// Must move input data
				memcpy(base_ptr + part1_size, base_ptr, file_size);
			// Sub-Case 2b - nothing has to be done with input data
			else
			{
				
			}

			log("Ext-e end");
			return true;
		}

		// Case 3 - Buffer is allocated in new place in memory, so we must move the input data
		map_reserved[found_pos] = req_size;

		auto readed_part = get<1>(bin_ptrs[bin_id]);

		uchar *base_ptr = get<0>(bin_ptrs[bin_id]) = buffer + found_pos;

		if (sorting_phases % 2 == 0)				// the result of sorting is in the same place as input
		{
			get<1>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			get<2>(bin_ptrs[bin_id]) = base_ptr;
			get<3>(bin_ptrs[bin_id]) = base_ptr + part1_size;
		}
		else
		{
			get<1>(bin_ptrs[bin_id]) = base_ptr;
			get<2>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			get<3>(bin_ptrs[bin_id]) = base_ptr;
		}
		get<4>(bin_ptrs[bin_id]) = base_ptr + part1_size;									// data
		get<5>(bin_ptrs[bin_id]) = get<4>(bin_ptrs[bin_id]) + out_buffer_size;
		if (kxmer_counter_size)
			get<6>(bin_ptrs[bin_id]) = base_ptr + kxmers_size;								//kxmers counter
		else
			get<6>(bin_ptrs[bin_id]) = nullptr;
		get<7>(bin_ptrs[bin_id]) = req_size;

		free_size += file_size;
		free_size -= req_size;

		// Make a copy of the readed part
		memcpy(get<1>(bin_ptrs[bin_id]), readed_part, file_size);

		// Remove the memory of the already copied part of data
		map_reserved.erase(present_pos);

		log("Ext-r end");
		return true;
	}


	void reserve(uint32 bin_id, uchar* &part, mba_t t)
	{		
		unique_lock<mutex> lck(mtx);
		if (t == mba_input_file)
			part = get<1>(bin_ptrs[bin_id]);
		else if (t == mba_input_array)
			part = get<2>(bin_ptrs[bin_id]);
		else if (t == mba_tmp_array)
			part = get<3>(bin_ptrs[bin_id]);
		else if (t == mba_suffix)
			part = get<4>(bin_ptrs[bin_id]);
		else if (t == mba_lut)
			part = get<5>(bin_ptrs[bin_id]);
		else if (t == mba_kxmer_counters)
			part = get<6>(bin_ptrs[bin_id]);
	}

	// Deallocate memory buffer - uchar*
	void free(uint32 bin_id, mba_t t)
	{		
		unique_lock<mutex> lck(mtx);
		if (t == mba_input_file)
			get<1>(bin_ptrs[bin_id]) = nullptr;
		else if (t == mba_input_array)
			get<2>(bin_ptrs[bin_id]) = nullptr;
		else if (t == mba_tmp_array)
			get<3>(bin_ptrs[bin_id]) = nullptr;
		else if (t == mba_suffix)
			get<4>(bin_ptrs[bin_id]) = nullptr;
		else if (t == mba_lut)
			get<5>(bin_ptrs[bin_id]) = nullptr;
		else if (t == mba_kxmer_counters)
			get<6>(bin_ptrs[bin_id]) = nullptr;

		if (!get<1>(bin_ptrs[bin_id]) && !get<2>(bin_ptrs[bin_id]) && !get<3>(bin_ptrs[bin_id]) && !get<4>(bin_ptrs[bin_id]) &&
			!get<5>(bin_ptrs[bin_id]) && !get<6>(bin_ptrs[bin_id]))
		{
			map_reserved.erase(get<0>(bin_ptrs[bin_id]) - buffer);

			log("Free");

			get<0>(bin_ptrs[bin_id]) = nullptr;
			free_size += get<7>(bin_ptrs[bin_id]);
			cv.notify_all();
		}
	}
};


class CTooLargeBinsQueue
{
	queue<int32, list<int32>> q;
	uint32 curr;
public:
	CTooLargeBinsQueue()
	{
		curr = 0;
	}

	~CTooLargeBinsQueue()
	{
	}


	bool get_next(int32& _bin_id)
	{
		if (q.empty())
			return false;
		_bin_id = q.front();
		q.pop();
		return true;
	}
	bool empty()
	{
		return q.empty();
	}
	void insert(int32 _bin_id)
	{
		q.push(_bin_id);
	}
};


class CBigBinPartQueue
{
	typedef std::tuple<int32, uchar*, uint64> data_t;
	typedef list<data_t> list_t;
	list_t l;	
	bool completed;
	mutable mutex mtx;
	condition_variable cv_pop;
public:
	void init()
	{
		completed = false;
	}
	CBigBinPartQueue()
	{
		init();
	}

	void push(int32 bin_id, uchar* data, uint64 size)
	{
		lock_guard<mutex> lck(mtx);
		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, data, size));
		if (was_empty)
			cv_pop.notify_all();
	}

	bool pop(int32& bin_id, uchar* &data, uint64& size)
	{
		unique_lock<mutex> lck(mtx);
		cv_pop.wait(lck, [this]{return !l.empty() || completed; });
		if (completed && l.empty())
			return false;
		
		bin_id = get<0>(l.front());
		data = get<1>(l.front());
		size = get<2>(l.front());
		l.pop_front();
		return true;
	}

	void mark_completed()
	{
		lock_guard<mutex> lck(mtx);
		completed = true;
		cv_pop.notify_all();
	}
};

class CBigBinKXmersQueue
{
	typedef std::tuple<int32, uchar*, uint64> data_t;
	typedef list<data_t> list_t;
	list_t l;
	uint32 n_writers;
	mutable mutex mtx;
	condition_variable cv_pop;	

	uint32 n_waiters;
	int32 current_id;
	condition_variable cv_push;	
public:
	CBigBinKXmersQueue(uint32 _n_writers)
	{
		n_waiters = 0;
		current_id = -1; //means queue is not initialized
		n_writers = _n_writers;
	}

	void push(int32 bin_id, uchar* data, uint64 size)
	{
		unique_lock<mutex> lck(mtx);
		++n_waiters;
		if (current_id == -1)
			current_id = bin_id;
		cv_push.wait(lck, [this, bin_id]{return bin_id == current_id || n_waiters == n_writers; });
		if (n_waiters == n_writers)
		{
			current_id = bin_id;
			cv_push.notify_all();
		}
		--n_waiters;
		
		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, data, size));
		if(was_empty)
			cv_pop.notify_all();
	}

	bool pop(int32& bin_id, uchar* &data, uint64& size)
	{
		unique_lock<mutex> lck(mtx);
		cv_pop.wait(lck, [this]{return !l.empty() || !n_writers; });
		if (l.empty() && !n_writers)
			return false;		
		bin_id = get<0>(l.front());
		data = get<1>(l.front());
		size = get<2>(l.front());
		l.pop_front();
		return true;
	}

	void mark_completed()
	{
		lock_guard<mutex> lck(mtx);
		--n_writers;
		if (!n_writers)
			cv_pop.notify_all();

		cv_push.notify_all();
	}

	~CBigBinKXmersQueue()
	{
	}

};

class CBigBinSortedPartQueue
{
	//bin_id, sub_bin_id,suff_buff, suff_buff_size, lut, lut_size, last_one_in_bin
	typedef std::tuple<int32, int32, uchar*, uint64, uint64*, uint64, bool> data_t;
	typedef list<data_t> list_t;
	list_t l;
	uint32 n_writers;
	mutable mutex mtx;
	condition_variable cv_pop;
public:
	CBigBinSortedPartQueue(uint32 _n_writers)
	{
		n_writers = _n_writers;
	}
	void push(int32 bin_id, int32 sub_bin_id, uchar* suff_buff, uint64 suff_buff_size, uint64* lut, uint64 lut_size, bool last_one_in_sub_bin)
	{
		lock_guard<mutex> lck(mtx);
		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, sub_bin_id, suff_buff, suff_buff_size, lut, lut_size, last_one_in_sub_bin));
		if (was_empty)
			cv_pop.notify_all();
	}
	bool pop(int32& bin_id, int32& sub_bin_id, uchar* &suff_buff, uint64& suff_buff_size, uint64* &lut, uint64 &lut_size, bool &last_one_in_sub_bin)
	{
		unique_lock<mutex> lck(mtx);
		cv_pop.wait(lck, [this]{return !n_writers || !l.empty(); });
		if (!n_writers && l.empty())
			return false;

		bin_id			= get<0>(l.front());
		sub_bin_id		= get<1>(l.front());
		suff_buff		= get<2>(l.front());
		suff_buff_size	= get<3>(l.front());
		lut				= get<4>(l.front());
		lut_size		= get<5>(l.front());
		last_one_in_sub_bin	= get<6>(l.front());

		l.pop_front();
		return true;
	}
	void mark_completed()
	{
		--n_writers;
		if (!n_writers)
			cv_pop.notify_all();
	}
};

class CBigBinKmerPartQueue
{
	typedef std::tuple<int32, uchar*, uint64, uchar*, uint64, uint64, uint64, uint64, uint64, bool> data_t;
	typedef list<data_t> list_t;
	list_t l;
	uint32 n_writers;
	mutable mutex mtx;
	condition_variable cv_pop;
	condition_variable cv_push;
	int32 curr_id;
	bool allow_next;
public:
	CBigBinKmerPartQueue(uint32 _n_writers)
	{
		n_writers = _n_writers;
		allow_next = true;
	}
	void push(int32 bin_id, uchar* suff_buff, uint64 suff_buff_size, uchar* lut, uint64 lut_size, uint64 n_unique, uint64 n_cutoff_min, uint64 n_cutoff_max, uint64 n_total, bool last_in_bin)
	{
		unique_lock<mutex> lck(mtx);	
		cv_push.wait(lck, [this, bin_id, lut_size]{return curr_id == bin_id || allow_next; });
		allow_next = false;
		if (last_in_bin)
		{
			allow_next = true;			
		}
		curr_id = bin_id;

		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, suff_buff, suff_buff_size, lut, lut_size, n_unique, n_cutoff_min, n_cutoff_max, n_total, last_in_bin));
		if (was_empty)
			cv_pop.notify_all();
		if (allow_next)
			cv_push.notify_all();
	}
	bool pop(int32& bin_id, uchar* &suff_buff, uint64& suff_buff_size, uchar* &lut, uint64& lut_size, uint64 &n_unique, uint64 &n_cutoff_min, uint64 &n_cutoff_max, uint64 &n_total, bool& last_in_bin)
	{
		unique_lock<mutex> lck(mtx);

		cv_pop.wait(lck, [this]{return !l.empty() || !n_writers; });
		if (!n_writers && l.empty())
			return false;
		bin_id = get<0>(l.front());
		suff_buff = get<1>(l.front());
		suff_buff_size = get<2>(l.front());
		lut = get<3>(l.front());
		lut_size = get<4>(l.front());
		n_unique = get<5>(l.front());
		n_cutoff_min = get<6>(l.front());
		n_cutoff_max = get<7>(l.front());
		n_total = get<8>(l.front());
		last_in_bin = get<9>(l.front());
		l.pop_front();
		return true;
	}

	void mark_completed()
	{
		lock_guard<mutex> lck(mtx);
		--n_writers;
		if(!n_writers)
			cv_pop.notify_all();
	}
};


class CBigBinDesc
{
	//lut_prefix_len, n_kmers, tmp_file_handle, string file_name, file_size
	typedef std::tuple<uint32, uint64, FILE*, string, uint64> elem_t;
	typedef map<int32, pair<int32, map<int32, elem_t>>> data_t;
	mutable mutex mtx;
	data_t m;
	int32 curr_id;
public:
	CBigBinDesc()
	{
		curr_id = -1;
	}
	void push(int32 bin_id, int32 sub_bin_id, uint32 lut_prefix_len, uint64 n_kmers, FILE* file, string desc, uint64 file_size)
	{
		lock_guard<mutex> lck(mtx);
		auto bin = m.find(bin_id);
		if (bin == m.end())
		{
			m[bin_id].first = -1;
			m[bin_id].second[sub_bin_id] = std::make_tuple(lut_prefix_len, n_kmers, file, desc, file_size);
		}
		else
		{
			auto sub_bin = bin->second.second.find(sub_bin_id);
			if (sub_bin == bin->second.second.end())
			{
				m[bin_id].second[sub_bin_id] = std::make_tuple(lut_prefix_len, n_kmers, file, desc, file_size);
			}
			else
			{				
				if(lut_prefix_len)
					get<0>(sub_bin->second) = lut_prefix_len;
				get<1>(sub_bin->second) += n_kmers;
				if (file)
				{
					get<2>(sub_bin->second) = file;
					get<3>(sub_bin->second) = desc;
				}
				get<4>(sub_bin->second) += file_size;
			}
		}		
	}

	bool get_n_sub_bins(int32 bin_id, uint32& size)
	{
		lock_guard<mutex> lck(mtx);
		auto e = m.find(bin_id);
		if (e == m.end())
			return false;

		size = (uint32)e->second.second.size();

		return true;
	}

	bool next_bin(int32& bin_id, uint32& size)
	{
		lock_guard<mutex> lck(mtx);		
		if (m.empty())
			return false;
		if (curr_id == -1)
		{
			curr_id = bin_id = m.begin()->first;
			size = (uint32)m.begin()->second.second.size();			
		}
		else
		{
			auto e = m.find(curr_id);

			e++;
			if (e == m.end())
				return false;
			curr_id = bin_id = e->first;
			size = (uint32)e->second.second.size();			
		}
		
		return true;
	}

	void reset_reading()
	{
		lock_guard<mutex> lck(mtx);
		curr_id = -1;
		for (auto& e : m)
			e.second.first = -1;			
	}

	bool next_sub_bin(int32 bin_id, int32& sub_bin_id, uint32& lut_prefix_len, uint64& n_kmers, FILE* &file, string& desc, uint64& file_size)
	{
		lock_guard<mutex> lck(mtx);
		auto& sub_bin = m.find(bin_id)->second;
		int32 curr_sub_bin_id = sub_bin.first;

		map<int32, elem_t>::iterator e;
		if (curr_sub_bin_id == -1)
			e = sub_bin.second.begin();
		else
		{
			e = sub_bin.second.find(curr_sub_bin_id);
			++e;
			if (e == sub_bin.second.end())
				return false;
		}
		sub_bin_id = sub_bin.first = e->first;
		lut_prefix_len = get<0>(e->second);
		n_kmers = get<1>(e->second);
		file = get<2>(e->second);
		desc = get<3>(e->second);		
		file_size = get<4>(e->second);
		return true;		
	}
};

class CCompletedBinsCollector
{
	list<int32> l;
	mutable mutex mtx;
	condition_variable cv_pop;
	uint32 n_writers;
public:
	CCompletedBinsCollector(uint32 _n_writers)
	{
		n_writers = _n_writers;
	}
	void push(int32 bin_id)
	{
		lock_guard<mutex> lck(mtx);
		bool was_empty = l.empty();
		l.push_back(bin_id);
		if (was_empty)
			cv_pop.notify_all();
	}

	bool pop(int32& bin_id)
	{
		unique_lock<mutex> lck(mtx);
		cv_pop.wait(lck, [this]{return !n_writers || !l.empty(); });
		if (!n_writers && l.empty())
			return false;

		bin_id = l.front();
		l.pop_front();
		return true;
	}
	void mark_completed()
	{
		lock_guard<mutex> lck(mtx);
		--n_writers;
		if (!n_writers)
			cv_pop.notify_all();
	}
};


class CDiskLogger
{
	uint64 current;
	uint64 max;
	mutable mutex mtx;
	
public:
	CDiskLogger()
	{
		current = max = 0;
	}
	void log_write(uint64 _size)
	{
		lock_guard<mutex> lck(mtx);
		current += _size;
		if (current > max)
			max = current;
	}
	void log_remove(uint64 _size)
	{
		lock_guard<mutex> lck(mtx);
		current -= _size;
	}
	uint64 get_max()
	{
		lock_guard<mutex> lck(mtx);
		return max;
	}
	uint64 get_current()
	{
		lock_guard<mutex> lck(mtx);
		return current;
	}
};

class CSortersManager
{
	int free_threads = 0;
	int max_sorters = 0;
	int working_with_additional;
	vector<int> n_sorters; // number of sorters working at the same time
	CBinQueue *bq;

	mutex mtx;
	condition_variable cv_get_next;	
public:
	CSortersManager(uint32 n_bins, uint32 n_threads, CBinQueue *_bq, int64 max_mem_size, const vector<pair<int32, int64>>& sorted_bins)
	{
		bq = _bq;
		n_sorters.resize(n_bins, 0);
		max_sorters = free_threads = n_threads;
		working_with_additional = 0;
		uint32 curr_sorters = 1;
		uint32 pos = 0;
		uint32 max_sorters = free_threads;

		while (curr_sorters < max_sorters)
		{
			if (pos >= sorted_bins.size())
				break;

			while (pos < sorted_bins.size())
			{
				if (sorted_bins[pos].second > max_mem_size / 2.0 / curr_sorters)
					n_sorters[sorted_bins[pos++].first] = /*max_sorters / */curr_sorters; //but in fact one more will be possible when available
				else
					break;
			}
			
			curr_sorters *= 2;
		}

		for (uint32 i = pos; i < sorted_bins.size(); ++i)
			n_sorters[sorted_bins[i].first] = max_sorters/*1*/;		
	}

	bool GetNext(int32 &bin_id, uchar *&part, uint64 &size, uint64 &n_rec, int& n_threads)
	{
		unique_lock<mutex> lck(mtx);
		
		bool no_more = false;
		bool poped = false;
		cv_get_next.wait(lck, [this, &bin_id, &part, &size, &n_rec, &no_more, &poped, &n_threads]
		{
			if (!poped)
			{
				poped = bq->pop_if_any(bin_id, part, size, n_rec);
				if (!poped)
					no_more = bq->completed();
			}
			if (no_more)
				return true;
			if (!poped)
				return false;

			int should_work_with_additional = max_sorters % n_sorters[bin_id];
			
			n_threads = max_sorters / n_sorters[bin_id];			
			if (working_with_additional < should_work_with_additional)
			{
				++n_threads;
			}
			return free_threads >= n_threads;
		});
		if (no_more)
			return false;
	
		free_threads -= n_threads;		

		if (max_sorters / n_sorters[bin_id] < n_threads)
			++working_with_additional;
		
		return true;
	}

	void ReturnThreads(uint32 n_threads, uint32 bin_id)
	{
		lock_guard<mutex> lck(mtx);		
		free_threads += n_threads;	
		if (max_sorters / n_sorters[bin_id] < (int)n_threads)
			--working_with_additional;		
		cv_get_next.notify_all();
	}

	void NotifyBQPush()
	{
		lock_guard<mutex> lck(mtx);	
		cv_get_next.notify_one();
	}
	void NotifyQueueCompleted()
	{
		lock_guard<mutex> lck(mtx);
		cv_get_next.notify_all();
	}
};

class CBamTaskManager
{
	mutex mtx;
	condition_variable cv;
	bool binary_reader_completed = false;
	bool ignore_rest = false;
	
	queue<tuple<uchar*, uint64, uint32, uint32>> bam_binary_part_queue; //data, size, id (of pack), file no

	bool splitters_preparer_is_working = false;
	
	//helper class
	class GunzippedQueue
	{
		map<uint32_t, queue<tuple<uchar*, uint64, uint32, uint32>>> _m; //data, size, id, file_no
		
		uint32_t current_id = 0;
		uint32_t last_id = (uint32_t)-1; //MAX_UINT32_T means last id is not known yet
		set<uint32_t> finished_ids;
	public:
		void Push(uchar* data, uint64 size, uint32 id, uint32 file_no)
		{
			_m[id].emplace(data, size, id, file_no); 
		}

		void NotifyIdFinished(uint32_t id)
		{
			finished_ids.insert(id);
		}

		void SetLastId(uint32_t value)
		{
			last_id = value;
		}

		enum class NextPackState { SUCCESS, NOT_YET_PRESENT, NO_MORE_PACKS };
		NextPackState GetNextPartState(uchar* &data, uint64 &size, uint32 &id, uint32 &file_no)
		{
			set<uint32_t>::iterator it;
			while (true)
			{
				auto& q = _m[current_id];
				if (!q.empty())
				{
					tie(data, size, id, file_no/*, test_id*/) = q.front();
					q.pop();

					return NextPackState::SUCCESS;
				}
				else if ((it = finished_ids.find(current_id)) != finished_ids.end()) //if current id is completed move to next id
				{
					finished_ids.erase(it);
					_m.erase(current_id);
					++current_id;
				}
				else if (last_id != (uint32_t)(-1) && current_id > last_id) //last id is known and current id is bigger than last -> completed
				{
					return NextPackState::NO_MORE_PACKS;
				}
				else
				{
					return NextPackState::NOT_YET_PRESENT; //pack is not yet present
				}
			}
		}

		void IgnoreRest(CMemoryPool* pmm_fastq)
		{
			for (auto& e : _m)
			{
				auto& q = e.second;
				while (!q.empty())
				{
					pmm_fastq->free(get<0>(q.front()));
					q.pop();
				}
			}
		}
	};

	GunzippedQueue gunzipped_queue;
public:
	
	enum TaskType { Gunzip, PrepareForSplitter}; 
	bool PushBinaryPack(uchar* data, uint64 size, uint32 id, uint32 file_no)
	{
		lock_guard<mutex> lck(mtx);
		if (ignore_rest)
			return false;
		bool was_empty = bam_binary_part_queue.empty();
		bam_binary_part_queue.emplace(data, size, id, file_no);
		if (was_empty)
			cv.notify_all();
		return true;
	}

	void NotifyBinaryReaderCompleted(uint32 last_id)
	{
		lock_guard<mutex> lck(mtx);
		gunzipped_queue.SetLastId(last_id);
		binary_reader_completed = true;
		cv.notify_all();
	}

	bool PushGunzippedPart(uchar* data, uint64 size, uint32 id, uint32 file_no)
	{
		lock_guard<mutex> lck(mtx);

		if (ignore_rest)
			return false;
		gunzipped_queue.Push(data, size, id, file_no);

		cv.notify_all();
		return true;
	}

	void NotifyIdFinished(uint32 id)
	{
		lock_guard<mutex> lck(mtx);
		gunzipped_queue.NotifyIdFinished(id);
		cv.notify_all();
	}

	
	void NotifySplitterPrepareTaskDone()
	{
		lock_guard<mutex> lck(mtx);
		splitters_preparer_is_working = false;
		cv.notify_all();
	}

	bool TakeNextPrepareForSplitterTaskIfExists(uchar* &data, uint64 &size, uint32 &id, uint32 &file_no)
	{
		lock_guard<mutex> lck(mtx);
		if (splitters_preparer_is_working)
			return false;
		auto state = gunzipped_queue.GetNextPartState(data, size, id, file_no);
		if (state == GunzippedQueue::NextPackState::SUCCESS)
		{
			splitters_preparer_is_working = true;
			return true;
		}
		return false;
	}

	bool PopTask(TaskType& taskType, uchar* &data, uint64 &size, uint32 &id, uint32 &file_no)
	{
		unique_lock<mutex> lck(mtx);
		GunzippedQueue::NextPackState state = GunzippedQueue::NextPackState::NOT_YET_PRESENT;
		cv.wait(lck, [this, &data, &size, &id, &file_no, &state]
		{
			if (ignore_rest)
				return true;
			if (!splitters_preparer_is_working) //only one thread can prepare parts for splitters
			{
				state = gunzipped_queue.GetNextPartState(data, size, id, file_no);
				if (state == GunzippedQueue::NextPackState::SUCCESS)
					return true;
			}
			if (!bam_binary_part_queue.empty())
				return true;

			if(binary_reader_completed && state == GunzippedQueue::NextPackState::NO_MORE_PACKS)
				return true;
			return false;			
		});
		
		if (ignore_rest)
			return false;

		if (state == GunzippedQueue::NextPackState::SUCCESS)
		{
			splitters_preparer_is_working = true;
			taskType = TaskType::PrepareForSplitter;
			return true;
		}

		if (!bam_binary_part_queue.empty())
		{
			taskType = TaskType::Gunzip;
			tie(data, size, id, file_no) = bam_binary_part_queue.front();
			bam_binary_part_queue.pop();
			return true;
		}

		return false;
	}

	void IgnoreRest(CMemoryPool* pmm_fastq, CMemoryPool* pmm_binary_file_reader)
	{
		lock_guard<mutex> lck(mtx);
		ignore_rest = true;
		while (!bam_binary_part_queue.empty())
		{
			pmm_binary_file_reader->free(get<0>(bam_binary_part_queue.front()));
			bam_binary_part_queue.pop();
		}
		gunzipped_queue.IgnoreRest(pmm_fastq);

		cv.notify_all();
	}

	struct SSplitterPrepareState
	{
		uint32 current_file_no = (uint32)(-1); //-1 means not started
		uchar* prev_part_data = nullptr;
		uint64 prev_part_size = 0;
	} splitter_prepare_state;
};

#endif

// ***** EOF
