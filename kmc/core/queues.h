/*
    This file is a part of KMC software distributed under GNU GPL 3 licence.
    The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

    Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

    Version: 2.2.0
    Date   : 2015-04-15
*/

#ifndef _QUEUES_H
#define _QUEUES_H

#include "kmc_typedefs.h"

#include "mem_disk_file.h"

#include <stdio.h>
#include <iostream>
#include <tuple>
#include <queue>
#include <list>
#include <map>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>


//************************************************************************************************************
class CInputFilesQueue {
	typedef std::string elem_t;
	typedef std::queue<elem_t, std::list<elem_t>> queue_t;

	queue_t q;
	bool is_completed;

	mutable std::mutex mtx;								// The std::mutex to synchronise on

  public:
	CInputFilesQueue(const std::vector<std::string> &file_names) {
		std::unique_lock<std::mutex> lck(mtx);

		for(std::vector<std::string>::const_iterator p = file_names.begin(); p != file_names.end(); ++p)
			q.push(*p);

		is_completed = false;
	};
	~CInputFilesQueue() {};

	bool empty() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty();
	}
	bool completed() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty() && is_completed;
	}
	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		is_completed = true;
	}
	bool pop(std::string &file_name) {
		std::lock_guard<std::mutex> lck(mtx);

		if(q.empty())
			return false;

		file_name = q.front();
		q.pop();
		return true;
	}
};

//************************************************************************************************************
class CPartQueue {
	typedef std::pair<uchar *, uint64> elem_t;
	typedef std::queue<elem_t, std::list<elem_t>> queue_t;

	queue_t q;
	bool is_completed;
	int n_readers;

	mutable std::mutex mtx;								// The std::mutex to synchronise on
	std::condition_variable cv_queue_empty;

  public:
	CPartQueue(int _n_readers) {
		std::unique_lock<std::mutex> lck(mtx);
		is_completed    = false;
		n_readers       = _n_readers;
	};
	~CPartQueue() {};

	bool empty() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty();
	}
	bool completed() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty() && !n_readers;
	}
	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		n_readers--;

		if(!n_readers)
			cv_queue_empty.notify_all();
	}
	void push(uchar *part, uint64 size) {
		std::unique_lock<std::mutex> lck(mtx);
		bool was_empty = q.empty();
		q.push(std::make_pair(part, size));

		if(was_empty)
			cv_queue_empty.notify_all();
	}
	bool pop(uchar *&part, uint64 &size) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !this->q.empty() || !this->n_readers;});

		if(q.empty())
			return false;

		part = q.front().first;
		size = q.front().second;
		q.pop();
		return true;
	}
};

//************************************************************************************************************
class CStatsPartQueue {
	typedef std::pair<uchar *, uint64> elem_t;
	typedef std::queue<elem_t, std::list<elem_t>> queue_t;

	queue_t q;

	mutable std::mutex mtx;
	std::condition_variable cv_queue_empty;
	int n_readers;
	int64 bytes_to_read;
  public:
	CStatsPartQueue(int _n_readers, int64 _bytes_to_read) {
		std::unique_lock<std::mutex> lck(mtx);
		n_readers = _n_readers;
		bytes_to_read = _bytes_to_read;
	}

	~CStatsPartQueue() {};

	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		n_readers--;

		if (!n_readers)
			cv_queue_empty.notify_all();
	}

	bool completed() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty() && !n_readers;
	}

	bool push(uchar *part, uint64 size) {
		std::unique_lock<std::mutex> lck(mtx);

		if (bytes_to_read <= 0)
			return false;

		bool was_empty = q.empty();
		q.push(std::make_pair(part, size));
		bytes_to_read -= size;

		if (was_empty)
			cv_queue_empty.notify_one();

		return true;
	}

	bool pop(uchar *&part, uint64 &size) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !this->q.empty() || !this->n_readers; });

		if (q.empty())
			return false;

		part = q.front().first;
		size = q.front().second;
		q.pop();
		return true;
	}


};

//************************************************************************************************************
class CBinPartQueue {
	typedef std::tuple<int32, uchar *, uint32, uint32> elem_t;
	typedef std::queue<elem_t, std::list<elem_t>> queue_t;
	queue_t q;

	int n_writers;
	bool is_completed;

	mutable std::mutex mtx;						// The std::mutex to synchronise on
	std::condition_variable cv_queue_empty;

  public:
	CBinPartQueue(int _n_writers) {
		std::lock_guard<std::mutex> lck(mtx);
		n_writers       = _n_writers;
		is_completed    = false;
	}
	~CBinPartQueue() {}

	bool empty() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty();
	}
	bool completed() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty() && !n_writers;
	}
	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		n_writers--;

		if(!n_writers)
			cv_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *part, uint32 true_size, uint32 alloc_size) {
		std::unique_lock<std::mutex> lck(mtx);
		bool was_empty = q.empty();
		q.push(std::make_tuple(bin_id, part, true_size, alloc_size));

		if(was_empty)
			cv_queue_empty.notify_all();
	}
	bool pop(int32 &bin_id, uchar *&part, uint32 &true_size, uint32 &alloc_size) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !q.empty() || !n_writers;});

		if(q.empty())
			return false;

		bin_id     = std::get<0>(q.front());
		part       = std::get<1>(q.front());
		true_size  = std::get<2>(q.front());
		alloc_size = std::get<3>(q.front());
		q.pop();
		return true;
	}
};

//************************************************************************************************************
class CBinDesc {
	typedef std::tuple<std::string, int64, uint64, uint32, uint32, CMemDiskFile*, uint64, uint64> desc_t;
	typedef std::map<int32, desc_t> map_t;

	map_t m;
	int32 bin_id;

	std::vector<int32> random_bins;

	mutable std::mutex mtx;

  public:
	CBinDesc() {
		std::lock_guard<std::mutex> lck(mtx);
		bin_id = -1;
	}
	~CBinDesc() {}

	void reset_reading() {
		std::lock_guard<std::mutex> lck(mtx);
		bin_id = -1;
	}

	bool empty() {
		std::lock_guard<std::mutex> lck(mtx);
		return m.empty();
	}

	void init_random() {
		std::lock_guard<std::mutex> lck(mtx);
		std::vector<std::pair<int32, int64>> bin_sizes;

		for (auto& p : m)
			bin_sizes.push_back(std::make_pair(p.first, std::get<2>(p.second)));

		sort(bin_sizes.begin(), bin_sizes.end(), [](const std::pair<int32, int64>& l,
		const std::pair<int32, int64>& r) {
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

	int32 get_next_random_bin() {
		std::lock_guard<std::mutex> lck(mtx);

		if (bin_id == -1)
			bin_id = 0;
		else
			++bin_id;

		if (bin_id >= (int32)m.size())
			return -1000;

		return random_bins[bin_id];
	}

	int32 get_next_bin() {
		std::lock_guard<std::mutex> lck(mtx);
		map_t::iterator p;

		if(bin_id == -1)
			p = m.begin();
		else {
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
	void insert(int32 bin_id, CMemDiskFile *file, std::string desc, int64 size, uint64 n_rec,
				uint64 n_plus_x_recs, uint64 n_super_kmers, uint32 buffer_size = 0, uint32 kmer_len = 0) {
		std::lock_guard<std::mutex> lck(mtx);
		map_t::iterator p = m.find(bin_id);

		if(p != m.end()) {
			if(desc != "") {
				std::get<0>(m[bin_id]) = desc;
				std::get<5>(m[bin_id]) = file;
			}

			std::get<1>(m[bin_id]) += size;
			std::get<2>(m[bin_id]) += n_rec;
			std::get<6>(m[bin_id]) += n_plus_x_recs;
			std::get<7>(m[bin_id]) += n_super_kmers;

			if(buffer_size) {
				std::get<3>(m[bin_id]) = buffer_size;
				std::get<4>(m[bin_id]) = kmer_len;
			}
		} else
			m[bin_id] = std::make_tuple(desc, size, n_rec, buffer_size, kmer_len, file, n_plus_x_recs, n_super_kmers);
	}
	void read(int32 bin_id, CMemDiskFile *&file, std::string &desc, uint64 &size, uint64 &n_rec,
			  uint64 &n_plus_x_recs, uint32 &buffer_size, uint32 &kmer_len) {
		std::lock_guard<std::mutex> lck(mtx);
		desc			= std::get<0>(m[bin_id]);
		file			= std::get<5>(m[bin_id]);
		size			= (uint64) std::get<1>(m[bin_id]);
		n_rec			= std::get<2>(m[bin_id]);
		buffer_size		= std::get<3>(m[bin_id]);
		kmer_len		= std::get<4>(m[bin_id]);
		n_plus_x_recs	= std::get<6>(m[bin_id]);
	}
	void read(int32 bin_id, CMemDiskFile *&file, std::string &desc, uint64 &size, uint64 &n_rec,
			  uint64 &n_plus_x_recs, uint64 &n_super_kmers) {
		std::lock_guard<std::mutex> lck(mtx);
		desc			= std::get<0>(m[bin_id]);
		file			= std::get<5>(m[bin_id]);
		size			= (uint64) std::get<1>(m[bin_id]);
		n_rec			= std::get<2>(m[bin_id]);
		n_plus_x_recs	= std::get<6>(m[bin_id]);
		n_super_kmers		= std::get<7>(m[bin_id]);
	}
};

//************************************************************************************************************
class CBinQueue {
	typedef std::tuple<int32, uchar *, uint64, uint64> elem_t;
	typedef std::queue<elem_t, std::list<elem_t>> queue_t;
	queue_t q;

	int n_writers;

	mutable std::mutex mtx;								// The std::mutex to synchronise on
	std::condition_variable cv_queue_empty;

  public:
	CBinQueue(int _n_writers) {
		std::lock_guard<std::mutex> lck(mtx);
		n_writers = _n_writers;
	}
	~CBinQueue() {}

	bool empty() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty();
	}
	bool completed() {
		std::lock_guard<std::mutex> lck(mtx);
		return q.empty() && !n_writers;
	}
	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		n_writers--;

		if(n_writers == 0)
			cv_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *part, uint64 size, uint64 n_rec) {
		std::lock_guard<std::mutex> lck(mtx);
		bool was_empty = q.empty();
		q.push(std::make_tuple(bin_id, part, size, n_rec));

		if(was_empty)
			cv_queue_empty.notify_all();
	}
	bool pop(int32 &bin_id, uchar *&part, uint64 &size, uint64 &n_rec) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !q.empty() || !n_writers;});

		if(q.empty())
			return false;

		bin_id = std::get<0>(q.front());
		part   = std::get<1>(q.front());
		size   = std::get<2>(q.front());
		n_rec  = std::get<3>(q.front());
		q.pop();
		return true;
	}
};

//************************************************************************************************************
class CKmerQueue {
	typedef std::tuple<int32, uchar*, uint64, uchar*, uint64, uint64, uint64, uint64, uint64> data_t;
	typedef std::list<data_t> list_t;

	int n_writers;
	mutable std::mutex mtx;								// The std::mutex to synchronise on
	std::condition_variable cv_queue_empty;

	list_t l;
	int32 n_bins;
  public:
	CKmerQueue(int32 _n_bins, int _n_writers) {
		std::lock_guard<std::mutex> lck(mtx);
		n_bins = _n_bins;
		n_writers = _n_writers;
	}
	~CKmerQueue() {
	}

	bool empty() {
		std::lock_guard<std::mutex> lck(mtx);
		return l.empty() && !n_writers;
	}
	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		n_writers--;

		if (!n_writers)
			cv_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *data, uint64 data_size, uchar *lut, uint64 lut_size, uint64 n_unique,
			  uint64 n_cutoff_min, uint64 n_cutoff_max, uint64 n_total) {
		std::lock_guard<std::mutex> lck(mtx);
		l.push_back(std::make_tuple(bin_id, data, data_size, lut, lut_size, n_unique, n_cutoff_min, n_cutoff_max,
									n_total));
		cv_queue_empty.notify_all();
	}
	bool pop(int32 &bin_id, uchar *&data, uint64 &data_size, uchar *&lut, uint64 &lut_size, uint64 &n_unique,
			 uint64 &n_cutoff_min, uint64 &n_cutoff_max, uint64 &n_total) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !l.empty() || !n_writers; });

		if (l.empty())
			return false;

		bin_id = std::get<0>(l.front());
		data = std::get<1>(l.front());
		data_size = std::get<2>(l.front());
		lut = std::get<3>(l.front());
		lut_size = std::get<4>(l.front());
		n_unique = std::get<5>(l.front());
		n_cutoff_min = std::get<6>(l.front());
		n_cutoff_max = std::get<7>(l.front());
		n_total = std::get<8>(l.front());
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

	mutable std::mutex mtx;								// The std::mutex to synchronise on
	std::condition_variable cv_memory_full;				// The condition to wait for

  public:
	CMemoryMonitor(uint64 _max_memory) {
		std::lock_guard<std::mutex> lck(mtx);
		max_memory    = _max_memory;
		memory_in_use = 0;
	}
	~CMemoryMonitor() {
	}

	void increase(uint64 n) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_memory_full.wait(lck, [this, n] {return memory_in_use + n <= max_memory;});
		memory_in_use += n;
	}
	void force_increase(uint64 n) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_memory_full.wait(lck, [this, n] {return memory_in_use + n <= max_memory || memory_in_use == 0;});
		memory_in_use += n;
	}
	void decrease(uint64 n) {
		std::lock_guard<std::mutex> lck(mtx);
		memory_in_use -= n;
		cv_memory_full.notify_all();
	}
	void info(uint64 &_max_memory, uint64 &_memory_in_use) {
		std::lock_guard<std::mutex> lck(mtx);
		_max_memory    = max_memory;
		_memory_in_use = memory_in_use;
	}
};

//************************************************************************************************************
class CMemoryPool {
	int64 total_size;
	int64 part_size;
	int64 n_parts_total;
	int64 n_parts_free;

	uchar *buffer, *raw_buffer;
	uint32 *stack;

	mutable std::mutex mtx;							// The std::mutex to synchronise on
	std::condition_variable cv;						// The condition to wait for

  public:
	CMemoryPool(int64 _total_size, int64 _part_size) {
		raw_buffer = NULL;
		buffer = NULL;
		stack  = NULL;
		prepare(_total_size, _part_size);
	}
	~CMemoryPool() {
		release();
	}

	void prepare(int64 _total_size, int64 _part_size) {
		release();
		n_parts_total = _total_size / _part_size;
		part_size     = (_part_size + 15) / 16 * 16;			// to allow std::mapping pointer to int*
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

		raw_buffer = NULL;
		buffer     = NULL;

		if(stack)
			delete[] stack;

		stack = NULL;
	}

	// Allocate memory buffer - uchar*
	void reserve(uchar* &part) {
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this] {return n_parts_free > 0;});
		part = buffer + stack[--n_parts_free]*part_size;
	}
	// Allocate memory buffer - char*
	void reserve(char* &part) {
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this] {return n_parts_free > 0;});
		part = (char*) (buffer + stack[--n_parts_free]*part_size);
	}
	// Allocate memory buffer - uint32*
	void reserve(uint32* &part) {
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this] {return n_parts_free > 0;});
		part = (uint32*) (buffer + stack[--n_parts_free]*part_size);
	}
	// Allocate memory buffer - uint64*
	void reserve(uint64* &part) {
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this] {return n_parts_free > 0;});
		part = (uint64*) (buffer + stack[--n_parts_free]*part_size);
	}
	// Allocate memory buffer - double*
	void reserve(double* &part) {
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this] {return n_parts_free > 0;});
		part = (double*) (buffer + stack[--n_parts_free]*part_size);
	}

	// Deallocate memory buffer - uchar*
	void free(uchar* part) {
		std::lock_guard<std::mutex> lck(mtx);
		stack[n_parts_free++] = (uint32) ((part - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - char*
	void free(char* part) {
		std::lock_guard<std::mutex> lck(mtx);
		stack[n_parts_free++] = (uint32) (((uchar*) part - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - uint32*
	void free(uint32* part) {
		std::lock_guard<std::mutex> lck(mtx);
		stack[n_parts_free++] = (uint32) ((((uchar *) part) - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - uint64*
	void free(uint64* part) {
		std::lock_guard<std::mutex> lck(mtx);
		stack[n_parts_free++] = (uint32) ((((uchar *) part) - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - double*
	void free(double* part) {
		std::lock_guard<std::mutex> lck(mtx);
		stack[n_parts_free++] = (uint32) ((((uchar *) part) - buffer) / part_size);
		cv.notify_all();
	}
};


class CMemoryBins {
	int64 total_size;
	int64 free_size;
	uint32 n_bins;

	bool use_strict_mem;

	typedef std::tuple<uchar*, uchar*, uchar*, uchar*, uchar*, uchar*, uchar*, int64> bin_ptrs_t;

  public:
	typedef enum { mba_input_file, mba_input_array, mba_tmp_array, mba_suffix, mba_kxmer_counters, mba_lut } mba_t;

  private:
	uchar *buffer, *raw_buffer;
	bin_ptrs_t *bin_ptrs;

	std::list<std::pair<uint64, uint64>> list_reserved;
	std::list<std::pair<uint32, uint64>> list_insert_order;

	mutable std::mutex mtx;							// The std::mutex to synchronise on
	std::condition_variable cv;						// The condition to wait for

  public:
	CMemoryBins(int64 _total_size, uint32 _n_bins, bool _use_strict_mem) {
		raw_buffer = NULL;
		buffer = NULL;
		bin_ptrs = NULL;
		use_strict_mem = _use_strict_mem;
		prepare(_total_size, _n_bins);
	}
	~CMemoryBins() {
		release();
	}

	int64 round_up_to_alignment(int64 x) {
		return (x + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;
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

		list_reserved.clear();
		list_insert_order.clear();
		list_reserved.push_back(std::make_pair(total_size, 0));		// guard
	}

	void release(void) {
		if (raw_buffer)
			::free(raw_buffer);

		raw_buffer = NULL;
		buffer = NULL;

		if (bin_ptrs)
			delete[] bin_ptrs;

		bin_ptrs = NULL;
	}

	// Prepare memory buffer for bin of given id
	bool init(uint32 bin_id, uint32 sorting_phases, int64 file_size, int64 kxmers_size, int64 out_buffer_size,
			  int64 kxmer_counter_size, int64 lut_size) {
		std::unique_lock<std::mutex> lck(mtx);
		int64 part1_size;
		int64 part2_size;

		if (sorting_phases % 2 == 0) {
			part1_size = kxmers_size + kxmer_counter_size;
			part2_size = std::max(std::max(file_size, kxmers_size), out_buffer_size + lut_size);
		} else {
			part1_size = std::max(kxmers_size + kxmer_counter_size, file_size);
			part2_size = std::max(kxmers_size, out_buffer_size + lut_size);
		}

		int64 req_size = part1_size + part2_size;

		if (use_strict_mem && req_size > total_size) {
			return false;
		}

		uint64 found_pos;
		uint64 last_found_pos;
		// Look for space to insert
		cv.wait(lck, [&]() -> bool{
			found_pos = total_size;

			if (!list_insert_order.empty()) {
				last_found_pos = list_insert_order.back().second;

				for (auto p = list_reserved.begin(); p != list_reserved.end(); ++p)
					if (p->first == last_found_pos) {
						uint64 last_end_pos = p->first + p->second;
						++p;

						if (last_end_pos + req_size <= p->first) {
							found_pos = last_end_pos;
							return true;
						} else
							break;
					}
			}

			uint64 prev_end_pos = 0;

			for (auto p = list_reserved.begin(); p != list_reserved.end(); ++p) {
				if (prev_end_pos + req_size <= p->first) {
					found_pos = prev_end_pos;
					return true;
				}

				prev_end_pos = p->first + p->second;
			}

			// Reallocate memory for buffer if necessary
			if (list_insert_order.empty() && req_size > (int64)list_reserved.back().first) {
				::free(raw_buffer);
				total_size = round_up_to_alignment(req_size);
				free_size = total_size;
				raw_buffer = (uchar*)malloc(total_size + ALIGNMENT);
				buffer = raw_buffer;

				while (((uint64)buffer) % ALIGNMENT)
					buffer++;

				list_reserved.back().first = total_size;
				found_pos = 0;
				return true;
			}

			return false;
		});
		// Reserve found free space
		list_insert_order.push_back(std::make_pair(bin_id, found_pos));

		for (auto p = list_reserved.begin(); p != list_reserved.end(); ++p)
			if (found_pos < p->first) {
				list_reserved.insert(p, std::make_pair(found_pos, req_size));
				break;
			}

		uchar *base_ptr = std::get<0>(bin_ptrs[bin_id]) = buffer + found_pos;

		if (sorting_phases % 2 == 0) {			// the result of sorting is in the same place as input
			std::get<1>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			std::get<2>(bin_ptrs[bin_id]) = base_ptr;
			std::get<3>(bin_ptrs[bin_id]) = base_ptr + part1_size;
		} else {
			std::get<1>(bin_ptrs[bin_id]) = base_ptr;
			std::get<2>(bin_ptrs[bin_id]) = base_ptr + part1_size;
			std::get<3>(bin_ptrs[bin_id]) = base_ptr;
		}

		std::get<4>(bin_ptrs[bin_id]) = base_ptr + part1_size;									// data
		std::get<5>(bin_ptrs[bin_id]) = std::get<4>(bin_ptrs[bin_id]) + out_buffer_size;

		if (kxmer_counter_size)
			std::get<6>(bin_ptrs[bin_id]) = base_ptr + kxmers_size;								//kxmers counter
		else
			std::get<6>(bin_ptrs[bin_id]) = NULL;

		free_size -= req_size;
		std::get<7>(bin_ptrs[bin_id]) = req_size;
		return true;
	}

	void reserve(uint32 bin_id, uchar* &part, mba_t t) {
		std::unique_lock<std::mutex> lck(mtx);

		if (t == mba_input_file)
			part = std::get<1>(bin_ptrs[bin_id]);
		else if (t == mba_input_array)
			part = std::get<2>(bin_ptrs[bin_id]);
		else if (t == mba_tmp_array)
			part = std::get<3>(bin_ptrs[bin_id]);
		else if (t == mba_suffix)
			part = std::get<4>(bin_ptrs[bin_id]);
		else if (t == mba_lut)
			part = std::get<5>(bin_ptrs[bin_id]);
		else if (t == mba_kxmer_counters)
			part = std::get<6>(bin_ptrs[bin_id]);
	}

	// Deallocate memory buffer - uchar*
	void free(uint32 bin_id, mba_t t) {
		std::unique_lock<std::mutex> lck(mtx);

		if (t == mba_input_file)
			std::get<1>(bin_ptrs[bin_id]) = NULL;
		else if (t == mba_input_array)
			std::get<2>(bin_ptrs[bin_id]) = NULL;
		else if (t == mba_tmp_array)
			std::get<3>(bin_ptrs[bin_id]) = NULL;
		else if (t == mba_suffix)
			std::get<4>(bin_ptrs[bin_id]) = NULL;
		else if (t == mba_lut)
			std::get<5>(bin_ptrs[bin_id]) = NULL;
		else if (t == mba_kxmer_counters)
			std::get<6>(bin_ptrs[bin_id]) = NULL;

		if (!std::get<1>(bin_ptrs[bin_id]) && !std::get<2>(bin_ptrs[bin_id]) && !std::get<3>(bin_ptrs[bin_id])
				&& !std::get<4>(bin_ptrs[bin_id]) && !std::get<5>(bin_ptrs[bin_id]) && !std::get<6>(bin_ptrs[bin_id])) {
			for (auto p = list_reserved.begin(); p != list_reserved.end() && p->second != 0; ++p) {
				if ((int64)p->first == std::get<0>(bin_ptrs[bin_id]) - buffer) {
					list_reserved.erase(p);
					break;
				}
			}

			for (auto p = list_insert_order.begin(); p != list_insert_order.end(); ++p)
				if (p->first == bin_id) {
					list_insert_order.erase(p);
					break;
				}

			std::get<0>(bin_ptrs[bin_id]) = NULL;
			free_size += std::get<7>(bin_ptrs[bin_id]);
			cv.notify_all();
		}
	}
};


class CTooLargeBinsQueue {
	std::queue<int32, std::list<int32>> q;
	uint32 curr;
  public:
	CTooLargeBinsQueue() {
		curr = 0;
	}

	~CTooLargeBinsQueue() {
	}


	bool get_next(int32& _bin_id) {
		if (q.empty())
			return false;

		_bin_id = q.front();
		q.pop();
		return true;
	}
	bool empty() {
		return q.empty();
	}
	void insert(int32 _bin_id) {
		q.push(_bin_id);
	}
};


class CBigBinPartQueue {
	typedef std::tuple<int32, uchar*, uint64> data_t;
	typedef std::list<data_t> list_t;
	list_t l;
	bool completed;
	mutable std::mutex mtx;
	std::condition_variable cv_pop;
  public:
	void init() {
		completed = false;
	}
	CBigBinPartQueue() {
		init();
	}

	void push(int32 bin_id, uchar* data, uint64 size) {
		std::lock_guard<std::mutex> lck(mtx);
		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, data, size));

		if (was_empty)
			cv_pop.notify_all();
	}

	bool pop(int32& bin_id, uchar* &data, uint64& size) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] {return !l.empty() || completed; });

		if (completed && l.empty())
			return false;

		bin_id = std::get<0>(l.front());
		data = std::get<1>(l.front());
		size = std::get<2>(l.front());
		l.pop_front();
		return true;
	}

	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		completed = true;
		cv_pop.notify_all();
	}
};

class CBigBinKXmersQueue {
	typedef std::tuple<int32, uchar*, uint64> data_t;
	typedef std::list<data_t> list_t;
	list_t l;
	uint32 n_writers;
	mutable std::mutex mtx;
	std::condition_variable cv_pop;

	uint32 n_waiters;
	int32 current_id;
	std::condition_variable cv_push;
  public:
	CBigBinKXmersQueue(uint32 _n_writers) {
		n_waiters = 0;
		current_id = -1; //means queue is not initialized
		n_writers = _n_writers;
	}

	void push(int32 bin_id, uchar* data, uint64 size) {
		std::unique_lock<std::mutex> lck(mtx);
		++n_waiters;

		if (current_id == -1)
			current_id = bin_id;

		cv_push.wait(lck, [this, bin_id] {return bin_id == current_id || n_waiters == n_writers; });

		if (n_waiters == n_writers) {
			current_id = bin_id;
			cv_push.notify_all();
		}

		--n_waiters;
		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, data, size));

		if(was_empty)
			cv_pop.notify_all();
	}

	bool pop(int32& bin_id, uchar* &data, uint64& size) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] {return !l.empty() || !n_writers; });

		if (l.empty() && !n_writers)
			return false;

		bin_id = std::get<0>(l.front());
		data = std::get<1>(l.front());
		size = std::get<2>(l.front());
		l.pop_front();
		return true;
	}

	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		--n_writers;

		if (!n_writers)
			cv_pop.notify_all();

		cv_push.notify_all();
	}

	~CBigBinKXmersQueue() {
	}

};

class CBigBinSortedPartQueue {
	//bin_id, sub_bin_id,suff_buff, suff_buff_size, lut, lut_size, last_one_in_bin
	typedef std::tuple<int32, int32, uchar*, uint64, uint64*, uint64, bool> data_t;
	typedef std::list<data_t> list_t;
	list_t l;
	uint32 n_writers;
	mutable std::mutex mtx;
	std::condition_variable cv_pop;
  public:
	CBigBinSortedPartQueue(uint32 _n_writers) {
		n_writers = _n_writers;
	}
	void push(int32 bin_id, int32 sub_bin_id, uchar* suff_buff, uint64 suff_buff_size, uint64* lut,
			  uint64 lut_size, bool last_one_in_sub_bin) {
		std::lock_guard<std::mutex> lck(mtx);
		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, sub_bin_id, suff_buff, suff_buff_size, lut, lut_size,
									last_one_in_sub_bin));

		if (was_empty)
			cv_pop.notify_all();
	}
	bool pop(int32& bin_id, int32& sub_bin_id, uchar* &suff_buff, uint64& suff_buff_size, uint64* &lut,
			 uint64 &lut_size, bool &last_one_in_sub_bin) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] {return !n_writers || !l.empty(); });

		if (!n_writers && l.empty())
			return false;

		bin_id			= std::get<0>(l.front());
		sub_bin_id		= std::get<1>(l.front());
		suff_buff		= std::get<2>(l.front());
		suff_buff_size	= std::get<3>(l.front());
		lut				= std::get<4>(l.front());
		lut_size		= std::get<5>(l.front());
		last_one_in_sub_bin	= std::get<6>(l.front());
		l.pop_front();
		return true;
	}
	void mark_completed() {
		--n_writers;

		if (!n_writers)
			cv_pop.notify_all();
	}
};

class CBigBinKmerPartQueue {
	typedef std::tuple<int32, uchar*, uint64, uchar*, uint64, uint64, uint64, uint64, uint64, bool> data_t;
	typedef std::list<data_t> list_t;
	list_t l;
	uint32 n_writers;
	mutable std::mutex mtx;
	std::condition_variable cv_pop;
	std::condition_variable cv_push;
	int32 curr_id;
	bool allow_next;
  public:
	CBigBinKmerPartQueue(uint32 _n_writers) {
		n_writers = _n_writers;
		allow_next = true;
	}
	void push(int32 bin_id, uchar* suff_buff, uint64 suff_buff_size, uchar* lut, uint64 lut_size, uint64 n_unique,
			  uint64 n_cutoff_min, uint64 n_cutoff_max, uint64 n_total, bool last_in_bin) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_push.wait(lck, [this, bin_id, lut_size] {return curr_id == bin_id || allow_next; });
		allow_next = false;

		if (last_in_bin) {
			allow_next = true;
		}

		curr_id = bin_id;
		bool was_empty = l.empty();
		l.push_back(std::make_tuple(bin_id, suff_buff, suff_buff_size, lut, lut_size, n_unique, n_cutoff_min,
									n_cutoff_max, n_total, last_in_bin));

		if (was_empty)
			cv_pop.notify_all();

		if (allow_next)
			cv_push.notify_all();
	}
	bool pop(int32& bin_id, uchar* &suff_buff, uint64& suff_buff_size, uchar* &lut, uint64& lut_size,
			 uint64 &n_unique, uint64 &n_cutoff_min, uint64 &n_cutoff_max, uint64 &n_total, bool& last_in_bin) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] {return !l.empty() || !n_writers; });

		if (!n_writers && l.empty())
			return false;

		bin_id = std::get<0>(l.front());
		suff_buff = std::get<1>(l.front());
		suff_buff_size = std::get<2>(l.front());
		lut = std::get<3>(l.front());
		lut_size = std::get<4>(l.front());
		n_unique = std::get<5>(l.front());
		n_cutoff_min = std::get<6>(l.front());
		n_cutoff_max = std::get<7>(l.front());
		n_total = std::get<8>(l.front());
		last_in_bin = std::get<9>(l.front());
		l.pop_front();
		return true;
	}

	void mark_completed() {
		std::lock_guard<std::mutex> lck(mtx);
		--n_writers;

		if(!n_writers)
			cv_pop.notify_all();
	}
};


class CBigBinDesc {
	//lut_prefix_len, n_kmers, tmp_file_handle, std::string file_name, file_size
	typedef std::tuple<uint32, uint32, FILE*, std::string, uint64> elem_t;
	typedef std::map<int32, std::pair<int32, std::map<int32, elem_t>>> data_t;
	mutable std::mutex mtx;
	data_t m;
	int32 curr_id;
  public:
	CBigBinDesc() {
		curr_id = -1;
	}
	void push(int32 bin_id, int32 sub_bin_id, uint32 lut_prefix_len, uint32 n_kmers, FILE* file, std::string desc,
			  uint64 file_size) {
		std::lock_guard<std::mutex> lck(mtx);
		auto bin = m.find(bin_id);

		if (bin == m.end()) {
			m[bin_id].first = -1;
			m[bin_id].second[sub_bin_id] = std::make_tuple(lut_prefix_len, n_kmers, file, desc, file_size);
		} else {
			auto sub_bin = bin->second.second.find(sub_bin_id);

			if (sub_bin == bin->second.second.end()) {
				m[bin_id].second[sub_bin_id] = std::make_tuple(lut_prefix_len, n_kmers, file, desc, file_size);
			} else {
				if(lut_prefix_len)
					std::get<0>(sub_bin->second) = lut_prefix_len;

				std::get<1>(sub_bin->second) += n_kmers;

				if (file) {
					std::get<2>(sub_bin->second) = file;
					std::get<3>(sub_bin->second) = desc;
				}

				std::get<4>(sub_bin->second) += file_size;
			}
		}
	}

	bool get_n_sub_bins(int32 bin_id, uint32& size) {
		std::lock_guard<std::mutex> lck(mtx);
		auto e = m.find(bin_id);

		if (e == m.end())
			return false;

		size = (uint32)e->second.second.size();
		return true;
	}

	bool next_bin(int32& bin_id, uint32& size) {
		std::lock_guard<std::mutex> lck(mtx);

		if (m.empty())
			return false;

		if (curr_id == -1) {
			curr_id = bin_id = m.begin()->first;
			size = (uint32)m.begin()->second.second.size();
		} else {
			auto e = m.find(curr_id);
			e++;

			if (e == m.end())
				return false;

			curr_id = bin_id = e->first;
			size = (uint32)e->second.second.size();
		}

		return true;
	}

	void reset_reading() {
		std::lock_guard<std::mutex> lck(mtx);
		curr_id = -1;

		for (auto& e : m)
			e.second.first = -1;
	}

	bool next_sub_bin(int32 bin_id, int32& sub_bin_id, uint32& lut_prefix_len, uint32& n_kmers, FILE* &file,
					  std::string& desc, uint64& file_size) {
		std::lock_guard<std::mutex> lck(mtx);
		auto& sub_bin = m.find(bin_id)->second;
		int32 curr_sub_bin_id = sub_bin.first;
		std::map<int32, elem_t>::iterator e;

		if (curr_sub_bin_id == -1)
			e = sub_bin.second.begin();
		else {
			e = sub_bin.second.find(curr_sub_bin_id);
			++e;

			if (e == sub_bin.second.end())
				return false;
		}

		sub_bin_id = sub_bin.first = e->first;
		lut_prefix_len = std::get<0>(e->second);
		n_kmers = std::get<1>(e->second);
		file = std::get<2>(e->second);
		desc = std::get<3>(e->second);
		file_size = std::get<4>(e->second);
		return true;
	}
};

class CCompletedBinsCollector {
	std::list<int32> l;
	mutable std::mutex mtx;
	std::condition_variable cv_pop;
	uint32 n_writers;
  public:
	CCompletedBinsCollector(uint32 _n_writers) {
		n_writers = _n_writers;
	}
	void push(int32 bin_id) {
		std::lock_guard<std::mutex> lck(mtx);
		bool was_empty = l.empty();
		l.push_back(bin_id);

		if (was_empty)
			cv_pop.notify_all();
	}

	bool pop(int32& bin_id) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] {return !n_writers || !l.empty(); });

		if (!n_writers && l.empty())
			return false;

		bin_id = l.front();
		l.pop_front();
		return true;
	}
	void mark_completed() {
		--n_writers;

		if (!n_writers)
			cv_pop.notify_all();
	}
};


class CDiskLogger {
	uint64 current;
	uint64 max;
	mutable std::mutex mtx;

  public:
	CDiskLogger() {
		current = max = 0;
	}
	void log_write(uint64 _size) {
		std::lock_guard<std::mutex> lck(mtx);
		current += _size;

		if (current > max)
			max = current;
	}
	void log_remove(uint64 _size) {
		std::lock_guard<std::mutex> lck(mtx);
		current -= _size;
	}
	uint64 get_max() {
		std::lock_guard<std::mutex> lck(mtx);
		return max;
	}
	uint64 get_current() {
		std::lock_guard<std::mutex> lck(mtx);
		return current;
	}

};
#endif

// ***** EOF
