/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _QUEUES_H_
#define _QUEUES_H_

#include "defs.h"
#include "bundle.h"
#include <mutex>
#include <vector>
#include <condition_variable>
#include <list>
#include <queue>



class CSufWriteQueue
{
	uint32 buf_size;
	uint32 max_inside;

	using elem_t = std::pair<uchar*, uint32>;
	std::list<elem_t> content;

	mutable std::mutex mtx;
	uint32 n_writers;
	std::condition_variable cv_pop, cv_push;
public:
	void init(uint32 _buf_size, uint32 _max_inside)
	{
		buf_size = _buf_size;
		max_inside = _max_inside;
		n_writers = 1;
	}

	void push(uchar* &buf, uint32 size)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_push.wait(lck, [this]{return content.size() < max_inside; });

		bool was_empty = content.empty();

		content.push_back(std::make_pair(buf, size));

		buf = new uchar[buf_size];

		if (was_empty)
			cv_pop.notify_all();
	}

	bool pop(uchar* &buf, uint32& size)
	{
		std::unique_lock<std::mutex> lck(mtx);
			cv_pop.wait(lck, [this]{return !content.empty() || !n_writers; });
		if (!n_writers && content.empty())
			return false;

		bool was_full = max_inside == content.size();

		buf = content.front().first;
		size = content.front().second;
		content.pop_front();

		if (was_full)
			cv_push.notify_all();

		return true;
	}


	void mark_completed()
	{
		std::lock_guard<std::mutex> lck(mtx);
		--n_writers;
		if (!n_writers)
			cv_pop.notify_all();
	}
};


template<unsigned SIZE> class CCircularQueue
{
	std::vector<CBundleData<SIZE>> buff;
	bool full, is_completed;
	int start, end;
	mutable std::mutex mtx;

	std::condition_variable cv_push;
	std::condition_variable cv_pop;
	bool forced_to_finish = false;

public:
	CCircularQueue(int size, uint32 bundle_size) : full(false), is_completed(false), start(0), end(0)
	{
		buff.reserve(size);
		for (int i = 0; i < size; ++i)
			buff.emplace_back(bundle_size);
	}
	CCircularQueue(int size) : buff(size), full(false), is_completed(false), start(0), end(0)
	{

	}

	bool push(CBundleData<SIZE>& bundle_data)
	{
		std::unique_lock<std::mutex> lck(mtx);

		cv_push.wait(lck, [this]{return !full || forced_to_finish; });
		
		if (forced_to_finish)
		{
			return false;
		}

		bool was_empty = start == end;

		std::swap(buff[end], bundle_data);
		bundle_data.Clear();
		end = (end + 1) % buff.size();

		if (end == start)
			full = true;

		if (was_empty)
			cv_pop.notify_all();

		return true;
	}

	bool pop(CBundleData<SIZE>& bundle_data)
	{
		std::unique_lock<std::mutex> lck(mtx);		
		cv_pop.wait(lck, [this]{ return start != end || full || is_completed || forced_to_finish; });

		if (forced_to_finish)
			return false;

		if (is_completed && !full && start == end)
			return false;

		bool was_full = full;		
		std::swap(buff[start], bundle_data);
		buff[start].Clear();
		start = (start + 1) % buff.size();
		full = false;
		if (was_full)
			cv_push.notify_all();
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

class CInputFilesQueue {
	typedef std::string elem_t;
	typedef std::queue<elem_t, std::list<elem_t>> queue_t;

	queue_t q;

	mutable std::mutex mtx;								// The mutex to synchronise on

public:
	CInputFilesQueue(const std::vector<std::string> &file_names) {
		std::unique_lock<std::mutex> lck(mtx);

		for (auto p = file_names.cbegin(); p != file_names.cend(); ++p)
			q.push(*p);

	};

	bool pop(std::string &file_name) {
		std::lock_guard<std::mutex> lck(mtx);

		if (q.empty())
			return false;

		file_name = q.front();
		q.pop();

		return true;
	}
};

class CMemoryPool {
	int64 total_size;
	int64 part_size;
	int64 n_parts_total;
	int64 n_parts_free;

	uchar *buffer, *raw_buffer;
	uint32 *stack;

	mutable std::mutex mtx;							// The mutex to synchronise on
	std::condition_variable cv;						// The condition to wait for

public:
	CMemoryPool(int64 _total_size, int64 _part_size) {
		raw_buffer = NULL;
		buffer = NULL;
		stack = NULL;
		prepare(_total_size, _part_size);
	}
	~CMemoryPool() {
		release();
	}

	void prepare(int64 _total_size, int64 _part_size) {
		release();

		n_parts_total = _total_size / _part_size;
		part_size = (_part_size + 15) / 16 * 16;			// to allow mapping pointer to int*
		n_parts_free = n_parts_total;

		total_size = n_parts_total * part_size;

		raw_buffer = new uchar[total_size + 64];
		buffer = raw_buffer;
		while (((uint64)buffer) % 64)
			buffer++;

		stack = new uint32[n_parts_total];
		for (uint32 i = 0; i < n_parts_total; ++i)
			stack[i] = i;
	}

	void release(void) {
		if (raw_buffer)
			delete[] raw_buffer;
		raw_buffer = NULL;
		buffer = NULL;

		if (stack)
			delete[] stack;
		stack = NULL;
	}

	// Allocate memory buffer - uchar*
	void reserve(uchar* &part)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0; });

		part = buffer + stack[--n_parts_free] * part_size;
	}
	// Allocate memory buffer - char*
	void reserve(char* &part)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0; });

		part = (char*)(buffer + stack[--n_parts_free] * part_size);
	}
	// Allocate memory buffer - uint32*
	void reserve(uint32* &part)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0; });

		part = (uint32*)(buffer + stack[--n_parts_free] * part_size);
	}
	// Allocate memory buffer - uint64*
	void reserve(uint64* &part)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0; });

		part = (uint64*)(buffer + stack[--n_parts_free] * part_size);
	}
	// Allocate memory buffer - double*
	void reserve(double* &part)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0; });

		part = (double*)(buffer + stack[--n_parts_free] * part_size);
	}

	// Deallocate memory buffer - uchar*
	void free(uchar* part)
	{
		std::lock_guard<std::mutex> lck(mtx);

		stack[n_parts_free++] = (uint32)((part - buffer) / part_size);

		cv.notify_all();
	}
	// Deallocate memory buffer - char*
	void free(char* part)
	{
		std::lock_guard<std::mutex> lck(mtx);

		stack[n_parts_free++] = (uint32)(((uchar*)part - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - uint32*
	void free(uint32* part)
	{
		std::lock_guard<std::mutex> lck(mtx);

		stack[n_parts_free++] = (uint32)((((uchar *)part) - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - uint64*
	void free(uint64* part)
	{
		std::lock_guard<std::mutex> lck(mtx);

		stack[n_parts_free++] = (uint32)((((uchar *)part) - buffer) / part_size);
		cv.notify_all();
	}
	// Deallocate memory buffer - double*
	void free(double* part)
	{
		std::lock_guard<std::mutex> lck(mtx);

		stack[n_parts_free++] = (uint32)((((uchar *)part) - buffer) / part_size);
		cv.notify_all();
	}
};

class CPartQueue {
	typedef std::pair<uchar *, uint64> elem_t;
	typedef std::queue<elem_t, std::list<elem_t>> queue_t;

	queue_t q;
	bool is_completed;
	int n_readers;

	mutable std::mutex mtx;								// The mutex to synchronise on
	std::condition_variable cv_queue_empty;

public:
	CPartQueue(int _n_readers) {
		std::unique_lock<std::mutex> lck(mtx);
		is_completed = false;
		n_readers = _n_readers;
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
		if (!n_readers)
			cv_queue_empty.notify_all();
	}
	void push(uchar *part, uint64 size) {
		std::unique_lock<std::mutex> lck(mtx);

		bool was_empty = q.empty();
		q.push(std::make_pair(part, size));

		if (was_empty)
			cv_queue_empty.notify_all();
	}
	bool pop(uchar *&part, uint64 &size) {
		std::unique_lock<std::mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this]{return !this->q.empty() || !this->n_readers; });

		if (q.empty())
			return false;

		std::tie(part, size) = q.front();
		q.pop();

		return true;
	}
};


struct CFilteringQueues
{
	CInputFilesQueue *input_files_queue;
	CPartQueue *input_part_queue, *filtered_part_queue;
	CMemoryPool *pmm_fastq_reader;
	CMemoryPool *pmm_fastq_filter;
};


#endif

