/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.2.3
Date   : 2023-12-08
*/

#ifndef _KFF_KMC2_UTILS_H
#define _KFF_KMC2_UTILS_H

#include "defs.h"
#include <mutex>
#include <thread>

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
		cv.wait(lck, [this] {return !n_tasks; });
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
		cv.wait(lck, [this] {return empty; });
		empty = false;
		cv.notify_all();
	}

	bool pop(CParentSubthreadDesc<SIZE>& _desc)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this] {return !empty || completed; });
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

#endif