/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#include "stdafx.h"

#ifdef WIN32
#include <windows.h>
#endif

#include <cstdio> // NULL
#include "thread_watch.h"


#ifdef WIN32

typedef struct
{
	ULARGE_INTEGER start;
	ULARGE_INTEGER stop;
} thread_watch_t;

class CThreadWatchImpl
{
	thread_watch_t timer_kernel, timer_user;
	LARGE_INTEGER frequency;
	double LIToSecs(LARGE_INTEGER & L);
public:
	CThreadWatchImpl();
	void StartTimer();
	void StopTimer();
	double GetElapsedTime();
};

double CThreadWatchImpl::LIToSecs(LARGE_INTEGER & L)
{
	return ((double)L.QuadPart / (double)frequency.QuadPart);
}

// **********************************************************
CThreadWatchImpl::CThreadWatchImpl()
{
	timer_kernel.start.QuadPart = 0;
	timer_kernel.stop.QuadPart = 0;
	timer_user.start.QuadPart = 0;
	timer_user.stop.QuadPart = 0;
	//	QueryPerformanceFrequency( &frequency );
	//	frequency = 1;		// 100ns
}

// **********************************************************
void CThreadWatchImpl::StartTimer()
{
	FILETIME CreationTime, ExitTime, KernelTime, UserTime;
	GetThreadTimes(GetCurrentThread(), &CreationTime, &ExitTime, &KernelTime, &UserTime);

	timer_kernel.start.LowPart = KernelTime.dwLowDateTime;
	timer_kernel.start.HighPart = KernelTime.dwHighDateTime;
	timer_user.start.LowPart = UserTime.dwLowDateTime;
	timer_user.start.HighPart = UserTime.dwHighDateTime;
}

// **********************************************************
void CThreadWatchImpl::StopTimer()
{
	FILETIME CreationTime, ExitTime, KernelTime, UserTime;
	GetThreadTimes(GetCurrentThread(), &CreationTime, &ExitTime, &KernelTime, &UserTime);
	//    QueryPerformanceCounter(&timer.stop);
	timer_kernel.stop.LowPart = KernelTime.dwLowDateTime;
	timer_kernel.stop.HighPart = KernelTime.dwHighDateTime;
	timer_user.stop.LowPart = UserTime.dwLowDateTime;
	timer_user.stop.HighPart = UserTime.dwHighDateTime;
}

// **********************************************************
double CThreadWatchImpl::GetElapsedTime()
{
	/*	LARGE_INTEGER time;
	time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
	return LIToSecs( time) ;*/
	LARGE_INTEGER time;

	time.QuadPart = (timer_kernel.stop.QuadPart - timer_kernel.start.QuadPart);
	time.QuadPart += (timer_user.stop.QuadPart - timer_user.start.QuadPart);

	return (double)time.QuadPart / 1e7;			// 100ns clock
}



CThreadWatch::CThreadWatch()
{
	pimpl = new CThreadWatchImpl;
}
void CThreadWatch::StartTimer()
{
	pimpl->StartTimer();
}
void CThreadWatch::StopTimer()
{
	pimpl->StopTimer();
}
double CThreadWatch::GetElapsedTime()
{
	return pimpl->GetElapsedTime();
}
CThreadWatch::~CThreadWatch()
{
	delete pimpl;
}
#else

// **********************************************************

// **********************************************************
// CThreadWatch
// **********************************************************
/*double CThreadWatch::LIToSecs( LARGE_INTEGER & L)
{
return 1.0;
}*/

// **********************************************************
CThreadWatch::CThreadWatch()
{
}

// **********************************************************
void CThreadWatch::StartTimer()
{
	rusage usage;
	getrusage(RUSAGE_THREAD, &usage);
	start_user = usage.ru_utime;
	start_kernel = usage.ru_stime;
}

// **********************************************************
void CThreadWatch::StopTimer()
{
	rusage usage;
	getrusage(RUSAGE_THREAD, &usage);
	stop_user = usage.ru_utime;
	stop_kernel = usage.ru_stime;
}

// **********************************************************
double CThreadWatch::GetElapsedTime()
{
	double ret = 0.0;

	ret += stop_user.tv_sec + stop_user.tv_usec / 1000000.0;
	ret += stop_kernel.tv_sec + stop_kernel.tv_usec / 1000000.0;
	ret -= start_user.tv_sec + start_user.tv_usec / 1000000.0;
	ret -= start_kernel.tv_sec + start_kernel.tv_usec / 1000000.0;

	return ret;
}

#endif