/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _THREAD_WATCH_H
#define _THREAD_WATCH_H

#ifdef WIN32

class CThreadWatchImpl;

// **********************************************************
class CThreadWatch
{	
	CThreadWatchImpl* pimpl;
public:
	CThreadWatch();
	void StartTimer();
	void StopTimer();
	double GetElapsedTime();
	~CThreadWatch();
};

#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

typedef timeval thread_watch_t;

class CThreadWatch
{
	thread_watch_t start_kernel, start_user;
	thread_watch_t stop_kernel, stop_user;

public:
	CThreadWatch();
	void StartTimer();
	void StopTimer();
	double GetElapsedTime();
};

#endif

#endif
// ***** EOF