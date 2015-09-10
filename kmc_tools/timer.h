/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>
class CTimer
{
	using time_p = std::chrono::system_clock::time_point;
	time_p _start, _end;
public:
	void start()
	{
		_start = std::chrono::high_resolution_clock::now();
	}
	double get_time()
	{
		auto time = std::chrono::high_resolution_clock::now() - _start;
		return static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(time).count());
	}
};



#endif 
