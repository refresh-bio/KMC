/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.0.0
Date   : 2019-05-19
*/

#ifndef _PERCENT_PROGRESS_H
#define _PERCENT_PROGRESS_H
#include "defs.h"
#include <string>
#include <iostream>
class CPercentProgress
{
	uint64 curr_val;
	uint64 max_val;	
	int32 curr_percent;
	std::string label;
	bool show_progress;
public:
	CPercentProgress(const std::string& label, bool show_progress)
		:
		curr_val(0),
		max_val(0),
		curr_percent(-1),
		label(label),
		show_progress(show_progress)
	{

	}
	void SetMaxVal(uint64 _max_val)
	{
		max_val = _max_val;
	}

	void NotifyProgress(uint64 val)
	{
		curr_val += val;
		int32 new_percent = 0;
		if (max_val)
			new_percent = static_cast<int32>((curr_val * 100) / max_val);
		else
			new_percent = 100;
		if (new_percent > curr_percent)
		{
			curr_percent = new_percent;
			if (show_progress)
			{
				std::cerr << "\r" << label << curr_percent << "%";
				std::cerr.flush();
			}
		}
	}
};


#endif

// ***** EOF