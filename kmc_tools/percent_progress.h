/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _PERCENT_PROGRESS_H
#define _PERCENT_PROGRESS_H

#include "defs.h"
#include <vector>
#include <string>
//************************************************************************************************************
// CPercentProgress - class to display progress of reading inputs
//************************************************************************************************************
class CPercentProgress
{
	bool ignore_rest = false;
	bool hide_progress = false;
	struct CDisplayItem
	{
		std::string name;
		uint64 cur_val = 0;
		uint64 max_val;
		uint32 cur_percent = 0;

		uint32 to_next_update;
		uint32 to_next_update_pattern;
	public:
		CDisplayItem(const std::string name, uint64 max_val);
	};
	std::vector<CDisplayItem> items;
	void display();
public:
	uint32 RegisterItem(const std::string& name, uint64 max_value);
	uint32 RegisterItem(uint64 max_value);
	void UpdateItem(uint32 id);
	void Complete(uint32 id);
	void UpdateItem(uint32 id, uint32 offset);
	void Hide(){ hide_progress = true; }
};

#endif

// ***** EOF