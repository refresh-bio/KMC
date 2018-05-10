/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include "percent_progress.h"
#include <iostream>
#include <string>
using namespace std;


/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
uint32 CPercentProgress::RegisterItem(const std::string& name, uint64 max_value)
{
	items.emplace_back(name, max_value);
	return static_cast<uint32>(items.size() - 1);
}

/*****************************************************************************************************************************/
uint32 CPercentProgress::RegisterItem(uint64 max_value)
{
	items.emplace_back("in" + std::to_string(items.size() + 1), max_value);
	display();
	return static_cast<uint32>(items.size() - 1);
}

/*****************************************************************************************************************************/
void CPercentProgress::UpdateItem(uint32 id)
{
	--items[id].to_next_update;
	if (!items[id].to_next_update)
	{
		items[id].to_next_update = items[id].to_next_update_pattern;
		UpdateItem(id, items[id].to_next_update_pattern);
	}
}

/*****************************************************************************************************************************/
void CPercentProgress::UpdateItem(uint32 id, uint32 offset)
{
	items[id].cur_val += offset;
	uint32 prev = items[id].cur_percent;
	if (items[id].max_val)
		items[id].cur_percent = static_cast<uint32>((items[id].cur_val * 100) / items[id].max_val);
	else
		items[id].cur_percent = 100;
	if (prev != items[id].cur_percent)
		display();
}

/*****************************************************************************************************************************/
void CPercentProgress::Complete(uint32 id)
{
	if (items[id].cur_percent != 100)
	{
		items[id].cur_percent = 100;
		display();		
	}
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
void CPercentProgress::display()
{
	if (hide_progress || ignore_rest)
		return;
	std::cerr << "\r";
	bool finished = true;
	for (auto& item : items)
	{
		std::cerr << item.name << ": " << item.cur_percent << "% ";
		if (item.cur_percent != 100)
			finished = false;
	}

	if (finished)
	{
		cerr << "\n";
		ignore_rest = true;
	}
	std::cerr.flush();
}

/*****************************************************************************************************************************/
CPercentProgress::CDisplayItem::CDisplayItem(const std::string name, uint64 max_val) : name(name), max_val(max_val)
{
	to_next_update_pattern = (uint32)MAX(1, max_val / 100);
	to_next_update = to_next_update_pattern;
}

// ***** EOF