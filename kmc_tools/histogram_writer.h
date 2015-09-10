/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _HISTOGRAM_WRITER_H
#define _HISTOGRAM_WRITER_H

#include "defs.h"
#include "config.h"
#include <vector>
#include <fstream>

template<typename KMCDB> class CHistogramWriter
{
	KMCDB& kmcdb;
	COutputDesc& output_desc;
	std::vector<uint32> counters;

public:
	CHistogramWriter(KMCDB& kmcdb) :kmcdb(kmcdb), output_desc(CConfig::GetInstance().output_desc)
	{

	}
	bool Process()
	{
		counters.resize(output_desc.cutoff_max + 1);
		uint32 counter;
		while (kmcdb.NextCounter(counter))
		{
			if (counter >= output_desc.cutoff_min && counter <= output_desc.cutoff_max)
				counters[counter]++;
		}
		std::ofstream file(output_desc.file_src);
		if (!file)
		{
			std::cout << "Error: cannot open file: " << output_desc.file_src << "\n";
			exit(1);
		}
		for (uint32 i = output_desc.cutoff_min; i <= output_desc.cutoff_max; ++i)
		{
			file << i << "\t" << counters[i] << "\n";
		}
		file.close();
		return true;
	}
};

#endif