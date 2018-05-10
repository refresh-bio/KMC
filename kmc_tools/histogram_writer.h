/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _HISTOGRAM_WRITER_H
#define _HISTOGRAM_WRITER_H

#include "defs.h"
#include "config.h"
#include <vector>
#include <fstream>


class CHistogramWriterBase
{
private:
	std::string& file_src;
	uint32 cutoff_max;
	uint32 cutoff_min;
	std::vector<uint64> counters;
protected:
	void Init()
	{
		counters.resize(cutoff_max + 1);
	}
	void ProcessCounter(uint32 counter)
	{
		if (counter >= cutoff_min && counter <= cutoff_max)
			counters[counter]++;
	}
	void Finish()
	{
		std::ofstream file(file_src);
		if (!file)
		{
			std::cerr << "Error: cannot open file: " << file_src << "\n";
			exit(1);
		}
		for (uint32 i = cutoff_min; i <= cutoff_max; ++i)
		{
			file << i << "\t" << counters[i] << "\n";
		}
		file.close();
	}
	CHistogramWriterBase(std::string& file_src, uint32 cutoff_max, uint32 cutoff_min):
		file_src(file_src),
		cutoff_max(cutoff_max),
		cutoff_min(cutoff_min)
	{
	}
};

template<typename KMCDB> class CHistogramWriter : public CHistogramWriterBase
{
	KMCDB& kmcdb;
		
public:
	CHistogramWriter(KMCDB& kmcdb) :
		CHistogramWriterBase(CConfig::GetInstance().output_desc.file_src, CConfig::GetInstance().output_desc.cutoff_max, CConfig::GetInstance().output_desc.cutoff_min),
		kmcdb(kmcdb)
	{

	}
	bool Process()
	{
		Init();
		uint32 counter;
		while (kmcdb.NextCounter(counter))
		{
			ProcessCounter(counter);
		}
		Finish();
		return true;
	}
};


class CHistogramWriterForTransform : public CHistogramWriterBase
{
public:
	CHistogramWriterForTransform(CTransformOutputDesc& output_desc) :
		CHistogramWriterBase(output_desc.file_src, output_desc.cutoff_max, output_desc.cutoff_min)
	{

	}

	void Init()
	{
		CHistogramWriterBase::Init();
	}

	void PutCounter(uint32 counter)
	{
		ProcessCounter(counter);
	}

	void Finish()
	{
		CHistogramWriterBase::Finish();
	}
};

#endif