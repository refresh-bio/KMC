/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include "fastq_writer.h"
#include <iostream>
using namespace std;


/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
CFastqWriter::CFastqWriter(CFilteringParams& Params, CFilteringQueues& Queues)
{
	output_src			= Params.output_src;
	filtered_part_queue = Queues.filtered_part_queue;
	pmm_fastq_filter	= Queues.pmm_fastq_filter;
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/
void CFastqWriter::Process()
{
	uchar* part;
	uint64 size;
	FILE* f = fopen(output_src.c_str(), "wb");
	if (!f)
	{
		cerr << "cannot open file :" << output_src;
		exit(1);
	}
	while (filtered_part_queue->pop(part, size))
	{
		if (fwrite(part, 1, size, f) != size)
		{
			cerr << "Error while writing to " << output_src << "\n";
			exit(1);
		}
		pmm_fastq_filter->free(part);
	}
	fclose(f);
}

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/
CWFastqWriter::CWFastqWriter(CFilteringParams& Params, CFilteringQueues& Queues)
	:writer(Params, Queues)
{

}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/
void CWFastqWriter::operator()()
{
	writer.Process();
}

// ***** EOF