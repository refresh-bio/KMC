/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _MEM_DISK_FILE_H
#define _MEM_DISK_FILE_H

#include "defs.h"
#include <string>
#include <stdio.h>
#include <vector>
using namespace std;


//************************************************************************************************************
// CMemDiskFile - wrapper for FILE* or memory equivalent
//************************************************************************************************************
class CMemDiskFile
{
	bool memory_mode;
	FILE* file;
	typedef pair<uchar*, uint64> elem_t;//buf,size
	typedef vector<elem_t> container_t;

	container_t container;
	string name;
public:
	CMemDiskFile(bool _memory_mode);
	void Open(const string& f_name);
	void Rewind();
	int Close();
	size_t Read(uchar * ptr, size_t size, size_t count);
	size_t Write(const uchar * ptr, size_t size, size_t count);
	void Remove();
};

#endif

// ***** EOF