/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.2
  Date   : 2023-03-09
*/

#ifndef _TMP_FILES_OWNER_H
#define _TMP_FILES_OWNER_H

#include "mem_disk_file.h"
#include <memory>

class CTmpFilesOwner
{
	std::vector<std::unique_ptr<CMemDiskFile>> files;
	bool memory_mode;
public:
	CTmpFilesOwner(uint32_t n_bins, bool memory_mode) :
		files(n_bins),
		memory_mode(memory_mode)
	{

	}
	void CreateInstances()
	{
		for(auto& uptr : files)
			uptr = std::make_unique<CMemDiskFile>(memory_mode);
	}
	
	CMemDiskFile* Get(uint32_t index)
	{
		return files[index].get();
	}

	void Release()
	{
		files.clear();
	}
};

#endif

// ***** EOF
