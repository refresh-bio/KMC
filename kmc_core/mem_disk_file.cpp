/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/

#include "mem_disk_file.h"
#include "critical_error_handler.h"
#include <sstream>
using namespace std;

//----------------------------------------------------------------------------------
// Constructor 
CMemDiskFile::CMemDiskFile(bool _memory_mode, bool _reopen_each_time) :
	memory_mode(_memory_mode),
	reopen_each_time(_reopen_each_time)
{
}

//----------------------------------------------------------------------------------
void CMemDiskFile::Open(const string& f_name)
{	
	if(memory_mode)
	{

	}
	else
	{
		file = fopen(f_name.c_str(), "wb+");

		if (!file)
		{
			std::ostringstream ostr;
			ostr << "Error: Cannot open temporary file " << f_name;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		if (reopen_each_time)
		{
			fclose(file);
			file = nullptr;
		}
		else
			setbuf(file, nullptr);
	}
	name = f_name;
}

//----------------------------------------------------------------------------------
void CMemDiskFile::Rewind()
{
	if (memory_mode)
	{

	}
	else
	{
		if (reopen_each_time)
			read_pos = 0;
		else
			rewind(file);
	}
}

//----------------------------------------------------------------------------------
int CMemDiskFile::Close()
{
	if(memory_mode)
	{
		for(auto& p : container)
		{
			delete[] p.first;
		}
		container.clear();
		return 0;
	}
	else
	{
		if (file)
		{
			auto ret = fclose(file);
			file = nullptr;
			return ret;
		}
		else
			return 0;
	}
}
//----------------------------------------------------------------------------------
void CMemDiskFile::Remove()
{
	if (!memory_mode)
		remove(name.c_str());
}
//----------------------------------------------------------------------------------
size_t CMemDiskFile::Read(uchar * ptr, size_t size, size_t count)
{
	if(memory_mode)
	{
		uint64 pos = 0;
		for(auto& p : container)
		{
			memcpy(ptr + pos, p.first, p.second);
			pos += p.second;
			delete[] p.first;
		}
		container.clear();
		return pos;
	}
	else
	{
		if (reopen_each_time)
		{
			file = fopen(name.c_str(), "rb");
			if (!file)
			{
				std::ostringstream ostr;
				ostr << "Error: Cannot open temporary file " << name;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			setbuf(file, nullptr);
			my_fseek(file, read_pos, SEEK_SET);
			auto res = fread(ptr, size, count, file);
			read_pos = my_ftell(file);
			fclose(file);
			file = nullptr;
			return res;
		}
		else
			return fread(ptr, size, count, file);
	}
}

//----------------------------------------------------------------------------------
size_t CMemDiskFile::Write(const uchar * ptr, size_t size, size_t count)
{
	if(memory_mode)
	{
		uchar *buf = new uchar[size * count];
		memcpy(buf, ptr, size * count);
		container.push_back(make_pair(buf, size * count));
		return size * count;
	}
	else
	{
		if (reopen_each_time)
		{
			file = fopen(name.c_str(), "ab");

			if (!file)
			{
				std::ostringstream ostr;
				ostr << "Error: Cannot open temporary file " << name;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			setbuf(file, nullptr);
			auto res = fwrite(ptr, size, count, file);
			fclose(file);
			file = nullptr;
			return res;
		}
		else
			return fwrite(ptr, size, count, file);
	}
}

//----------------------------------------------------------------------------------
CMemDiskFile::~CMemDiskFile()
{
	Close();
	Remove();
}

// ***** EOF
