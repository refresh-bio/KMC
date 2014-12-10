/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 2.0
  Date   : 2014-07-04
*/

#include "mmer.h"
#include "kmc_file.h"
#include <iostream>


uint64 CKMCFile::part_size = 1 << 25;

// ----------------------------------------------------------------------------------
// Open files *.kmc_pre & *.kmc_suf, read them to RAM, close files. 
// The file *.kmc_suf is opened for random access
// IN	: file_name - the name of kmer_counter's output
// RET	: true		- if successful
// ----------------------------------------------------------------------------------
bool CKMCFile::OpenForRA(const std::string &file_name)
{
	uint64 size;
	size_t result;

	if(file_pre || file_suf)
		return false;

	if(!OpenASingleFile(file_name + ".kmc_pre", file_pre, size, (char *)"KMCP"))
		return false;

	ReadParamsFrom_prefix_file_buf(size);

	fclose(file_pre);
	file_pre = NULL;
		
	if(!OpenASingleFile(file_name + ".kmc_suf", file_suf, size, (char *)"KMCS"))
		return false;

	sufix_file_buf = new uchar[size];
	result = fread (sufix_file_buf, 1, size, file_suf);
	if(result == 0)
		return false;

	fclose(file_suf);
	file_suf = NULL;

	is_opened = opened_for_RA;
	prefix_index = 0;
	sufix_number = 0;
	return true;
}

//----------------------------------------------------------------------------------
// Open files *kmc_pre & *.kmc_suf, read *.kmc_pre to RAM, close *kmc.pre
// *.kmc_suf is buffered
// IN	: file_name - the name of kmer_counter's output
// RET	: true		- if successful
//----------------------------------------------------------------------------------
bool CKMCFile::OpenForListing(const std::string &file_name)
{
	uint64 size;
	size_t result;

	if(is_opened)
		return false;
	
	if(file_pre || file_suf)
		return false;

	if(!OpenASingleFile(file_name + ".kmc_pre", file_pre, size, (char *)"KMCP"))
		return false;

	ReadParamsFrom_prefix_file_buf(size);
	fclose(file_pre);
	file_pre = NULL;

	end_of_file = total_kmers == 0;

	if(!OpenASingleFile(file_name + ".kmc_suf", file_suf, size, (char *)"KMCS"))
		return false;

	sufix_file_buf = new uchar[part_size];
	result = fread (sufix_file_buf, 1, part_size, file_suf);
	if(result == 0)
		return false;

	is_opened = opened_for_listing;
	prefix_index = 0;
	sufix_number = 0;
	index_in_partial_buf = 0;
	return true;
}
//----------------------------------------------------------------------------------
CKMCFile::CKMCFile()
{
	file_pre = NULL;	
	file_suf = NULL;

	prefix_file_buf = NULL;
	sufix_file_buf = NULL;
	signature_map = NULL;

	is_opened = closed;
	end_of_file = false;
};
//----------------------------------------------------------------------------------	
CKMCFile::~CKMCFile()
{
	if(file_pre)
		fclose(file_pre);	
	if(file_suf)
		fclose(file_suf);
	if(prefix_file_buf)
		delete [] prefix_file_buf;
	if(sufix_file_buf)
		delete [] sufix_file_buf;
	if (signature_map)
		delete[] signature_map;
};
//----------------------------------------------------------------------------------	
// Open a file, recognize its size and check its marker. Auxiliary function.
// IN	: file_name - the name of a file to open
// RET	: true		- if successful
//----------------------------------------------------------------------------------
bool CKMCFile::OpenASingleFile(const std::string &file_name, FILE *&file_handler, uint64 &size, char marker[])
{
	char _marker[4];
	size_t result;

 	if((file_handler = my_fopen(file_name.c_str(), "rb")) == NULL)
			return false;
	
	my_fseek(file_handler, 0, SEEK_END);
	size = my_ftell(file_handler);					//the size of a whole file
	
	my_fseek(file_handler, -4, SEEK_CUR);
	result = fread (_marker, 1, 4, file_handler);
	if(result == 0)
		return false;

	size = size - 4;							//the size of the file without the terminal marker
	if (strncmp (marker, _marker, 4) != 0)
	{
		fclose(file_handler);	
		file_handler = NULL;
		return false;
	}

	rewind (file_handler);
	result = fread (_marker, 1, 4, file_handler);
	if(result == 0)
		return false;

	size = size - 4;							//the size of the file without initial and terminal markers 

	if (strncmp (marker, _marker, 4) != 0)
	{
		fclose(file_handler);	
		file_handler = NULL;
		return false;
	}

	return true;
};
//-------------------------------------------------------------------------------------
// Recognize current parameters from kmc_databese. Auxiliary function.
// IN	: the size of the file *.kmc_pre, without initial and terminal markers 
// RET	: true - if succesfull
//----------------------------------------------------------------------------------
bool CKMCFile::ReadParamsFrom_prefix_file_buf(uint64 &size)
{
	size_t result;

	my_fseek(file_pre, -8, SEEK_END);
	
	int64 header_offset;
	header_offset = fgetc(file_pre);
		
	size = size - 4;	//file size without the size of header_offset (and without 2 markers)

	my_fseek(file_pre, (0LL - (header_offset + 8)), SEEK_END);
	result = fread(&kmer_length, 1, sizeof(uint32), file_pre);
	result = fread(&mode, 1, sizeof(uint32), file_pre);
	result = fread(&counter_size, 1, sizeof(uint32), file_pre);
	result = fread(&lut_prefix_length, 1, sizeof(uint32), file_pre);
	result = fread(&signature_len, 1, sizeof(uint32), file_pre);
	result = fread(&min_count, 1, sizeof(uint32), file_pre);
	original_min_count = min_count;
	result = fread(&max_count, 1, sizeof(uint32), file_pre);
	original_max_count = max_count;
	result = fread(&total_kmers, 1, sizeof(uint64), file_pre);

	signature_map_size = ((1 << (2 * signature_len)) + 1);
	uint64 lut_area_size_in_bytes = size - (signature_map_size * sizeof(uint32) + header_offset + 8);
	single_LUT_size = 1 << (2 * lut_prefix_length);
	uint64 last_data_index = lut_area_size_in_bytes / sizeof(uint64);

	rewind(file_pre);
	my_fseek(file_pre, +4, SEEK_CUR);
	prefix_file_buf_size = (lut_area_size_in_bytes + 8) / sizeof(uint64);		//reads without 4 bytes of a header_offset (and without markers)		
	prefix_file_buf = new uint64[prefix_file_buf_size];
	result = fread(prefix_file_buf, 1, (size_t)(lut_area_size_in_bytes + 8), file_pre);
	if (result == 0)
		return false;
	prefix_file_buf[last_data_index] = total_kmers + 1;

	signature_map = new uint32[signature_map_size];
	result = fread(signature_map, 1, signature_map_size * sizeof(uint32), file_pre);
	if (result == 0)
		return false;

	sufix_size = (kmer_length - lut_prefix_length) / 4;		 
	
	sufix_rec_size = sufix_size + counter_size;	

	return true;
}
//------------------------------------------------------------------------------------------
// Check if kmer exists. 
// IN : kmer  - kmer
// OUT: count - kmer's counter if kmer exists
// RET: true  - if kmer exists
//------------------------------------------------------------------------------------------
bool CKMCFile::CheckKmer(CKmerAPI &kmer, float &count)
{
	if(is_opened != opened_for_RA)
		return false;
	if(end_of_file)
		return false;
	
	uint32 signature = kmer.get_signature(signature_len);
	
	uint32 bin_start_pos = signature_map[signature];
	bin_start_pos *= single_LUT_size;

	//recognize a prefix:
	uint64 pattern_prefix_value = kmer.kmer_data[0];

	uint32 pattern_offset = (sizeof(pattern_prefix_value) * 8) - (lut_prefix_length * 2) - (kmer.byte_alignment * 2);
	
	pattern_prefix_value = pattern_prefix_value >> pattern_offset;  //complements with 0
	if(pattern_prefix_value >= prefix_file_buf_size)
		return false;
	//look into the array with data

	int64 index_start = *(prefix_file_buf + bin_start_pos + pattern_prefix_value);			
	int64 index_stop = *(prefix_file_buf + bin_start_pos + pattern_prefix_value + 1) - 1;	
 
	uchar *sufix_byte_ptr; 
	uint64 sufix = 0;
	
										//sufix_offset is always 56
	uint32 sufix_offset = 56;			// the ofset of a sufix is for shifting the sufix towards MSB, to compare the sufix with a pattern
										// Bytes of a pattern to search are always shifted towards MSB
	
	uint32 row_index = 0;				// the number of a current row in an array kmer_data
	
	bool found = false;

	//binary search:

	while (index_start <= index_stop) 
	{
		int64 mid_index = (index_start + index_stop) / 2; 
		sufix_byte_ptr = &sufix_file_buf[mid_index * sufix_rec_size];

		uint64 pattern = 0;
	  
		pattern_offset = (lut_prefix_length + kmer.byte_alignment ) * 2;		
	  
		for(uint32 a = 0; a < sufix_size; a ++)		//check byte by byte
		{
			pattern = kmer.kmer_data[row_index];
			pattern = pattern << pattern_offset;
			pattern = pattern & 0xff00000000000000;
			
			sufix = sufix_byte_ptr[a];
			sufix = sufix << sufix_offset;
		
			if(pattern != sufix)					
				break;

			pattern_offset += 8;
			
			if (pattern_offset == 64)				//the end of a word
			{
					pattern_offset = 0;
					row_index++;
			}
		}

		if(pattern == sufix)
		{
		  found = true;
		  break;
		}
		if( sufix < pattern )
			index_start = mid_index + 1;
		else
			index_stop = mid_index - 1;
	}
	
	if(found)
	{
		sufix_byte_ptr += sufix_size;
		uint32 int_counter;
	
		int_counter = *sufix_byte_ptr;

		for(uint32 b = 1; b < counter_size; b ++)
		{
			uint32 aux = 0x000000ff & *(sufix_byte_ptr + b);

			aux = aux << 8 * ( b);
			int_counter = aux | int_counter;
		}
	
		if(mode == 0)
			count = (float)int_counter;
		else
			memcpy(&count, &int_counter, counter_size);
		
		if((count >= min_count) && (count <= max_count))
			return true;
		else
			return false;
	}
	return false;
}

//-----------------------------------------------------------------------------------------------
// Check if end of file
// RET: true - all kmers are listed
//-----------------------------------------------------------------------------------------------
bool CKMCFile::Eof(void)
{
	return end_of_file;	
}
//-----------------------------------------------------------------------------------------------
// Read next kmer
// OUT: kmer - next kmer
// OUT: count - kmer's counter
// RET: true - if not EOF
//-----------------------------------------------------------------------------------------------
bool CKMCFile::ReadNextKmer(CKmerAPI &kmer, float &count)
{
	uint32 int_counter;

	if(is_opened != opened_for_listing)
		return false;
	do
	{
		if(end_of_file)
			return false;
		
		if(sufix_number == prefix_file_buf[prefix_index + 1]) 
		{
			prefix_index++;
						
			while (prefix_file_buf[prefix_index] == prefix_file_buf[prefix_index + 1])
				prefix_index++;
		}
	
		uint32 off = (sizeof(prefix_index) * 8) - (lut_prefix_length * 2) - kmer.byte_alignment * 2;
			
		uint64 temp_prefix = prefix_index << off;	// shift prefix towards MSD
		
		kmer.kmer_data[0] = temp_prefix;			// store prefix in an object CKmerAPI

		for(uint32 i = 1; i < kmer.no_of_rows; i++)
			kmer.kmer_data[i] = 0;

		//read sufix:
		uint32 row_index = 0;
 		uint64 suf = 0;
	
		off = off - 8;
				
 		for(uint32 a = 0; a < sufix_size; a ++)
		{
			if(index_in_partial_buf == part_size)
				Reload_sufix_file_buf();
						
			suf = sufix_file_buf[index_in_partial_buf++];
			suf = suf << off;
			kmer.kmer_data[row_index] = kmer.kmer_data[row_index] | suf;

			if (off == 0)				//the end of a word in kmer_data
			{
					off = 56;
					row_index++;
			}
			else
					off -=8;
		}
	
		//read counter:
		if(index_in_partial_buf == part_size)
			Reload_sufix_file_buf();
		
		int_counter = sufix_file_buf[index_in_partial_buf++];

		for(uint32 b = 1; b < counter_size; b++)
		{
			if(index_in_partial_buf == part_size)
				Reload_sufix_file_buf();
			
			uint32 aux = 0x000000ff & sufix_file_buf[index_in_partial_buf++];
			aux = aux << 8 * ( b);
			int_counter = aux | int_counter;
		}
	
		if(mode == 0)
			count = (float)int_counter;
		else
			memcpy(&count, &int_counter, counter_size);
	
		sufix_number++;
	
		if(sufix_number == total_kmers)
			end_of_file = true;
	}
	while((count < min_count) || (count > max_count));

	return true;
}
//-------------------------------------------------------------------------------
// Reload a contents of an array "sufix_file_buf" for listing mode. Auxiliary function.
//-------------------------------------------------------------------------------
void CKMCFile::Reload_sufix_file_buf()
{
		fread (sufix_file_buf, 1, (size_t) part_size, file_suf);
		index_in_partial_buf = 0;
};
//-------------------------------------------------------------------------------
// Release memory and close files in case they were opened 
// RET: true - if files have been readed
//-------------------------------------------------------------------------------
bool CKMCFile::Close()
{
	if(is_opened)
	{
		if(file_pre)
		{
			fclose(file_pre);	
			file_pre = NULL;
		}
		if(file_suf)
		{
			fclose(file_suf);
			file_suf = NULL;
		}
	
		is_opened = closed;
		end_of_file = false;
		delete [] prefix_file_buf;
		prefix_file_buf = NULL;
		delete [] sufix_file_buf;
		sufix_file_buf = NULL;
		delete[] signature_map;
		signature_map = NULL;

		return true;
	}
	else
		return false;
};
//----------------------------------------------------------------------------------
// Set initial values to enable listing kmers from the begining. Only in listing mode
// RET: true - if a file has been opened for listing
//----------------------------------------------------------------------------------
bool CKMCFile::RestartListing(void)
{
	if(is_opened == opened_for_listing)
	{
		
		my_fseek ( file_suf , 4 , SEEK_SET );
		fread (sufix_file_buf, 1, (size_t) part_size, file_suf);

		prefix_index = 0;
		sufix_number = 0;
		index_in_partial_buf = 0;

		end_of_file = total_kmers == 0;

		return true;
	}
	return false;
		
};
//----------------------------------------------------------------------------------------
// Set the minimal value for a counter. Kmers with counters below this theshold are ignored
// IN	: x - minimal value for a counter
// RET	: true - if successful 
//----------------------------------------------------------------------------------------
bool CKMCFile::SetMinCount(uint32 x)
{
	if((original_min_count <= x) && (x < max_count))
	{
		min_count = x;
		return true;
	} 
	else
		return false;
}

//----------------------------------------------------------------------------------------
// Return a value of min_count. Kmers with counters below this theshold are ignored 
// RET	: a value of min_count
//----------------------------------------------------------------------------------------
uint32 CKMCFile::GetMinCount(void)
{
	return min_count;
};

//----------------------------------------------------------------------------------------
// Set the maximal value for a counter. Kmers with counters above this theshold are ignored
// IN	: x - maximal value for a counter
// RET	: true - if successful 
//----------------------------------------------------------------------------------------
bool CKMCFile::SetMaxCount(uint32 x)
{
	if((original_max_count >= x) && (x > min_count))
	{
		max_count = x;
		return true; 
	}
	else
		return false;
}


//----------------------------------------------------------------------------------------
// Return a value of max_count. Kmers with counters above this theshold are ignored 
// RET	: a value of max_count
//----------------------------------------------------------------------------------------
uint32 CKMCFile::GetMaxCount(void)
{
	return max_count;
}

//----------------------------------------------------------------------------------------
// Set original (readed from *.kmer_pre) values for min_count and max_count
//----------------------------------------------------------------------------------------
void CKMCFile::ResetMinMaxCounts(void)
{
	min_count = original_min_count;
	max_count = original_max_count;
} 

//----------------------------------------------------------------------------------------
// Return the length of kmers
// RET	: the length of kmers
//----------------------------------------------------------------------------------------
uint32 CKMCFile::KmerLength(void)
{
	return kmer_length;			
}

//----------------------------------------------------------------------------------------
// Check if kmer exists
// IN	: kmer - kmer
// RET	: true if kmer exists
//----------------------------------------------------------------------------------------
bool CKMCFile::IsKmer(CKmerAPI &kmer)
{
	float _count;
	if(CheckKmer(kmer, _count))
		return true;
	else
		return false;
}

//-----------------------------------------------------------------------------------------
// Check the total number of kmers between current min_count and max_count
// RET	: total number of kmers or 0 if a database has not been opened
//-----------------------------------------------------------------------------------------
uint64 CKMCFile::KmerCount(void)
{
	if(is_opened)
		if((min_count == original_min_count) && (max_count == original_max_count))
			return total_kmers;
		else
		{
			uint32 count;
			uint32 int_counter;
			uint64 aux_kmerCount = 0;

			if(is_opened == opened_for_RA)
			{
				uchar *ptr = sufix_file_buf;
				
				for(uint64 i = 0; i < total_kmers; i++)		
				{
					ptr += sufix_size;
					int_counter = *ptr;
					ptr++;

					for(uint32 b = 1; b < counter_size; b ++)
					{
						uint32 aux = 0x000000ff & *(ptr);
						aux = aux << 8 * ( b);
						int_counter = aux | int_counter;
						ptr++;
					}
					
					if(mode == 0)
						count = int_counter;
					else
						memcpy(&count, &int_counter, counter_size);
	
					if((count >= min_count) && (count <= max_count))
						aux_kmerCount++;
				}
			}
			else //opened_for_listing
			{
				CKmerAPI kmer(kmer_length);
				float count;
				RestartListing();
				for(uint64 i = 0; i < total_kmers; i++)		
				{
					ReadNextKmer(kmer, count);
					if((count >= min_count) && (count <= max_count))
						aux_kmerCount++;
				}
				RestartListing();
			}
			return aux_kmerCount;
		}
	else
		return 0 ;
}
//---------------------------------------------------------------------------------
// Get current parameters from kmer_database
// OUT	:	_kmer_length	- the length of kmers
//			_mode			- mode
//			_counter_size	- the size of a counter in bytes 
//			_lut_prefix_length - the number of prefix's symbols cut from kmers 
//			_min_count		- the minimal number of kmer's appearances 
//			_max_count		- the maximal number of kmer's appearances
//			_total_kmers	- the total number of kmers
// RET	: true if kmer_database has been opened
//---------------------------------------------------------------------------------
bool CKMCFile::Info(uint32 &_kmer_length, uint32 &_mode, uint32 &_counter_size, uint32 &_lut_prefix_length, uint32 &_signature_len, uint32 &_min_count, uint32 &_max_count, uint64 &_total_kmers)
{
	if(is_opened)
	{
		_kmer_length = kmer_length;
		_mode = mode;
		_counter_size = counter_size;
		_lut_prefix_length = lut_prefix_length;
		_signature_len = signature_len;
		_min_count = min_count;
		_max_count = max_count;
		_total_kmers = total_kmers;
		return true;
	}
	return false;
};

// ***** EOF
