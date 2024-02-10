/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.4
  Date   : 2024-02-09
*/

#include "mmer.h"
#include "kmc_file.h"
#include <tuple>
#include <string>

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

	if (file_pre || file_suf)
		return false;

	if (!OpenASingleFile(file_name + ".kmc_pre", file_pre, size, (char *)"KMCP"))
		return false;

	ReadParamsFrom_prefix_file_buf(size, open_mode::opened_for_RA);

	if (!OpenASingleFile(file_name + ".kmc_suf", file_suf, size, (char *)"KMCS"))
		return false;

	sufix_file_buf = new uchar[size];
	result = fread(sufix_file_buf, 1, size, file_suf);
	if (result != size)
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

	if (is_opened)
		return false;
	
	if (file_pre || file_suf)
		return false;

	if (!OpenASingleFile(file_name + ".kmc_pre", file_pre, size, (char *)"KMCP"))
		return false;

	ReadParamsFrom_prefix_file_buf(size, open_mode::opened_for_listing);

	end_of_file = total_kmers == 0;

	if (!OpenASingleFile(file_name + ".kmc_suf", file_suf, size, (char *)"KMCS"))
		return false;

	sufix_file_buf = new uchar[part_size];

	suffix_file_total_to_read = size;
	suf_file_left_to_read = suffix_file_total_to_read;
	auto to_read = MIN(suf_file_left_to_read, part_size);
	auto readed = fread(sufix_file_buf, 1, to_read, file_suf);
	if (readed != to_read)
	{
		std::cerr << "Error: some error while reading suffix file\n";
		return false;
	}

	suf_file_left_to_read -= readed;

	is_opened = opened_for_listing;
	prefix_index = 0;
	sufix_number = 0;
	index_in_partial_buf = 0;
	return true;
}
//----------------------------------------------------------------------------------
// Open files *kmc_pre & *.kmc_suf
// Following k-mers will be read in bins which order is defined by the second file
// IN	: file_name - the name of kmer_counter's output
// RET	: true		- if successful
//----------------------------------------------------------------------------------
bool CKMCFile::OpenForListingWithBinOrder(const std::string& file_name, const std::string& bin_order_file_name)
{
	uint64 size_pre_file;
	uint64 size_suf_file;

	if (is_opened)
		return false;

	if (file_pre || file_suf)
		return false;

	if (!OpenASingleFile(file_name + ".kmc_pre", file_pre, size_pre_file, (char*)"KMCP"))
		return false;

	end_of_file = total_kmers == 0;

	if (!OpenASingleFile(file_name + ".kmc_suf", file_suf, size_suf_file, (char*)"KMCS"))
		return false;

	sufix_file_buf = new uchar[part_size];

	ReadParamsFrom_prefix_file_buf(size_pre_file, open_mode::opened_for_listing_with_bin_order, bin_order_file_name);

	is_opened = opened_for_listing_with_bin_order;

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
}
//----------------------------------------------------------------------------------	
CKMCFile::~CKMCFile()
{
	if (file_pre)
		fclose(file_pre);
	if (file_suf)
		fclose(file_suf);
	if (prefix_file_buf)
		delete[] prefix_file_buf;
	if (sufix_file_buf)
		delete[] sufix_file_buf;
	if (signature_map)
		delete[] signature_map;
}
//----------------------------------------------------------------------------------	
// Open a file, recognize its size and check its marker. Auxiliary function.
// IN	: file_name - the name of a file to open
// RET	: true		- if successful
//----------------------------------------------------------------------------------
bool CKMCFile::OpenASingleFile(const std::string &file_name, FILE *&file_handler, uint64 &size, char marker[])
{
	char _marker[4];
	size_t result;

	if ((file_handler = my_fopen(file_name.c_str(), "rb")) == NULL)
		return false;

	my_fseek(file_handler, 0, SEEK_END);
	size = my_ftell(file_handler);					//the size of a whole file

	my_fseek(file_handler, -4, SEEK_CUR);
	result = fread(_marker, 1, 4, file_handler);
	if (result == 0)
		return false;

	size = size - 4;							//the size of the file without the terminal marker
	if (strncmp(marker, _marker, 4) != 0)
	{
		fclose(file_handler);
		file_handler = NULL;
		return false;
	}

	rewind(file_handler);
	result = fread(_marker, 1, 4, file_handler);
	if (result == 0)
		return false;

	size = size - 4;							//the size of the file without initial and terminal markers 

	if (strncmp(marker, _marker, 4) != 0)
	{
		fclose(file_handler);
		file_handler = NULL;
		return false;
	}

	return true;
}
//-------------------------------------------------------------------------------------
// Recognize current parameters from kmc_databese. Auxiliary function.
// IN	: the size of the file *.kmc_pre, without initial and terminal markers 
// RET	: true - if succesfull
//----------------------------------------------------------------------------------
bool CKMCFile::ReadParamsFrom_prefix_file_buf(uint64 &size, open_mode _open_mode, const std::string& bin_order_file_name)
{
	size_t prev_pos = my_ftell(file_pre);
	my_fseek(file_pre, -12, SEEK_END);
	size_t result;

	result = fread(&kmc_version, sizeof(uint32), 1, file_pre);
	if (kmc_version != 0 && kmc_version != 0x200) //only this versions are supported, 0 = kmc1, 0x200 = kmc2
		return false;
	my_fseek(file_pre, prev_pos, SEEK_SET);

	if (kmc_version == 0x200)
	{
		my_fseek(file_pre, -8, SEEK_END);
		
		int64 header_offset;
		header_offset = fgetc(file_pre);
		
		size = size - 4;	//file size without the size of header_offset (and without 2 markers)

		my_fseek(file_pre, (0LL - (header_offset + 8)), SEEK_END);
		result = fread(&kmer_length, 1, sizeof(uint32), file_pre);
		result = fread(&mode, 1, sizeof(uint32), file_pre);
		if (mode != 0)
		{
			std::cerr << "Error: Quake quake compatible counters are not supported anymore\n";
			return false;
		}
		result = fread(&counter_size, 1, sizeof(uint32), file_pre);
		result = fread(&lut_prefix_length, 1, sizeof(uint32), file_pre);
		result = fread(&signature_len, 1, sizeof(uint32), file_pre);
		result = fread(&min_count, 1, sizeof(uint32), file_pre);
		original_min_count = min_count;
		uint32 max_count_uint32;
		result = fread(&max_count_uint32, 1, sizeof(uint32), file_pre);
		max_count = max_count_uint32;
		original_max_count = max_count;
		result = fread(&total_kmers, 1, sizeof(uint64), file_pre);
		result = fread(&both_strands, 1, 1, file_pre);
		both_strands = !both_strands;

		signature_map_size = ((1 << (2 * signature_len)) + 1);
		uint64 lut_area_size_in_bytes = size - (signature_map_size * sizeof(uint32)+header_offset + 8);
		single_LUT_size = 1 << (2 * lut_prefix_length);
		uint64 last_data_index = lut_area_size_in_bytes / sizeof(uint64);

		n_bins = last_data_index / single_LUT_size;

		signature_map = new uint32[signature_map_size];

		fseek(file_pre, 4 + lut_area_size_in_bytes + 8, SEEK_SET);
		result = fread(signature_map, 1, signature_map_size * sizeof(uint32), file_pre);
		if (result == 0)
			return false;

		sufix_size = (kmer_length - lut_prefix_length) / 4;

		sufix_rec_size = sufix_size + counter_size;

		if(_open_mode == opened_for_RA)
		{
			rewind(file_pre);
			my_fseek(file_pre, +4, SEEK_CUR);
			prefix_file_buf_size = (lut_area_size_in_bytes + 8) / sizeof(uint64);		//reads without 4 bytes of a header_offset (and without markers)
			prefix_file_buf = new uint64[prefix_file_buf_size];
			result = fread(prefix_file_buf, 1, (size_t)(lut_area_size_in_bytes + 8), file_pre);
			if (result == 0)
				return false;

			prefix_file_buf[last_data_index] = total_kmers + 1; //I think + 1 if wrong, but due to the implementation of binary search it does not matter, it was here in kmc 0.3 and I leave it this way just in case...

			result = fread(signature_map, 1, signature_map_size * sizeof(uint32), file_pre);
			if (result == 0)
				return false;

			fclose(file_pre);
			file_pre = nullptr;
		}
		else if (_open_mode == opened_for_listing)
		{
			prefixFileBufferForListingMode = std::make_unique<CPrefixFileBufferForListingMode>(file_pre, last_data_index, lut_prefix_length, false, total_kmers);
		}
		else if (_open_mode == opened_for_listing_with_bin_order)
		{
			ordered_bin_reading = std::make_unique<OrderedBinReading>(
				file_pre, 4,
				file_suf, 4,
				bin_order_file_name,
				n_bins,
				signature_len,
				signature_map,
				signature_map_size,
				single_LUT_size,
				sufix_rec_size);
		}
		else
			throw std::runtime_error("unknown _open_mode");

		return true;
	}
	else if (kmc_version == 0)
	{
		if (_open_mode == opened_for_listing_with_bin_order)
			throw std::runtime_error("opened_for_listing_with_bin_order open mode may be used only for KMC database in KMC2 format");

		my_fseek(file_pre, -8, SEEK_END);

		uint64 header_offset;
		header_offset = fgetc(file_pre);

		size = size - 4;

		my_fseek(file_pre, (0LL - (header_offset + 8)), SEEK_END);
		result = fread(&kmer_length, 1, sizeof(uint32), file_pre);
		result = fread(&mode, 1, sizeof(uint32), file_pre);
		if (mode != 0)
		{
			std::cerr << "Error: Quake quake compatible counters are not supported anymore\n";
			return false;
		}
		result = fread(&counter_size, 1, sizeof(uint32), file_pre);
		result = fread(&lut_prefix_length, 1, sizeof(uint32), file_pre);
		result = fread(&min_count, 1, sizeof(uint32), file_pre);
		original_min_count = min_count;

		uint32 max_count_lo;
		result = fread(&max_count_lo, 1, sizeof(uint32), file_pre);
		max_count = max_count_lo;
		original_max_count = max_count;
		result = fread(&total_kmers, 1, sizeof(uint64), file_pre);
		result = fread(&both_strands, 1, 1, file_pre);
		both_strands = !both_strands;

		uint32 max_count_hi;
		result = fread(&max_count_hi, 1, sizeof(uint32), file_pre);
		max_count += (uint64)max_count_hi << 32;
		original_max_count = max_count;

		prefix_file_buf_size = (1ull << (2 * lut_prefix_length)) + 1;
		uint64 last_data_index = prefix_file_buf_size - 1;

		if (_open_mode == opened_for_RA)
		{
			prefix_file_buf = new uint64[prefix_file_buf_size];
			fseek(file_pre, 4, SEEK_SET);
			result = fread(prefix_file_buf, 1, (size_t)(prefix_file_buf_size * sizeof(uint64)), file_pre);
			if (result == 0)
				return false;

			prefix_file_buf[last_data_index] = total_kmers + 1; //I think + 1 if wrong, but due to the implementation of binary search it does not matter, it was here in kmc 0.3 and I leave it this way just in case...

			fclose(file_pre);
			file_pre = nullptr;
		}
		else
		{
			prefixFileBufferForListingMode = std::make_unique<CPrefixFileBufferForListingMode>(file_pre, last_data_index, lut_prefix_length, true, total_kmers);
		}

		sufix_size = (kmer_length - lut_prefix_length) / 4;

		sufix_rec_size = sufix_size + counter_size;

		return true;

	}
	return false;
}

//-----------------------------------------------------------------------------------------------
uint32_t CKMCFile::GetNBins() const
{
	if (is_opened != opened_for_listing_with_bin_order)
		return 0;
	return ordered_bin_reading->get_n_bins();
}

//-----------------------------------------------------------------------------------------------
bool CKMCFile::StartBin()
{
	if (is_opened != opened_for_listing_with_bin_order)
		return false;
	return ordered_bin_reading->next_bin();
}

//-----------------------------------------------------------------------------------------------
bool CKMCFile::StartBin(uint32_t bin_id)
{
	if (is_opened != opened_for_listing_with_bin_order)
		return false;
	ordered_bin_reading->start_bin(bin_id);
	return true;
}

//-----------------------------------------------------------------------------------------------
bool CKMCFile::ReadNextKmerFromBin(CKmerAPI& kmer, uint64& count)
{
	if (is_opened != opened_for_listing_with_bin_order)
		return false;

	uint32 off = (sizeof(prefix_index) * 8) - (lut_prefix_length * 2) - kmer.byte_alignment * 2;
	return ordered_bin_reading->read_from_cur_bin(kmer, count, off, sufix_size, counter_size, min_count, max_count);
}

//------------------------------------------------------------------------------------------
// Check if kmer exists. 
// IN : kmer  - kmer
// OUT: count - kmer's counter if kmer exists
// RET: true  - if kmer exists
//------------------------------------------------------------------------------------------
bool CKMCFile::CheckKmer(CKmerAPI &kmer, uint32 &count)
{
	if(is_opened != opened_for_RA)
		return false;
	if(end_of_file)
		return false;
	
	//recognize a prefix:
	uint64 pattern_prefix_value = kmer.kmer_data[0];

	uint32 pattern_offset = (sizeof(pattern_prefix_value)* 8) - (lut_prefix_length * 2) - (kmer.byte_alignment * 2);
	int64 index_start = 0, index_stop = 0;

	pattern_prefix_value = pattern_prefix_value >> pattern_offset;  //complements with 0
	if (pattern_prefix_value >= prefix_file_buf_size)
		return false;

	if (kmc_version == 0x200)
	{
		uint32 signature = kmer.get_signature(signature_len);
		uint32 bin_start_pos = signature_map[signature];
		bin_start_pos *= single_LUT_size;				
		//look into the array with data
		index_start = *(prefix_file_buf + bin_start_pos + pattern_prefix_value);
		index_stop = *(prefix_file_buf + bin_start_pos + pattern_prefix_value + 1) - 1;
	}
	else if (kmc_version == 0)
	{
		//look into the array with data
		index_start = prefix_file_buf[pattern_prefix_value];
		index_stop = prefix_file_buf[pattern_prefix_value + 1] - 1;
	}
	uint64 tmp_count ;
	bool res = BinarySearch(index_start, index_stop, kmer, tmp_count, pattern_offset);
	count = (uint32)tmp_count;
	return res;
}

//------------------------------------------------------------------------------------------
// Check if kmer exists. 
// IN : kmer  - kmer
// OUT: count - kmer's counter if kmer exists
// RET: true  - if kmer exists
//------------------------------------------------------------------------------------------
bool CKMCFile::CheckKmer(CKmerAPI &kmer, uint64 &count)
{
	if (is_opened != opened_for_RA)
		return false;
	if (end_of_file)
		return false;

	//recognize a prefix:
	uint64 pattern_prefix_value = kmer.kmer_data[0];

	uint32 pattern_offset = (sizeof(pattern_prefix_value)* 8) - (lut_prefix_length * 2) - (kmer.byte_alignment * 2);
	int64 index_start = 0, index_stop = 0;

	pattern_prefix_value = pattern_prefix_value >> pattern_offset;  //complements with 0
	if (pattern_prefix_value >= prefix_file_buf_size)
		return false;

	if (kmc_version == 0x200)
	{
		uint32 signature = kmer.get_signature(signature_len);
		uint32 bin_start_pos = signature_map[signature];
		bin_start_pos *= single_LUT_size;
		//look into the array with data
		index_start = *(prefix_file_buf + bin_start_pos + pattern_prefix_value);
		index_stop = *(prefix_file_buf + bin_start_pos + pattern_prefix_value + 1) - 1;
	}
	else if (kmc_version == 0)
	{
		//look into the array with data
		index_start = prefix_file_buf[pattern_prefix_value];
		index_stop = prefix_file_buf[pattern_prefix_value + 1] - 1;
	}
	return BinarySearch(index_start, index_stop, kmer, count, pattern_offset);
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
bool CKMCFile::ReadNextKmer(CKmerAPI &kmer, uint32 &count)
{
	if(is_opened != opened_for_listing)
		return false;
	do
	{
		if(end_of_file)
			return false;
	
		uint32 off = (sizeof(prefix_index) * 8) - (lut_prefix_length * 2) - kmer.byte_alignment * 2;
		uint64 prefix = prefixFileBufferForListingMode->GetPrefix(sufix_number);

		uint64 temp_prefix = prefix << off;	// shift prefix towards MSD. "& prefix_mask" necessary for kmc2 db format
		
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
		if (counter_size == 0)
			count = 1;
		else
		{
			
			if (index_in_partial_buf == part_size)
				Reload_sufix_file_buf();

			count = sufix_file_buf[index_in_partial_buf++];

			for (uint32 b = 1; b < counter_size; b++)
			{
				if (index_in_partial_buf == part_size)
					Reload_sufix_file_buf();

				uint32 aux = 0x000000ff & sufix_file_buf[index_in_partial_buf++];
				aux = aux << 8 * (b);
				count = aux | count;
			}
		}
		sufix_number++;
	
		if(sufix_number == total_kmers)
			end_of_file = true;
	}
	while ((counter_size != 0) && ((count < min_count) || (count > max_count))); //do not applay filtering if counter_size == 0 as it does not make sense

	return true;
}

//-----------------------------------------------------------------------------------------------
// Read next kmer
// OUT: kmer - next kmer
// OUT: count - kmer's counter
// RET: true - if not EOF
//-----------------------------------------------------------------------------------------------
bool CKMCFile::ReadNextKmer(CKmerAPI &kmer, uint64 &count)
{
	if (is_opened != opened_for_listing)
		return false;
	do
	{
		if (end_of_file)
			return false;

		uint32 off = (sizeof(prefix_index)* 8) - (lut_prefix_length * 2) - kmer.byte_alignment * 2;

		uint64_t prefix = prefixFileBufferForListingMode->GetPrefix(sufix_number);
		uint64 temp_prefix = prefix << off;	// shift prefix towards MSD. "& prefix_mask" necessary for kmc2 db format

		kmer.kmer_data[0] = temp_prefix;			// store prefix in an object CKmerAPI

		for (uint32 i = 1; i < kmer.no_of_rows; i++)
			kmer.kmer_data[i] = 0;

		//read sufix:
		uint32 row_index = 0;
		uint64 suf = 0;

		off = off - 8;

		for (uint32 a = 0; a < sufix_size; a++)
		{
			if (index_in_partial_buf == part_size)
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
				off -= 8;
		}

		//read counter:
		if (counter_size == 0)
			count = 1;
		else
		{
			if (index_in_partial_buf == part_size)
				Reload_sufix_file_buf();

			count = sufix_file_buf[index_in_partial_buf++];

			for (uint32 b = 1; b < counter_size; b++)
			{
				if (index_in_partial_buf == part_size)
					Reload_sufix_file_buf();

				uint64 aux = 0x000000ff & sufix_file_buf[index_in_partial_buf++];
				aux = aux << 8 * (b);
				count = aux | count;
			}
		}
		sufix_number++;

		if (sufix_number == total_kmers)
			end_of_file = true;

	} while ((counter_size != 0) && ((count < min_count) || (count > max_count))); //do not applay filtering if counter_size == 0 as it does not make sense

	return true;
}

//-------------------------------------------------------------------------------
// Reload a contents of an array "sufix_file_buf" for listing mode. Auxiliary function.
//-------------------------------------------------------------------------------
void CKMCFile::Reload_sufix_file_buf()
{
	auto to_read = MIN(suf_file_left_to_read, part_size);
	auto readed = fread(sufix_file_buf, 1, (size_t)to_read, file_suf);
	suf_file_left_to_read -= readed;
	if (readed != to_read)
	{
		//TODO: what about this exit? maybe turn into exception
		std::cerr << "Error: some error while reading suffix file\n";
		exit(1);
	}
	index_in_partial_buf = 0;
}
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
}
//----------------------------------------------------------------------------------
// Set initial values to enable listing kmers from the begining. Only in listing mode
// RET: true - if a file has been opened for listing
//----------------------------------------------------------------------------------
bool CKMCFile::RestartListing(void)
{
	if(is_opened == opened_for_listing)
	{
		my_fseek(file_suf , 4 , SEEK_SET);
		suf_file_left_to_read = suffix_file_total_to_read;
		auto to_read = MIN(suf_file_left_to_read, part_size);
		auto readed = fread(sufix_file_buf, 1, to_read, file_suf);
		if (readed != to_read)
		{
			std::cerr << "Error: some error while reading suffix file\n";
			return false;
		}

		suf_file_left_to_read -= readed;
		prefix_index = 0;
		sufix_number = 0;
		index_in_partial_buf = 0;

		end_of_file = total_kmers == 0;

		return true;
	}
	return false;
		
}
//----------------------------------------------------------------------------------------
// Set the minimal value for a counter. Kmers with counters below this theshold are ignored
// IN	: x - minimal value for a counter
// RET	: true - if successful 
//----------------------------------------------------------------------------------------
bool CKMCFile::SetMinCount(uint32 x)
{
	if((original_min_count <= x) && (x <= max_count))
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
}

//----------------------------------------------------------------------------------------
// Set the maximal value for a counter. Kmers with counters above this theshold are ignored
// IN	: x - maximal value for a counter
// RET	: true - if successful 
//----------------------------------------------------------------------------------------
bool CKMCFile::SetMaxCount(uint32 x)
{
	if((original_max_count >= x) && (x >= min_count))
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
uint64 CKMCFile::GetMaxCount(void)
{
	return max_count;
}

//----------------------------------------------------------------------------------------
// Return true if KMC was run without -b switch
// RET	: a value of both_strands
//----------------------------------------------------------------------------------------
bool CKMCFile::GetBothStrands(void)
{
	return both_strands;
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
	uint32 _count;
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
			uint64 aux_kmerCount = 0;

			if(is_opened == opened_for_RA)
			{
				uchar *ptr = sufix_file_buf;
				
				for(uint64 i = 0; i < total_kmers; i++)		
				{
					ptr += sufix_size;

					if (counter_size == 0)
						count = 1;
					else
					{
						count = *ptr;
						ptr++;

						for(uint32 b = 1; b < counter_size; b ++)
						{
							uint32 aux = 0x000000ff & *(ptr);
							aux = aux << 8 * ( b);
							count = aux | count;
							ptr++;
						}
					}
					//applay filtering only if counter_size != 0
					if((counter_size == 0) || ((count >= min_count) && (count <= max_count)))
						aux_kmerCount++;
				}
			}
			else //opened_for_listing
			{
				CKmerAPI kmer(kmer_length);
				uint32 count;
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
bool CKMCFile::Info(uint32 &_kmer_length, uint32 &_mode, uint32 &_counter_size, uint32 &_lut_prefix_length, uint32 &_signature_len, uint32 &_min_count, uint64 &_max_count, uint64 &_total_kmers)
{
	if(is_opened)
	{
		_kmer_length = kmer_length;
		_mode = mode;
		_counter_size = counter_size;
		_lut_prefix_length = lut_prefix_length;
		if (kmc_version == 0x200)
			_signature_len = signature_len;
		else
			_signature_len = 0; //for kmc1 there is no signature_len
		_min_count = min_count;
		_max_count = max_count;
		_total_kmers = total_kmers;
		return true;
	}
	return false;
}

// Get current parameters from kmer_database
bool CKMCFile::Info(CKMCFileInfo& info)
{
	if (is_opened)
	{
		info.kmer_length = kmer_length;
		info.mode = mode;
		info.counter_size = counter_size;
		info.lut_prefix_length = lut_prefix_length;
		if (kmc_version == 0x200)
			info.signature_len = signature_len;
		else
			info.signature_len = 0; //for kmc1 there is no signature_len
		info.min_count = min_count;
		info.max_count = max_count;
		info.total_kmers = total_kmers;
		info.both_strands = both_strands;
		return true;
	}
	return false;
}


//---------------------------------------------------------------------------------
// Get counters from read
// OUT	:	counters    	- vector of counters of each k-mer in read (of size read_len - kmer_len + 1), if some k-mer is invalid (i.e. contains 'N') the counter is equal to 0
// IN   :   read			- 
// RET	:   true if success, false if k > read length or some failure 
//---------------------------------------------------------------------------------
bool CKMCFile::GetCountersForRead(const std::string& read, std::vector<uint32>& counters)
{
	if (is_opened != opened_for_RA)
		return false;

	if (read.length() < kmer_length)
	{
		counters.clear();
		return false;
	}

	if (kmc_version == 0x200)
	{		
		if (both_strands)
			return GetCountersForRead_kmc2_both_strands(read, counters);
		else
			return GetCountersForRead_kmc2(read, counters);
	}
	else if (kmc_version == 0)
	{
		if (both_strands)
			return GetCountersForRead_kmc1_both_strands(read,counters);
		else
			return GetCountersForRead_kmc1(read, counters);
	}
	else
		return false; //never should be here
}

//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
uint32 CKMCFile::count_for_kmer_kmc1(CKmerAPI& kmer)
{
	//recognize a prefix:

	uint64 pattern_prefix_value = kmer.kmer_data[0];

	uint32 pattern_offset = (sizeof(pattern_prefix_value)* 8) - (lut_prefix_length * 2) - (kmer.byte_alignment * 2);

	pattern_prefix_value = pattern_prefix_value >> pattern_offset;  //complements with 0
	if (pattern_prefix_value >= prefix_file_buf_size)
		return false;
	//look into the array with data

	int64 index_start = prefix_file_buf[pattern_prefix_value];
	int64 index_stop = prefix_file_buf[pattern_prefix_value + 1] - 1;

	uint64 counter = 0;
	if (BinarySearch(index_start, index_stop, kmer, counter, pattern_offset))
		return (uint32)counter;
	return 0;
}

//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
uint32 CKMCFile::count_for_kmer_kmc2(CKmerAPI& kmer, uint32 bin_start_pos)
{
	//recognize a prefix:
	uint64 pattern_prefix_value = kmer.kmer_data[0];

	uint32 pattern_offset = (sizeof(pattern_prefix_value)* 8) - (lut_prefix_length * 2) - (kmer.byte_alignment * 2);

	pattern_prefix_value = pattern_prefix_value >> pattern_offset;  //complements with 0
	if (pattern_prefix_value >= prefix_file_buf_size)
		return false;
	//look into the array with data

	int64 index_start = *(prefix_file_buf + bin_start_pos + pattern_prefix_value);
	int64 index_stop = *(prefix_file_buf + bin_start_pos + pattern_prefix_value + 1) - 1;

	uint64 counter = 0;
	if (BinarySearch(index_start, index_stop, kmer, counter, pattern_offset))
		return (uint32)counter;
	return 0;
}

//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
bool CKMCFile::GetCountersForRead_kmc1_both_strands(const std::string& read, std::vector<uint32>& counters)
{
	uint32 read_len = static_cast<uint32>(read.length());
	counters.resize(read.length() - kmer_length + 1);
	std::string transformed_read = read;
	for (char& c : transformed_read)
		c = CKmerAPI::num_codes[(uchar)c];

	uint32 i = 0;
	CKmerAPI kmer(kmer_length), kmer_rev(kmer_length);
	uint32 pos = 0;
	uint32 rev_pos = kmer_length - 1;

	uint32 counters_pos = 0;

	while (i + kmer_length - 1 < read_len)
	{
		bool contains_N = false;
		while (i < read_len && pos < kmer_length)
		{
			if (CKmerAPI::num_codes[(uchar)read[i]] < 0)
			{
				pos = 0;
				rev_pos = kmer_length - 1;
				kmer.clear();
				kmer_rev.clear();
				++i;
				uint32 wrong_kmers = MIN(i - counters_pos, static_cast<uint32>(counters.size()) - counters_pos);
				fill_n(counters.begin() + counters_pos, wrong_kmers, 0);
				counters_pos += wrong_kmers;
				contains_N = true;
				break;
			}
			else
			{
				kmer_rev.insert2bits(rev_pos--, 3 - CKmerAPI::num_codes[(uchar)read[i]]);
				kmer.insert2bits(pos++, CKmerAPI::num_codes[(uchar)read[i++]]);
				
			}
		}
		if (contains_N)
			continue;
		if (pos == kmer_length)
		{
			if(kmer < kmer_rev)
				counters[counters_pos++] = count_for_kmer_kmc1(kmer);
			else
				counters[counters_pos++] = count_for_kmer_kmc1(kmer_rev);
		}
		else
			break;

		while (i < read_len)
		{
			if (CKmerAPI::num_codes[(uchar)read[i]] < 0)
			{
				pos = 0;
				break;
			}
			kmer_rev.SHR_insert2bits(3 - CKmerAPI::num_codes[(uchar)read[i]]);
			kmer.SHL_insert2bits(CKmerAPI::num_codes[(uchar)read[i++]]);
			if(kmer < kmer_rev)
				counters[counters_pos++] = count_for_kmer_kmc1(kmer);
			else
				counters[counters_pos++] = count_for_kmer_kmc1(kmer_rev);
		}
	}
	if (counters_pos < counters.size())
	{
		fill_n(counters.begin() + counters_pos, counters.size() - counters_pos, 0);
		counters_pos = static_cast<uint32>(counters.size());
	}
	return true;
}

//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
bool CKMCFile::GetCountersForRead_kmc1(const std::string& read, std::vector<uint32>& counters)
{	
	uint32 read_len = static_cast<uint32>(read.length());
	counters.resize(read.length() - kmer_length + 1);
	std::string transformed_read = read;
	for (char& c : transformed_read)
		c = CKmerAPI::num_codes[(uchar)c];

	uint32 i = 0;
	CKmerAPI kmer(kmer_length);
	uint32 pos = 0;
	
	uint32 counters_pos = 0;

	while (i + kmer_length - 1 < read_len)
	{
		bool contains_N = false;
		while (i < read_len && pos < kmer_length)
		{
			if (CKmerAPI::num_codes[(uchar)read[i]] < 0)
			{
				pos = 0;
				kmer.clear();
				++i;
				uint32 wrong_kmers = MIN(i - counters_pos, static_cast<uint32>(counters.size()) - counters_pos);
				fill_n(counters.begin() + counters_pos, wrong_kmers, 0);
				counters_pos += wrong_kmers;
				contains_N = true;
				break;
			}
			else
				kmer.insert2bits(pos++, CKmerAPI::num_codes[(uchar)read[i++]]);
		}
		if (contains_N)
			continue;
		if (pos == kmer_length)
		{			
			counters[counters_pos++] = count_for_kmer_kmc1(kmer);
		}
		else
			break;

		while (i < read_len)
		{
			if (CKmerAPI::num_codes[(uchar)read[i]] < 0)
			{
				pos = 0;
				break;
			}
			kmer.SHL_insert2bits(CKmerAPI::num_codes[(uchar)read[i++]]);
			counters[counters_pos++] = count_for_kmer_kmc1(kmer);
		}
	}
	if (counters_pos < counters.size())
	{
		fill_n(counters.begin() + counters_pos, counters.size() - counters_pos, 0);
		counters_pos = static_cast<uint32>(counters.size());
	}
	return true;
}

//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
void CKMCFile::GetSuperKmers(const std::string& transformed_read, super_kmers_t& super_kmers)
{
	uint32 i = 0;
	uint32 len = 0; //length of super k-mer
	uint32 signature_start_pos;
	CMmer current_signature(signature_len), end_mmer(signature_len);

	while (i + kmer_length - 1 < transformed_read.length())
	{
		bool contains_N = false;
		//building first signature after 'N' or at the read beginning
		for (uint32 j = 0; j < signature_len; ++j, ++i)
		{
			if (transformed_read[i] < 0)//'N'
			{
				contains_N = true;
				break;
			}
		}
		//signature must be shorter than k-mer so if signature contains 'N', k-mer will contains it also
		if (contains_N)
		{
			++i;
			continue;
		}
		len = signature_len;
		signature_start_pos = i - signature_len;
		current_signature.insert(transformed_read.c_str() + signature_start_pos);
		end_mmer.set(current_signature);

		for (; i < transformed_read.length(); ++i)
		{
			if (transformed_read[i] < 0)//'N'
			{
				if (len >= kmer_length)
				{
					super_kmers.push_back(std::make_tuple(i - len, len, signature_map[current_signature.get()]));
				}
				len = 0;
				++i;
				break;
			}
			end_mmer.insert(transformed_read[i]);
			if (end_mmer < current_signature)//signature at the end of current k-mer is lower than current
			{
				if (len >= kmer_length)
				{
					super_kmers.push_back(std::make_tuple(i - len, len, signature_map[current_signature.get()]));
					len = kmer_length - 1;
				}
				current_signature.set(end_mmer);
				signature_start_pos = i - signature_len + 1;
			}
			else if (end_mmer == current_signature)
			{
				current_signature.set(end_mmer);
				signature_start_pos = i - signature_len + 1;
			}
			else if (signature_start_pos + kmer_length - 1 < i)//need to find new signature
			{
				super_kmers.push_back(std::make_tuple(i - len, len, signature_map[current_signature.get()]));
				len = kmer_length - 1;
				//looking for new signature
				++signature_start_pos;
				//building first signature in current k-mer
				end_mmer.insert(transformed_read.c_str() + signature_start_pos);
				current_signature.set(end_mmer);
				for (uint32 j = signature_start_pos + signature_len; j <= i; ++j)
				{
					end_mmer.insert(transformed_read[j]);
					if (end_mmer <= current_signature)
					{
						current_signature.set(end_mmer);
						signature_start_pos = j - signature_len + 1;
					}
				}
			}
			++len;
		}
	}
	if (len >= kmer_length)//last one in read
	{
		super_kmers.push_back(std::make_tuple(i - len, len, signature_map[current_signature.get()]));
	}
}

//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
bool CKMCFile::GetCountersForRead_kmc2_both_strands(const std::string& read, std::vector<uint32>& counters)
{
	counters.resize(read.length() - kmer_length + 1);
	std::string transformed_read = read;
	for (char& c : transformed_read)
		c = CKmerAPI::num_codes[(uchar)c];

	super_kmers_t super_kmers;
	GetSuperKmers(transformed_read, super_kmers);

	uint32 counters_pos = 0;
	if (super_kmers.empty())
	{
		fill_n(counters.begin(), counters.size(), 0);
		return true;
	}

	CKmerAPI kmer(kmer_length), rev_kmer(kmer_length);

	uint32 last_end = 0;

	//'N' somewhere in first k-mer
	if (std::get<0>(super_kmers.front()) > 0)
	{
		fill_n(counters.begin(), std::get<0>(super_kmers.front()), 0);
		last_end = std::get<0>(super_kmers.front());
		counters_pos = std::get<0>(super_kmers.front());
	}
	for (auto& super_kmer : super_kmers)
	{
		//'N's between super k-mers
		if (last_end < std::get<0>(super_kmer))
		{
			uint32 gap = std::get<0>(super_kmer) -last_end;
			fill_n(counters.begin() + counters_pos, kmer_length + gap - 1, 0);
			counters_pos += kmer_length + gap - 1;
		}
		last_end = std::get<0>(super_kmer) +std::get<1>(super_kmer);

		kmer.from_binary(transformed_read.c_str() + std::get<0>(super_kmer));
		rev_kmer.from_binary_rev(transformed_read.c_str() + std::get<0>(super_kmer));

		uint32 bin_start_pos = std::get<2>(super_kmer) * single_LUT_size;
		if(kmer < rev_kmer)
			counters[counters_pos++] = count_for_kmer_kmc2(kmer, bin_start_pos);
		else
			counters[counters_pos++] = count_for_kmer_kmc2(rev_kmer, bin_start_pos);

		for (uint32 i = std::get<0>(super_kmer) +kmer_length; i < std::get<0>(super_kmer) +std::get<1>(super_kmer); ++i)
		{
			kmer.SHL_insert2bits(transformed_read[i]);
			rev_kmer.SHR_insert2bits(3 - transformed_read[i]);
			if(kmer < rev_kmer)
				counters[counters_pos++] = count_for_kmer_kmc2(kmer, bin_start_pos);
			else
				counters[counters_pos++] = count_for_kmer_kmc2(rev_kmer, bin_start_pos);
		}
	}
	//'N's at the end of read
	if (counters_pos < counters.size())
	{
		fill_n(counters.begin() + counters_pos, counters.size() - counters_pos, 0);
		counters_pos = static_cast<uint32>(counters.size());
	}

	return true;
}


//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
bool CKMCFile::GetCountersForRead_kmc2(const std::string& read, std::vector<uint32>& counters)
{	
	counters.resize(read.length() - kmer_length + 1);
	std::string transformed_read = read;
	for (char& c : transformed_read)
		c = CKmerAPI::num_codes[(uchar)c];
	
	super_kmers_t super_kmers;
	GetSuperKmers(transformed_read, super_kmers);
	
	uint32 counters_pos = 0;
	if (super_kmers.empty())
	{
		fill_n(counters.begin(), counters.size(), 0);
		return true;
	}

	CKmerAPI kmer(kmer_length);

	uint32 last_end = 0;

	//'N' somewhere in first k-mer
	if (std::get<0>(super_kmers.front()) > 0)
	{
		fill_n(counters.begin(), std::get<0>(super_kmers.front()), 0);
		last_end = std::get<0>(super_kmers.front());
		counters_pos = std::get<0>(super_kmers.front());
	}
	for (auto& super_kmer : super_kmers)
	{
		//'N's between super k-mers
		if (last_end < std::get<0>(super_kmer))
		{
			uint32 gap = std::get<0>(super_kmer) -last_end;
			fill_n(counters.begin() + counters_pos, kmer_length + gap - 1, 0);
			counters_pos += kmer_length + gap - 1;
		}
		last_end = std::get<0>(super_kmer) + std::get<1>(super_kmer);
		
		kmer.from_binary(transformed_read.c_str() + std::get<0>(super_kmer));

		uint32 bin_start_pos = std::get<2>(super_kmer) * single_LUT_size;
		counters[counters_pos++] = count_for_kmer_kmc2(kmer, bin_start_pos);

		for (uint32 i = std::get<0>(super_kmer) +kmer_length; i < std::get<0>(super_kmer) +std::get<1>(super_kmer); ++i)
		{
			kmer.SHL_insert2bits(transformed_read[i]);
			counters[counters_pos++] = count_for_kmer_kmc2(kmer, bin_start_pos);
		}
	}
	//'N's at the end of read
	if (counters_pos < counters.size())
	{
		fill_n(counters.begin() + counters_pos, counters.size() - counters_pos, 0);
		counters_pos = static_cast<uint32>(counters.size());
	}

	return true;
}


//---------------------------------------------------------------------------------
// Auxiliary function.
//---------------------------------------------------------------------------------
bool CKMCFile::BinarySearch(int64 index_start, int64 index_stop, const CKmerAPI& kmer, uint64& counter, uint32 pattern_offset)
{
	if (index_start >= static_cast<int64>(total_kmers))
		return false;
	uchar *sufix_byte_ptr = nullptr;
	uint64 sufix = 0;

	//sufix_offset is always 56
	uint32 sufix_offset = 56;			// the offset of a sufix is for shifting the sufix towards MSB, to compare the sufix with a pattern
	// Bytes of a pattern to search are always shifted towards MSB

	uint32 row_index = 0;				// the number of a current row in an array kmer_data

	bool found = false;

	while (index_start <= index_stop)
	{
		int64 mid_index = (index_start + index_stop) / 2;
		sufix_byte_ptr = &sufix_file_buf[mid_index * sufix_rec_size];

		uint64 pattern = 0;

		pattern_offset = (lut_prefix_length + kmer.byte_alignment) * 2;

		row_index = 0;
		for (uint32 a = 0; a < sufix_size; a++)		//check byte by byte
		{
			pattern = kmer.kmer_data[row_index];
			pattern = pattern << pattern_offset;
			pattern = pattern & 0xff00000000000000;

			sufix = sufix_byte_ptr[a];
			sufix = sufix << sufix_offset;

			if (pattern != sufix)
				break;

			pattern_offset += 8;

			if (pattern_offset == 64)				//the end of a word
			{
				pattern_offset = 0;
				row_index++;
			}
		}

		if (pattern == sufix)
		{
			found = true;
			break;
		}
		if (sufix < pattern)
			index_start = mid_index + 1;
		else
			index_stop = mid_index - 1;
	}

	if (found)
	{
		if (counter_size == 0)
			counter = 1;
		else
		{
			sufix_byte_ptr += sufix_size;
			counter = *sufix_byte_ptr;

			for (uint32 b = 1; b < counter_size; b++)
			{
				uint64 aux = 0x000000ff & *(sufix_byte_ptr + b);

				aux = aux << 8 * (b);
				counter = aux | counter;
			}
		}
		//applay filtering only if counter_size != 0
		return (counter_size == 0) || ((counter >= min_count) && (counter <= max_count));
	}
	return false;
}


// ***** EOF
