/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz and Agnieszka Debudaj-Grabysz

  Version: 2.1
  Date   : 2014-12-18
*/

#ifndef _KMER_API_H
#define _KMER_API_H


#include "kmer_defs.h"
#include <string>
#include <iostream>
#include "mmer.h"
class CKMCFile;

class CKmerAPI 
{
protected:

	uint64 *kmer_data;				// An array to store kmer's data. On 64 bits 32 symbols can be stored
									// Data are shifted to let sufix's symbols to start with a border of a byte

	
	uint32 kmer_length;				// Kmer's length, in symbols
	uchar byte_alignment;			// A number of "empty" symbols placed before prefix to let sufix's symbols to start with a border of a byte

	uint32 no_of_rows;				// A number of 64-bits words allocated for kmer_data 	

	friend class CKMCFile;

public:
	static const char char_codes[];
	static char num_codes[256];
	struct _si  
    {   _si()
        {
            for(int i = 0; i < 256; i++)
                num_codes[i] = -1;
			num_codes['A'] = num_codes['a'] = 0;
			num_codes['C'] = num_codes['c'] = 1;
			num_codes['G'] = num_codes['g'] = 2;
			num_codes['T'] = num_codes['t'] = 3;
        }
    } static _init;


// ----------------------------------------------------------------------------------
// The constructor creates kmer for the number of symbols equal to length. 
// The array kmer_data has the size of ceil((length + byte_alignment) / 32))
// IN	: length - a number of symbols of a kmer
// ----------------------------------------------------------------------------------
	inline CKmerAPI(uint32 length = 0)
	{
		if(length)
		{
			if(length % 4)
				byte_alignment = 4 - (length % 4);	
			else
				byte_alignment = 0;

			no_of_rows = (((length + byte_alignment) % 32) ? (length + byte_alignment) / 32 + 1 : (length + byte_alignment) / 32); 
			//no_of_rows = (int)ceil((double)(length + byte_alignment) / 32);
			kmer_data = new uint64[no_of_rows];

			memset(kmer_data, 0, sizeof(*kmer_data) * no_of_rows);
		}
		else
		{
			kmer_data = NULL;
			no_of_rows = 0;
			byte_alignment = 0;
		}
		kmer_length = length;
	};
//-----------------------------------------------------------------------
// The destructor
//-----------------------------------------------------------------------
	inline ~CKmerAPI()
	{
		if (kmer_data != NULL)
			delete [] kmer_data;
	};

//-----------------------------------------------------------------------
// The copy constructor
//-----------------------------------------------------------------------
	inline CKmerAPI(const CKmerAPI &kmer)
	{
		kmer_length = kmer.kmer_length;
		byte_alignment = kmer.byte_alignment;
		no_of_rows = kmer.no_of_rows;
		
		kmer_data = new uint64[no_of_rows];
			
		for(uint32 i = 0; i < no_of_rows; i++)
			kmer_data[i] = kmer.kmer_data[i];

	};

//-----------------------------------------------------------------------
// The operator =
//-----------------------------------------------------------------------	
	inline CKmerAPI& operator=(const CKmerAPI &kmer)
	{
		if(kmer.kmer_length != kmer_length)		
		{
			if(kmer_length && kmer_data)
				delete [] kmer_data;
		
			kmer_length = kmer.kmer_length;
			byte_alignment = kmer.byte_alignment;
			no_of_rows = kmer.no_of_rows;
		
			kmer_data = new uint64[no_of_rows];
		}

		for(uint32 i = 0; i < no_of_rows; i++)
			kmer_data[i] = kmer.kmer_data[i];

		return *this;
	};

//-----------------------------------------------------------------------
// The operator ==
//-----------------------------------------------------------------------
	inline bool operator==(const CKmerAPI &kmer)
	{
			if(kmer.kmer_length != kmer_length)
				return false;

			for(uint32 i = 0; i < no_of_rows; i++)
				if(kmer.kmer_data[i] != kmer_data[i])
					return false;

			return true;

	};

//-----------------------------------------------------------------------
// Operator < . If arguments differ in length a result is undefined
//-----------------------------------------------------------------------
	inline bool operator<(const CKmerAPI &kmer)
	{
			if(kmer.kmer_length != kmer_length)
				return false;					

			for(uint32 i = 0; i < no_of_rows; i++)
				if(kmer.kmer_data[i] > kmer_data[i])
					return true;
				else
					if(kmer.kmer_data[i] < kmer_data[i])
						return false;
				
			return false;
	};

//-----------------------------------------------------------------------
// Return a symbol of a kmer from an indicated position (numbered form 0).
// The symbol is returned as an ASCI character A/C/G/T
// IN	: pos - a position of a symbol
// RET	: symbol - a symbol placed on a position pos
//-----------------------------------------------------------------------
	inline char get_asci_symbol(unsigned int pos)
	{
		if(pos >= kmer_length)
			return 0;
		
		uint32 current_row = (pos + byte_alignment) / 32;
		uint32 current_pos = ((pos + byte_alignment) % 32) * 2;
		uint64 mask = 0xc000000000000000 >> current_pos;
		uint64 symbol = kmer_data[current_row] & mask;
		symbol = symbol >> (64 - current_pos - 2);
		return char_codes[symbol];
	
	};

	//-----------------------------------------------------------------------
	// Return a symbol of a kmer from an indicated position (numbered form 0)
	// The symbol is returned as a numerical value 0/1/2/3
	// IN	: pos - a position of a symbol
	// RET	: symbol - a symbol placed on a position pos
	//-----------------------------------------------------------------------
	inline uchar get_num_symbol(unsigned int pos)
	{
		if (pos >= kmer_length)
			return 0;

		uint32 current_row = (pos + byte_alignment) / 32;
		uint32 current_pos = ((pos + byte_alignment) % 32) * 2;
		uint64 mask = 0xc000000000000000 >> current_pos;
		uint64 symbol = kmer_data[current_row] & mask;
		symbol = symbol >> (64 - current_pos - 2);
		uchar* byte_ptr = reinterpret_cast<uchar*>(&symbol);
		return *byte_ptr;

	};

//-----------------------------------------------------------------------
// Convert kmer into string (an alphabet ACGT)
// RET	: string kmer
//-----------------------------------------------------------------------
	inline std::string to_string()
	{
		std::string string_kmer;
		uchar *byte_ptr;
		uchar c;
		uchar temp_byte_alignment = byte_alignment;
		uint32 cur_string_size = 0;

		string_kmer.reserve(kmer_length + 1);
		string_kmer.resize(kmer_length + 1);
		

		for(uint32 row_counter = 0; row_counter < no_of_rows; row_counter++)
		{
			byte_ptr = reinterpret_cast<uchar*>(&kmer_data[row_counter]);

			byte_ptr += 7;					// shift a pointer towards a MSB
		
			for(uint32 i = 0; (i < kmer_length) && (i < 32); i += 4)		// 32 symbols of any "row" in kmer_data
			{				
				if((i == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0xc0 & *byte_ptr;			//11000000
					c = c >> 6;
					string_kmer[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}

				if((i == 0) && temp_byte_alignment)				// check if a  byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0x30 & *byte_ptr;			//00110000
					c = c >> 4;
					string_kmer[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}
				
				if((i == 0) && temp_byte_alignment)				// check if a  byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0x0c & *byte_ptr;			//00001100
					c = c >> 2;
					string_kmer[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}
																// no need to check byte alignment as its length is at most 3 
				c = 0x03 & *byte_ptr;			//00000011
				string_kmer[cur_string_size++] = char_codes[c];
				if (cur_string_size == kmer_length) break;

				byte_ptr--;
			}
		}
		string_kmer[cur_string_size] = '\0';
		return string_kmer;
	};
//-----------------------------------------------------------------------
// Convert kmer into string (an alphabet ACGT). The function assumes enough memory was allocated
// OUT	: str - string kmer. 
//-----------------------------------------------------------------------
	inline void to_string(char *str)
	{
		uchar *byte_ptr;
		uchar c;
		uchar temp_byte_alignment = byte_alignment;
		uchar cur_string_size = 0;
		

		for(uint32 row_counter = 0; row_counter < no_of_rows; row_counter++)
		{
			byte_ptr = reinterpret_cast<uchar*>(&kmer_data[row_counter]);

			byte_ptr += 7;					// shift a pointer towards a MSB
		
			for(uint32 i = 0; (i < kmer_length) && (i < 32); i += 4)		// 32 symbols of any "row" in kmer_data
			{
				
				if((i == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0xc0 & *byte_ptr;			//11000000
					c = c >> 6;
					str[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}

				if((i == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0x30 & *byte_ptr;			//00110000
					c = c >> 4;
					str[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}
				
				if((i == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0x0c & *byte_ptr;			//00001100
					c = c >> 2;
					str[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}
																// no need to check byte alignment as its length is at most 3  
				c = 0x03 & *byte_ptr;			//00000011
				str[cur_string_size++] = char_codes[c];
				if (cur_string_size == kmer_length) break;

				byte_ptr--;
			}
		}
		str[cur_string_size] = '\0';
	};

//-----------------------------------------------------------------------
// Convert kmer into string (an alphabet ACGT)
// OUT 	: str - string kmer
//-----------------------------------------------------------------------
	 inline void to_string(std::string &str) 
	 {
		uchar *byte_ptr;
		uchar c;
		uchar temp_byte_alignment = byte_alignment;
		uint32 cur_string_size = 0;

		str.reserve(kmer_length + 1);
		str.resize(kmer_length + 1);
		

		for(uint32 row_counter = 0; row_counter < no_of_rows; row_counter++)
		{
			byte_ptr = reinterpret_cast<uchar*>(&kmer_data[row_counter]);

			byte_ptr += 7;					// shift a pointer towards a MSB
		
			for(uint32 i = 0; (i < kmer_length) && (i < 32); i += 4)		// 32 symbols of any "row" in kmer_data
			{
				
				if((i == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0xc0 & *byte_ptr;			//11000000
					c = c >> 6;
					str[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}

				if((i == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0x30 & *byte_ptr;			//00110000
					c = c >> 4;
					str[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}
				
				if((i == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
					temp_byte_alignment--;						
				else
				{
					c = 0x0c & *byte_ptr;			//00001100
					c = c >> 2;
					str[cur_string_size++] = char_codes[c];
					if (cur_string_size == kmer_length) break;
				}
																// no need to check byte alignment as its length is at most 3 
				c = 0x03 & *byte_ptr;			//00000011
				str[cur_string_size++] = char_codes[c];
				if (cur_string_size == kmer_length) break;

				byte_ptr--;
			}
		}
		str[cur_string_size] = '\0';
	 };

//-----------------------------------------------------------------------
// Convert a string of an alphabet ACGT into a kmer of a CKmerAPI
// IN	: kmer_string	- a string of an alphabet ACGT
// RET	: true			- if succesfull
//-----------------------------------------------------------------------
	 inline bool from_string(std::string kmer_string)
	 {
		 unsigned char c_char;
		 uchar c_binary;
		 uchar temp_byte_alignment;

		 for (uint32 ii = 0; ii < kmer_string.size(); ++ii)
		 {
		 if (num_codes[(uchar)kmer_string[ii]] == -1)
			 return false;
		 }
		

		if(kmer_length != kmer_string.size())	
		{
			if(kmer_length && kmer_data)
				delete [] kmer_data;
		
			kmer_length = (uint32)kmer_string.size();

			if(kmer_length % 4)
				byte_alignment = 4 - (kmer_length % 4);	
			else
				byte_alignment = 0;

			
			if(kmer_length != 0)
			{
				no_of_rows = (((kmer_length + byte_alignment) % 32) ? (kmer_length + byte_alignment) / 32 + 1 : (kmer_length + byte_alignment) / 32);
				//no_of_rows = (int)ceil((double)(kmer_length + byte_alignment) / 32);
				kmer_data = new uint64[no_of_rows];
				//memset(kmer_data, 0, sizeof(*kmer_data) * no_of_rows);
			}
		}

		memset(kmer_data, 0, sizeof(*kmer_data) * no_of_rows);
		temp_byte_alignment = byte_alignment;
		uint32 i = 0;
		uint32 i_in_string = 0;
		uchar *byte_ptr;

		for(uint32 row_index = 0; row_index < no_of_rows; row_index++)
		{
			byte_ptr = reinterpret_cast<uchar*>(&kmer_data[row_index]);
			byte_ptr += 7;					// shift a pointer towards a MSB

			while(i < kmer_length)
			{
				if((i_in_string == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
				{															
					temp_byte_alignment--;						
					i++;
				}
				else
				{
					c_char = kmer_string[i_in_string];
					c_binary = num_codes[c_char];
					c_binary = c_binary << 6;		//11000000
					*byte_ptr = *byte_ptr | c_binary;
					i++;
					i_in_string++;
					if(i_in_string == kmer_length) break;
				}

				if((i_in_string == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
				{
					temp_byte_alignment--;						
					i++;
				}
				else
				{
					c_char = kmer_string[i_in_string];
					c_binary = num_codes[c_char];
					c_binary = c_binary << 4;
					*byte_ptr = *byte_ptr | c_binary;
					i++;
					i_in_string++;
					if(i_in_string == kmer_length) break;
				}

				//!!!if((i == 0) && temp_byte_alignment)	//poprawka zg3oszona przez Maaka D3ugosza			// check if a byte_alignment placed before a prefix is to be skipped
				if ((i_in_string == 0) && temp_byte_alignment)				// check if a byte_alignment placed before a prefix is to be skipped
				{
					temp_byte_alignment--;						
					i++;
				}
				else
				{
					c_char = kmer_string[i_in_string];
					c_binary = num_codes[c_char];
					c_binary = c_binary << 2;
					*byte_ptr = *byte_ptr | c_binary;
					i++;
					i_in_string++;
					if(i_in_string == kmer_length) break;
				}

				c_char = kmer_string[i_in_string];
				c_binary = num_codes[c_char];
				*byte_ptr = *byte_ptr | c_binary;
				i++;
				i_in_string++;
				if(i_in_string == kmer_length) break;

				if( i % 32 == 0) 
					break; //check if a new "row" is to be started
				byte_ptr--;
			}
		};
		return true;
	}
//-----------------------------------------------------------------------
// Counts a signature of an existing kmer
// IN	: sig_len	- the length of a signature
// RET	: signature value
//-----------------------------------------------------------------------
	 uint32 get_signature(uint32 sig_len)
	 {
		 uchar symb;
		 CMmer cur_mmr(sig_len);
		 
		 for(uint32 i = 0; i < sig_len; ++i)
		 {
			 symb = get_num_symbol(i);
			 cur_mmr.insert(symb);
		 }
		 CMmer min_mmr(cur_mmr);
		 for (uint32 i = sig_len; i < kmer_length; ++i)
		 {
			 symb = get_num_symbol(i);
			 cur_mmr.insert(symb);
			 
			 if (cur_mmr < min_mmr)
				 min_mmr = cur_mmr;
		 }
		 return min_mmr.get();
	 }
	
};


#endif

// ***** EOF
