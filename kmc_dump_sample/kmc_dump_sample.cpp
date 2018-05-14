/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  This file demonstrates the example usage of kmc_api software. 
  It reads kmer_counter's output and prints kmers to an output file.

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include <iostream>
#include "../kmc_api/kmc_file.h"

void print_info(void);

//----------------------------------------------------------------------------------
// Check if --help or --version was used
bool help_or_version(int argc, char** argv)
{
	const std::string version = "--version";
	const std::string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

int _tmain(int argc, char* argv[])
{
	if (argc == 1 || help_or_version(argc, argv))
	{
		print_info();
		return 0;
	}

	CKMCFile kmer_data_base;
	int32 i;
	uint32 min_count_to_set = 0;
	uint32 max_count_to_set = 0;
	std::string input_file_name;
	std::string output_file_name;

	FILE * out_file;
	//------------------------------------------------------------
	// Parse input parameters
	//------------------------------------------------------------
	if(argc < 3)
	{
		print_info();
		return EXIT_FAILURE;
	}

	for(i = 1; i < argc; ++i)
	{
		if(argv[i][0] == '-')
		{	
			if(strncmp(argv[i], "-ci", 3) == 0)
				min_count_to_set = atoi(&argv[i][3]);
			else if(strncmp(argv[i], "-cx", 3) == 0)
					max_count_to_set = atoi(&argv[i][3]);
		}
		else
			break;
	}

	if(argc - i < 2)
	{ 
		print_info();
		return EXIT_FAILURE;
	}

	input_file_name = std::string(argv[i++]);
	output_file_name = std::string(argv[i]);

	if((out_file = fopen (output_file_name.c_str(),"wb")) == NULL)
	{
		print_info();
		return EXIT_FAILURE;
	}

	setvbuf(out_file, NULL ,_IOFBF, 1 << 24);

	//------------------------------------------------------------------------------
	// Open kmer database for listing and print kmers within min_count and max_count
	//------------------------------------------------------------------------------

	if (!kmer_data_base.OpenForListing(input_file_name))
	{
		print_info();
		return EXIT_FAILURE ;
	}
	else
	{
		uint32 _kmer_length;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		uint64 _total_kmers;

		kmer_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

		CKmerAPI kmer_object(_kmer_length);
		
		if(min_count_to_set)
		if (!(kmer_data_base.SetMinCount(min_count_to_set)))
				return EXIT_FAILURE;
		if(max_count_to_set)
		if (!(kmer_data_base.SetMaxCount(max_count_to_set)))
				return EXIT_FAILURE;	

		
		std::string str;
		if (_mode) //quake compatible mode
		{
			float counter;			
			while (kmer_data_base.ReadNextKmer(kmer_object, counter))
			{
				kmer_object.to_string(str);	
				fprintf(out_file, "%s\t%f\n", str.c_str(), counter);			
			}
		}
		else 
		{
			uint32 counter;			
			while (kmer_data_base.ReadNextKmer(kmer_object, counter))
			{
				kmer_object.to_string(str);
				fprintf(out_file, "%s\t%u\n", str.c_str(), counter);
			}
		}
		
	
		fclose(out_file);
		kmer_data_base.Close();
	}

	return EXIT_SUCCESS; 
}
// -------------------------------------------------------------------------
// Print execution options 
// -------------------------------------------------------------------------
void print_info(void)
{
	std::cout << "KMC dump ver. " << KMC_VER << " (" << KMC_DATE << ")\n"
			  << "\nUsage:\nkmc_dump [options] <kmc_database> <output_file>\n"
			  << "Parameters:\n"
			  << "<kmc_database> - kmer_counter's output\n"
			  << "Options:\n"
			  << "-ci<value> - exclude k-mers occurring less than <value> times\n"
			  << "-cx<value> - exclude k-mers occurring more of than <value> times\n";
};

// ***** EOF
