/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#include "stdafx.h"
#include "parameters_parser.h"
#include <iostream>
using namespace std;


uint32 CParametersParser::replace_zero(uint32 val, const char* param_name, uint32 value_to_set_if_zero)
{
	if (val == 0)
	{
		cout << "Warning: min value for " << param_name << " is " << value_to_set_if_zero << ". Your value will be converted to " << value_to_set_if_zero << "\n";
		return value_to_set_if_zero;
	}
	return val;
}


void CParametersParser::parse_int_or_float(bool& force_float, bool& force_int, float& float_value, uint32& int_val, const char* param_name)
{
	if (strchr(argv[pos] + 3, '.'))
	{
		float_value = (float)atof(argv[pos++] + 3);
		if (float_value > 1.0f || float_value < 0.0f)
		{
			cout << " Error: wrong value for fastq input parameter: "<< param_name <<"\n";
			exit(1);
		}
		if (force_int)
		{
			cout << "Error: both -ci, -cx must be specified as real number [0;1] or as integer \n";
			exit(1);
		}
		force_float = true;
		config.filtering_params.use_float_value = true;
	}
	else
	{
		int_val = atoi(argv[pos++] + 3);
		if (force_float)
		{
			cout << "Error: both -ci, -cx must be specified as real number [0;1] or as integer \n";
			exit(1);
		}
		force_int = true;
		config.filtering_params.use_float_value = false;
	}
}

void CParametersParser::parse_global_params()
{
	//defaults
	config.avaiable_threads = thread::hardware_concurrency();

	//override defaults if specified
	for( ; pos < argc && argv[pos][0] == '-' ; ++pos)
	{
		if (strncmp(argv[pos], "-t", 2) == 0)
		{
			config.avaiable_threads = atoi(argv[pos] + 2);
			continue;
		}
		if (argv[pos][1] == 'v')
		{
			config.verbose = true;
			continue;
		}
		if (strncmp(argv[pos], "-hp", 3) == 0)
		{
			config.percent_progress.Hide();
			continue;
		}
	}
}

void CParametersParser::read_input_fastq_desc()
{
	if (pos >= argc)
	{
		cout << "Error: Input fastq files(s) missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cout << "Error: Input fastq file(s) required, but " << argv[pos] << " found\n";
		exit(1);
	}
	string input_file_name = argv[pos++];
	if (input_file_name[0] != '@')
		config.filtering_params.input_srcs.push_back(input_file_name);
	else
	{
		ifstream in(input_file_name.c_str() + 1);
		if (!in.good())
		{
			cout << "Error: No " << input_file_name.c_str() + 1 << " file\n";
			exit(1);
		}
		string s;
		while (getline(in, s))
		{
			if (s != "")
				config.filtering_params.input_srcs.push_back(s);
		}
		in.close();
	}

	bool force_float = false;
	bool force_int = false;
	
	for (int i = 0; i < 3 && pos < argc; ++i)
	{
		if(argv[pos][0] != '-')
			break;
		if (strncmp(argv[pos], "-ci", 3) == 0)
		{
			parse_int_or_float(force_float, force_int, config.filtering_params.f_min_kmers, config.filtering_params.n_min_kmers, "-ci");
		}
		else if (strncmp(argv[pos], "-cx", 3) == 0)
		{
			parse_int_or_float(force_float, force_int, config.filtering_params.f_max_kmers, config.filtering_params.n_max_kmers, "-cx");
		}
		else if (strncmp(argv[pos], "-f", 2) == 0)
		{
			switch (argv[pos++][2])
			{
			case 'a':
				config.filtering_params.input_file_type = CFilteringParams::file_type::fasta;
				break;
			case 'q':
				config.filtering_params.input_file_type = CFilteringParams::file_type::fastq;
				break;
			default:
				cout << "Error: unknow parameter " << argv[pos - 1] << "\n";
				exit(1);
				break;
			}
		}
	}
	config.filtering_params.output_file_type = config.filtering_params.input_file_type;
}


void CParametersParser::read_output_fastq_desc()
{
	if (pos >= argc)
	{
		cout << "Error: Output fastq source missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cout << "Error: Output fastq source required, but " << argv[pos] << "found\n";
		exit(1);
	}
	config.filtering_params.output_src = argv[pos++];

	while (pos < argc && argv[pos][0] == '-')
	{
		if (strncmp(argv[pos], "-f", 2) == 0)
		{
			switch (argv[pos][2])
			{
			case 'q':
				config.filtering_params.output_file_type = CFilteringParams::file_type::fastq;
				break;
			case 'a':
				config.filtering_params.output_file_type = CFilteringParams::file_type::fasta;
				break;
			default:
				cout << "Error: unknown parameter " << argv[pos] << "\n";
				exit(1);
				break;
			}
			if (config.filtering_params.input_file_type == CFilteringParams::file_type::fasta && config.filtering_params.output_file_type == CFilteringParams::file_type::fastq)
			{
				cout << "Error: cannot set -fq for output when -fa is set for input\n";
				exit(1);
			}
		}
		else
		{
			cout << "Error: Unknown parameter: " << argv[pos] << "\n";
			exit(1);
		}
		++pos;
	}
}


void CParametersParser::read_dump_params()
{
	while (pos < argc && argv[pos][0] == '-')
	{
		if (strncmp(argv[pos], "-s", 2) == 0)
		{
			config.dump_params.sorted_output = true;
		}
		else
		{
			cout << "Warning: Unknow parameter for dump operation: " << argv[pos] << "\n";			
		}
		++pos;
	}

}

void CParametersParser::read_input_desc()
{
	if (pos >= argc)
	{
		cout << "Error: Input database source missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cout << "Error: Input database source required, but " << argv[pos] << "found\n";
		exit(1);
	}
	CInputDesc desc(argv[pos++]);
	config.input_desc.push_back(desc);
	for (int i = 0; i < 2 && pos < argc; ++i)
	{
		if (strncmp(argv[pos], "-", 1) != 0)
			break;
		if (strncmp(argv[pos], "-ci", 3) == 0)
		{
			config.input_desc.back().cutoff_min = replace_zero(atoi(argv[pos++] + 3), "-ci", 1);
		}
		else if (strncmp(argv[pos], "-cx", 3) == 0)
		{
			config.input_desc.back().cutoff_max = replace_zero(atoi(argv[pos++] + 3), "-cx", 1);
		}
		else
		{
			cout << "Error: Unknow parameter: " << argv[pos];
			exit(1);
		}
	}
}

void CParametersParser::read_output_desc()
{
	if (pos >= argc)
	{
		cout << "Error: Output database source missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cout << "Error: Output database source required, but " << argv[pos] << "found\n";
		exit(1);
	}
	config.output_desc.file_src = argv[pos++];
	for (int i = 0; i < 2 && pos < argc; ++i)
	{
		if (strncmp(argv[pos], "-", 1) != 0)
			break;
		if (strncmp(argv[pos], "-ci", 3) == 0)
		{
			config.output_desc.cutoff_min = replace_zero(atoi(argv[pos++] + 3), "-ci", 1);
		}
		else if (strncmp(argv[pos], "-cx", 3) == 0)
		{
			config.output_desc.cutoff_max = replace_zero(atoi(argv[pos++] + 3), "-cx", 1);
		}
		else if (strncmp(argv[pos], "-cs", 3) == 0)
		{
			config.output_desc.counter_max = replace_zero(atoi(argv[pos++] + 3), "-cs", 1); 
		}
		else
		{
			cout << "Error: Unknow parameter: " << argv[pos];
			exit(1);
		}
	}
}

void CParametersParser::Usage()
{
	CUsageDisplayerFactory disp(CConfig::GetInstance().mode);
	disp.GetUsageDisplayer().Display();
}

CParametersParser::CParametersParser(int argc, char** argv) :argc(argc), argv(argv), config(CConfig::GetInstance())
{
	pos = 0;
	if (argc < 2)
	{
		Usage();
		exit(1);
	}
}

void CParametersParser::Parse()
{
	pos = 1;
	parse_global_params();
	if (strcmp(argv[pos], "intersect") == 0)
	{
		config.mode = CConfig::Mode::INTERSECTION;
	}
	else if (strcmp(argv[pos], "kmers_subtract") == 0)
	{
		config.mode = CConfig::Mode::KMERS_SUBTRACT;
	}
	else if (strcmp(argv[pos], "counters_subtract") == 0)
	{
		config.mode = CConfig::Mode::COUNTERS_SUBTRACT;
	}
	else if (strcmp(argv[pos], "union") == 0)
	{
		config.mode = CConfig::Mode::UNION;
	}
	else if (strcmp(argv[pos], "complex") == 0)
	{
		config.mode = CConfig::Mode::COMPLEX;
	}
	else if (strcmp(argv[pos], "sort") == 0)
	{
		config.mode = CConfig::Mode::SORT;
	}
	else if (strcmp(argv[pos], "reduce") == 0)
	{
		config.mode = CConfig::Mode::REDUCE;
	}
	else if (strcmp(argv[pos], "compact") == 0)
	{
		config.mode = CConfig::Mode::COMPACT;
	}
	else if (strcmp(argv[pos], "histogram") == 0)
	{
		config.mode = CConfig::Mode::HISTOGRAM;
	}
	else if (strcmp(argv[pos], "dump") == 0)
	{
		config.mode = CConfig::Mode::DUMP;
	}
	else if (strcmp(argv[pos], "compare") == 0)
	{
		config.mode = CConfig::Mode::COMPARE;
	}
	else if (strcmp(argv[pos], "filter") == 0)
	{
		config.mode = CConfig::Mode::FILTER;
	}
	else
	{
		cout << "Error: Unknow mode: " << argv[pos] << "\n";
		Usage();
		exit(1);
	}

	if (argc == 2)
	{
		Usage();
		exit(1);
	}

	pos++;
	if (config.mode == CConfig::Mode::INTERSECTION || config.mode == CConfig::Mode::KMERS_SUBTRACT || config.mode == CConfig::Mode::UNION || config.mode == CConfig::Mode::COUNTERS_SUBTRACT)
	{
		read_input_desc(); //first input
		read_input_desc(); //second input
		read_output_desc(); //output
	}
	else if (config.mode == CConfig::Mode::FILTER)
	{
		read_input_desc(); //kmc db
		read_input_fastq_desc(); //fastq input
		read_output_fastq_desc();
	}
	else if (config.mode == CConfig::Mode::COMPLEX)
	{
		if (strncmp(argv[2], "-", 1) == 0)
		{
			cout << "Error: operations description file expected but " << argv[2] << " found\n";
			exit(1);
		}		
		complex_parser = make_unique<CParser>(argv[pos]);
		complex_parser->ParseInputs();
	}
	else if (config.mode == CConfig::Mode::DUMP)
	{
		read_dump_params();
		read_input_desc();
		read_output_desc();

	}
	else if (config.mode == CConfig::Mode::SORT || config.mode == CConfig::Mode::HISTOGRAM || config.mode == CConfig::Mode::REDUCE || config.mode == CConfig::Mode::COMPACT)
	{
		read_input_desc();
		read_output_desc();
		if (config.mode == CConfig::Mode::COMPACT)
		{
			if (config.output_desc.counter_max)
				cout << "Warning: -cs can not be specified for compact operation, value specified will be ignored\n";
			config.output_desc.counter_max = 1;
		}
	}
	else if (config.mode == CConfig::Mode::COMPARE)
	{
		read_input_desc();
		read_input_desc();
	}
}

bool CParametersParser::validate_input_dbs()
{
	config.headers.push_back(CKMC_header(config.input_desc.front().file_src));

	uint32 kmer_len = config.headers.front().kmer_len;
	uint32 mode = config.headers.front().mode;
	if (mode == 1)
	{
		cout << "Error: quality counters are not supported in kmc tools\n"; 
		return false;
	}
	for (uint32 i = 1; i < config.input_desc.size(); ++i)
	{
		config.headers.push_back(CKMC_header(config.input_desc[i].file_src));
		CKMC_header& h = config.headers.back();
		if (h.mode != mode)
		{
			cout << "Error: quality/direct based counters conflict!\n"; 
			return false;
		}
		if (h.kmer_len != kmer_len)
		{
			cout << "Database " << config.input_desc.front().file_src << " contains " << kmer_len << "-mers, but database " << config.input_desc[i].file_src << " contains " << h.kmer_len << "-mers\n";
			return false;
		}
	}
	config.kmer_len = kmer_len;


	//update cutoff_min and coutoff_max if it was not set with parameters
	for (uint32 i = 0; i < config.input_desc.size(); ++i)
	{
		if (config.input_desc[i].cutoff_min == 0)
			config.input_desc[i].cutoff_min = config.headers[i].min_count;
		if (config.input_desc[i].cutoff_max == 0)
			config.input_desc[i].cutoff_max = config.headers[i].max_count;
	}

	//update output description if it was not set with parameters
	if (config.output_desc.cutoff_min == 0)
	{
		uint32 min_cutoff_min = config.input_desc.front().cutoff_min;
		for (uint32 i = 0; i < config.input_desc.size(); ++i)
		{
			if (config.input_desc[i].cutoff_min < min_cutoff_min)
				min_cutoff_min = config.input_desc[i].cutoff_min;
		}
		config.output_desc.cutoff_min = min_cutoff_min;
		if (config.verbose)
			cout << "-ci was not specified for output. It will be set to " << min_cutoff_min << "\n";
	}

	if (config.output_desc.cutoff_max == 0)
	{
		if (config.mode == CConfig::Mode::HISTOGRAM) //for histogram default value differs
		{
			config.output_desc.cutoff_max = MIN(config.headers.front().max_count, MIN(HISTOGRAM_MAX_COUNTER_DEFAULT, (uint32)((1ull << (8 * config.headers.front().counter_size)) - 1)));
		}
		else
		{
			uint32 max_cutoff_max = config.input_desc.front().cutoff_max;
			for (uint32 i = 0; i < config.input_desc.size(); ++i)
			{
				if (config.input_desc[i].cutoff_max > max_cutoff_max)
					max_cutoff_max = config.input_desc[i].cutoff_max;
			}
			config.output_desc.cutoff_max = max_cutoff_max;
		}
		
		if (config.verbose)
			cout << "-cx was not specified for output. It will be set to " << config.output_desc.cutoff_max << "\n";
	}
	if (config.output_desc.counter_max == 0)
	{
		uint32 max_counter_max = config.headers.front().counter_size;
		for (uint32 i = 0; i < config.headers.size(); ++i)
		{
			if (config.headers[i].counter_size> max_counter_max)
				max_counter_max = config.headers[i].counter_size;
		}

		max_counter_max = (uint32)((1ull << (max_counter_max << 3)) - 1);
		config.output_desc.counter_max = max_counter_max;
		if (config.verbose)
			cout << "-cs was not specified for output. It will be set to " << max_counter_max << "\n";
	}
	return true;
}
	
void CParametersParser::SetThreads()
{
	uint32 threads_left = config.avaiable_threads;
	//threads distribution: as many as possible for kmc2 database input, 1 thread for main thread which make operations calculation
	vector<reference_wrapper<CInputDesc>> kmc2_desc;
	
	if (!config.Is1ArgOper()) 
		threads_left = MAX(1, threads_left - 1);
	
	for (uint32 i = 0; i < config.headers.size(); ++i)
	{
		if (config.headers[i].IsKMC2())
		{
			kmc2_desc.push_back(ref(config.input_desc[i]));
		}
	}
	if (kmc2_desc.size())
	{
		uint32 per_signle_kmc2_input = MAX(1, (uint32)(threads_left / kmc2_desc.size()));
		uint32 per_last_kmc2_input = MAX(1, (uint32)((threads_left + kmc2_desc.size() - 1) / kmc2_desc.size()));

		for (uint32 i = 0; i < kmc2_desc.size() - 1; ++i)
			kmc2_desc[i].get().threads = per_signle_kmc2_input;

		kmc2_desc.back().get().threads = per_last_kmc2_input;
	}
}



// ***** EOF