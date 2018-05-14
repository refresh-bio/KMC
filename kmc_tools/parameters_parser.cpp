/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include "parameters_parser.h"
#include <iostream>
using namespace std;


uint32 CParametersParser::replace_zero(uint32 val, const char* param_name, uint32 value_to_set_if_zero)
{
	if (val == 0)
	{
		cerr << "Warning: min value for " << param_name << " is " << value_to_set_if_zero << ". Your value will be converted to " << value_to_set_if_zero << "\n";
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
			cerr << "Error: wrong value for fastq input parameter: "<< param_name <<"\n";
			exit(1);
		}
		if (force_int)
		{
			cerr << "Error: both -ci, -cx must be specified as real number [0;1] or as integer \n";
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
			cerr << "Error: both -ci, -cx must be specified as real number [0;1] or as integer \n";
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
			if (strlen(argv[pos]) < 3)
			{
				std::cerr << "Error: -t require value\n";
				exit(1);
			}
			config.avaiable_threads = atoi(argv[pos] + 2);
			continue;
		}
		else if (argv[pos][1] == 'v')
		{
			config.verbose = true;
			continue;
		}
		else if (strncmp(argv[pos], "-hp", 3) == 0)
		{
			config.percent_progress.Hide();
			continue;
		}
		else
		{
			std::cerr << "Error: unknown global option " << argv[pos] << "\n";
			exit(1);
		}
	}
}

void CParametersParser::read_input_fastq_desc()
{
	if (pos >= argc)
	{
		cerr << "Error: Input fastq files(s) missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cerr << "Error: Input fastq file(s) required, but " << argv[pos] << " found\n";
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
			cerr << "Error: No " << input_file_name.c_str() + 1 << " file\n";
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
				cerr << "Error: unknow parameter " << argv[pos - 1] << "\n";
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
		cerr << "Error: Output fastq source missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cerr << "Error: Output fastq source required, but " << argv[pos] << "found\n";
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
				cerr << "Error: unknown parameter " << argv[pos] << "\n";
				exit(1);
				break;
			}
			if (config.filtering_params.input_file_type == CFilteringParams::file_type::fasta && config.filtering_params.output_file_type == CFilteringParams::file_type::fastq)
			{
				cerr << "Error: cannot set -fq for output when -fa is set for input\n";
				exit(1);
			}
		}
		else
		{
			cerr << "Error: Unknown parameter: " << argv[pos] << "\n";
			exit(1);
		}
		++pos;
	}
}

void CParametersParser::read_filter_params()
{
	while (pos < argc && argv[pos][0] == '-')
	{
		if (strncmp(argv[pos], "-t", 2) == 0)
		{
			config.filtering_params.filter_mode = CFilteringParams::FilterMode::trim;
		}
		else if (strncmp(argv[pos], "-hm", 3) == 0)
		{
			config.filtering_params.filter_mode = CFilteringParams::FilterMode::hard_mask;
		}
		else
		{
			cerr << "Warning: Unknow parameter for filter operation: " << argv[pos] << "\n";
		}
		++pos;
	}
}

void CParametersParser::read_check_params()
{
	if (pos >= argc)
	{
		std::cerr << "Error: check operation require k-mer to check\n";
		exit(1);
	}
	config.check_params.kmer = argv[pos++];
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
			cerr << "Warning: Unknow parameter for dump operation: " << argv[pos] << "\n";			
		}
		++pos;
	}

}


void CParametersParser::read_operation_type()
{
	if (strncmp(argv[pos], "-oc", 3) == 0)
	{
		char* mode = argv[pos] + 3;
		if (strcmp(mode, "min") == 0)
		{
			config.counter_op_type = CounterOpType::MIN;
		}
		else if (strcmp(mode, "max") == 0)
		{
			config.counter_op_type = CounterOpType::MAX;
		}
		else if (strcmp(mode, "sum") == 0)
		{
			config.counter_op_type = CounterOpType::SUM;
		}
		else if (strcmp(mode, "diff") == 0)
		{
			config.counter_op_type = CounterOpType::DIFF;
		}
		else if (strcmp(mode, "left") == 0)
		{
			config.counter_op_type = CounterOpType::FROM_DB1;
		}
		else if (strcmp(mode, "right") == 0)
		{
			config.counter_op_type = CounterOpType::FROM_DB2;
		}
		else
		{
			cerr << "Error: unknown counter calculation mode: " << mode << ". Allowed values: min, max, sum, diff, left, right\n";
			exit(1);
		}
		++pos;

		if (config.mode == CConfig::Mode::KMERS_SUBTRACT)
		{
			cerr << "Error: counter calculation mode not allowed for kmers_subtract\n";
			exit(1);
		}
	}
}

void CParametersParser::read_input_desc()
{
	if (pos >= argc)
	{
		cerr << "Error: Input database source missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cerr << "Error: Input database source required, but " << argv[pos] << "found\n";
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
			cerr << "Error: Unknow parameter: " << argv[pos];
			exit(1);
		}
	}
}

bool CParametersParser::read_output_for_transform()
{
	if (pos >= argc)
		return false;

	CTransformOutputDesc::OpType op_type;
	uint64 counter_value = 0; //for set_counts only
	if (strcmp(argv[pos], "sort") == 0)
	{
		op_type = CTransformOutputDesc::OpType::SORT;	
	}
	else if (strcmp(argv[pos], "reduce") == 0)
	{
		op_type = CTransformOutputDesc::OpType::REDUCE;
	}
	else if (strcmp(argv[pos], "compact") == 0)
	{
		op_type = CTransformOutputDesc::OpType::COMPACT;
	}
	else if (strcmp(argv[pos], "histogram") == 0)
	{
		op_type = CTransformOutputDesc::OpType::HISTOGRAM;
	}
	else if (strcmp(argv[pos], "dump") == 0)
	{
		op_type = CTransformOutputDesc::OpType::DUMP;
	}
	else if (strcmp(argv[pos], "set_counts") == 0)
	{
		op_type = CTransformOutputDesc::OpType::SET_COUNTS;
		++pos;
		if (pos >= argc)
		{
			cerr << "Error: set_counts operation requires count value\n";
			Usage();
			exit(1);
		}
		
		if (strncmp(argv[pos], "-", 1) == 0)
		{
			cerr << "Error: Count value expected, but " << argv[pos] << " found\n";
			exit(1);
		}

		char* ptr_end = nullptr;
		counter_value = strtoull(argv[pos], &ptr_end, 10);
		auto len = strlen(argv[pos]);
		if (argv[pos] + len != ptr_end)
		{
			cerr << "Error: Count value expected, but " << argv[pos] << " found\n";
			exit(1);
		}		

		if (counter_value > numeric_limits<uint32>::max())
		{
			cerr << "Error: currenlty kmc_tools supports counter values up to " << numeric_limits<uint32>::max() << "\n";
			exit(1);
		}

	}
	else
	{
		cerr << "Error: unknown operation: " << argv[pos] << "\n";
		Usage();
		exit(1);
	}

	config.transform_output_desc.emplace_back(op_type);

	++pos;

	//read op paramts
	while (pos < argc)
	{
		if (strncmp(argv[pos], "-", 1) != 0)
			break;
		if (strncmp(argv[pos], "-s", 3) == 0)
		{
			if (config.transform_output_desc.back().op_type == CTransformOutputDesc::OpType::DUMP)
			{
				config.transform_output_desc.back().sorted_output = true;
			}
			else
			{
				cerr << "Error: -s parameter allowed only for dump operation\n";
				Usage();
				exit(1);
			}
		}
		else
		{
			cerr << "Error: unknown operation parameter: " << argv[pos] <<"\n";
			exit(1);
		}
		++pos;
	}

	if (pos >= argc)
	{
		cerr << "Error: Output database path missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cerr << "Error: Output database path required, but " << argv[pos] << " found\n";
		exit(1);
	}

	config.transform_output_desc.back().file_src = argv[pos++];

	if (op_type == CTransformOutputDesc::OpType::SET_COUNTS)
		config.transform_output_desc.back().counter_value = counter_value;

	//read database params
	while (pos < argc)
	{
		if (strncmp(argv[pos], "-", 1) != 0)
			break;
		if (strncmp(argv[pos], "-ci", 3) == 0)
		{
			config.transform_output_desc.back().cutoff_min = replace_zero(atoi(argv[pos++] + 3), "-ci", 1);
		}
		else if (strncmp(argv[pos], "-cx", 3) == 0)
		{
			config.transform_output_desc.back().cutoff_max = replace_zero(atoi(argv[pos++] + 3), "-cx", 1);
		}
		else if (strncmp(argv[pos], "-cs", 3) == 0)
		{
			config.transform_output_desc.back().counter_max = replace_zero(atoi(argv[pos++] + 3), "-cs", 1);
		}
		else
		{
			cerr << "Error: Unknown parameter: " << argv[pos];
			Usage();
			exit(1);
		}
	}
	if (op_type == CTransformOutputDesc::OpType::COMPACT)
	{
		if (config.transform_output_desc.back().counter_max)
			cerr << "Warning: -cs can not be specified for compact operation, value specified will be ignored\n";
		config.transform_output_desc.back().counter_max = 1;
	}	
	if (op_type == CTransformOutputDesc::OpType::SET_COUNTS)
	{
		auto& tmp = config.transform_output_desc.back();
		if (tmp.counter_max || tmp.cutoff_max || tmp.cutoff_min)
			cerr << "Warning: -cs, -cx, -ci cannot be specified for set_counts operation, values will be ignored\n";
		tmp.counter_max = tmp.cutoff_max = numeric_limits<uint32>::max();
		tmp.cutoff_min = 1;
	}

	return true;
}

bool CParametersParser::read_output_desc_for_simple()
{
	if (pos >= argc)
		return false;

	//get op name
	CSimpleOutputDesc::OpType op_type;
	if (strcmp(argv[pos], "intersect") == 0)
	{
		op_type = CSimpleOutputDesc::OpType::INTERSECT;
	}
	else if (strcmp(argv[pos], "kmers_subtract") == 0)
	{
		op_type = CSimpleOutputDesc::OpType::KMERS_SUBTRACTION;
	}
	else if (strcmp(argv[pos], "counters_subtract") == 0)
	{
		op_type = CSimpleOutputDesc::OpType::COUNTERS_SUBTRACTION;
	}
	else if (strcmp(argv[pos], "union") == 0)
	{
		op_type = CSimpleOutputDesc::OpType::UNION;
	}
	else if (strcmp(argv[pos], "reverse_kmers_subtract") == 0)
	{
		op_type = CSimpleOutputDesc::OpType::REVERSE_KMERS_SUBTRACTION;
	}
	else if (strcmp(argv[pos], "reverse_counters_subtract") == 0)
	{
		op_type = CSimpleOutputDesc::OpType::REVERSE_COUNTERS_SUBTRACTION;
	}
	else
	{
		cerr << "Error: unknown operation: " << argv[pos] << "\n";
		Usage();
		exit(1);
	}

	++pos;
	if (pos >= argc)
	{
		cerr << "Error: Output database path missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cerr << "Error: Output database path required, but " << argv[pos] << " found\n";
		exit(1);
	}
	config.simple_output_desc.emplace_back(op_type);
	config.simple_output_desc.back().file_src = argv[pos++];

	while (pos < argc)
	{
		if (strncmp(argv[pos], "-", 1) != 0)
			break;
		if (strncmp(argv[pos], "-ci", 3) == 0)
		{
			config.simple_output_desc.back().cutoff_min = replace_zero(atoi(argv[pos++] + 3), "-ci", 1);
		}
		else if (strncmp(argv[pos], "-cx", 3) == 0)
		{
			config.simple_output_desc.back().cutoff_max = replace_zero(atoi(argv[pos++] + 3), "-cx", 1);
		}
		else if (strncmp(argv[pos], "-cs", 3) == 0)
		{
			config.simple_output_desc.back().counter_max = replace_zero(atoi(argv[pos++] + 3), "-cs", 1);
		}
		else if (strncmp(argv[pos], "-oc", 3) == 0)
		{
			if (op_type == CSimpleOutputDesc::OpType::KMERS_SUBTRACTION || op_type == CSimpleOutputDesc::OpType::REVERSE_KMERS_SUBTRACTION)
			{
				cerr << "Error: -oc not allowed for kmers_subtract and reverse_kmers_subtract as it doesn't make sense (equal k-mers form both input will not be present in output)\n";
				exit(1);
			}
			char* mode = argv[pos] + 3;
			if (strcmp(mode, "min") == 0)
			{
				config.simple_output_desc.back().counter_op = CounterOpType::MIN;
			}
			else if (strcmp(mode, "max") == 0)
			{
				config.simple_output_desc.back().counter_op = CounterOpType::MAX;
			}
			else if (strcmp(mode, "sum") == 0)
			{
				config.simple_output_desc.back().counter_op = CounterOpType::SUM;
			}
			else if (strcmp(mode, "diff") == 0)
			{
				config.simple_output_desc.back().counter_op = CounterOpType::DIFF;
			}
			else if (strcmp(mode, "left") == 0)
			{
				config.simple_output_desc.back().counter_op = CounterOpType::FROM_DB1;
			}
			else if (strcmp(mode, "right") == 0)
			{
				config.simple_output_desc.back().counter_op = CounterOpType::FROM_DB2;
			}
			else
			{
				cerr << "Error: unknown counter calculation mode: " << mode << ". Allowed values: min, max, sum, diff, left, right\n";
				exit(1);
			}
			++pos;
		}
		else
		{
			cerr << "Error: Unknow parameter: " << argv[pos];
			Usage();
			exit(1);
		}

	}
	return true;
}

void CParametersParser::read_output_desc()
{
	if (pos >= argc)
	{
		cerr << "Error: Output database source missed\n";
		exit(1);
	}
	if (strncmp(argv[pos], "-", 1) == 0)
	{
		cerr << "Error: Output database source required, but " << argv[pos] << "found\n";
		exit(1);
	}
	
	config.output_desc.file_src = argv[pos++];
	for (int i = 0; i < 3 && pos < argc; ++i)
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
			cerr << "Error: Unknow parameter: " << argv[pos];
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
		config.counter_op_type = CounterOpType::MIN;
	}
	else if (strcmp(argv[pos], "kmers_subtract") == 0)
	{
		config.mode = CConfig::Mode::KMERS_SUBTRACT;
	}
	else if (strcmp(argv[pos], "counters_subtract") == 0)
	{
		config.mode = CConfig::Mode::COUNTERS_SUBTRACT;
		config.counter_op_type = CounterOpType::DIFF;
	}
	else if (strcmp(argv[pos], "union") == 0)
	{
		config.mode = CConfig::Mode::UNION;
		config.counter_op_type = CounterOpType::SUM;
	}
	else if (strcmp(argv[pos], "complex") == 0)
	{
		config.mode = CConfig::Mode::COMPLEX;
	}
	else if (strcmp(argv[pos], "transform") == 0)
	{
		config.mode = CConfig::Mode::TRANSFORM;
	}
	else if (strcmp(argv[pos], "simple") == 0)
	{
		config.mode = CConfig::Mode::SIMPLE_SET;
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
	else if (strcmp(argv[pos], "info") == 0)
	{
		config.mode = CConfig::Mode::INFO;
	}
	else if (strcmp(argv[pos], "check") == 0)
	{
		config.mode = CConfig::Mode::CHECK;
	}
	else
	{
		cerr << "Error: Unknow mode: " << argv[pos] << "\n";
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
		read_operation_type();
		read_input_desc(); //first input
		read_input_desc(); //second input
		read_output_desc(); //output
	}
	else if (config.mode == CConfig::Mode::FILTER)
	{
		read_filter_params();
		read_input_desc(); //kmc db
		read_input_fastq_desc(); //fastq input
		read_output_fastq_desc();
		if (config.filtering_params.use_float_value && config.filtering_params.filter_mode != CFilteringParams::FilterMode::normal)
		{
			cerr << "Error: trim (-t) and soft mask (-hm) are not compatibile with float values of cut off (-ci -cx)\n";
			exit(1);
		}
	}
	else if (config.mode == CConfig::Mode::COMPLEX)
	{
		if (strncmp(argv[pos], "-", 1) == 0)
		{
			cerr << "Error: operations description file expected but " << argv[2] << " found\n";
			exit(1);
		}		
		complex_parser = make_unique<CParser>(argv[pos]);
		complex_parser->ParseInputs();
		complex_parser->ParseOutput();
	}
	else if (config.mode == CConfig::Mode::TRANSFORM)
	{
		read_input_desc();
		while (read_output_for_transform())
			;

		if (config.transform_output_desc.size() == 0)
		{
			cerr << "Error: output missed\n";
			Usage();
			exit(1);
		}

	}
	else if (config.mode == CConfig::Mode::SIMPLE_SET)
	{
		read_input_desc(); //first input
		read_input_desc(); //second input
		while (read_output_desc_for_simple())
			;
		if (config.simple_output_desc.size() == 0)
		{
			cerr << "Error: output missed\n";
			Usage();
			exit(1);
		}
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
				cerr << "Warning: -cs can not be specified for compact operation, value specified will be ignored\n";
			config.output_desc.counter_max = 1;
		}
		if (config.mode == CConfig::Mode::HISTOGRAM)
		{
			if (config.output_desc.cutoff_max)
				std::cerr << "Warning: -cx not allowed for histogram output, value specified will be ignored\n";
			if(config.output_desc.cutoff_min)
				std::cerr << "Warning: -ci not allowed for histogram output, value specified will be ignored\n";
			if (config.output_desc.counter_max)
				std::cerr << "Warning: -cs not allowed for histogram output, value specified will be ignored\n";


			config.output_desc.cutoff_max = config.input_desc.front().cutoff_max;
			config.output_desc.cutoff_min = config.input_desc.front().cutoff_min;
			config.output_desc.counter_max = 0;

		}
	}
	else if (config.mode == CConfig::Mode::INFO)
	{
		read_input_desc();
	}
	else if (config.mode == CConfig::Mode::CHECK)
	{
		read_input_desc();
		read_check_params();		
	}
	else if (config.mode == CConfig::Mode::COMPARE)
	{
		read_input_desc();
		read_input_desc();
	}
}

uint32 CParametersParser::get_min_cutoff_min()
{
	uint32 min_cutoff_min = config.input_desc.front().cutoff_min;
	for (uint32 i = 0; i < config.input_desc.size(); ++i)
	{
		if (config.input_desc[i].cutoff_min < min_cutoff_min)
			min_cutoff_min = config.input_desc[i].cutoff_min;
	}
	return min_cutoff_min;
}

uint32 CParametersParser::get_max_cutoff_max()
{
	uint32 max_cutoff_max = config.input_desc.front().cutoff_max;
	for (uint32 i = 0; i < config.input_desc.size(); ++i)
	{
		if (config.input_desc[i].cutoff_max > max_cutoff_max)
			max_cutoff_max = config.input_desc[i].cutoff_max;
	}
	return max_cutoff_max;
}

uint32 CParametersParser::get_max_counter_max()
{
	uint32 max_counter_max = config.headers.front().counter_size;
	for (uint32 i = 0; i < config.headers.size(); ++i)
	{
		if (config.headers[i].counter_size> max_counter_max)
			max_counter_max = config.headers[i].counter_size;
	}
	max_counter_max = (uint32)((1ull << (max_counter_max << 3)) - 1);
	return max_counter_max;
}

bool CParametersParser::validate_input_dbs()
{
	config.headers.push_back(CKMC_header(config.input_desc.front().file_src));

	uint32 kmer_len = config.headers.front().kmer_len;
	uint32 mode = config.headers.front().mode;
	if (mode == 1)
	{
		cerr << "Error: quality counters are not supported in kmc tools\n"; 
		return false;
	}
	for (uint32 i = 1; i < config.input_desc.size(); ++i)
	{
		config.headers.push_back(CKMC_header(config.input_desc[i].file_src));
		CKMC_header& h = config.headers.back();
		if (h.mode != mode)
		{
			cerr << "Error: quality/direct based counters conflict!\n";
			return false;
		}
		if (h.kmer_len != kmer_len)
		{
			cerr << "Database " << config.input_desc.front().file_src << " contains " << kmer_len << "-mers, but database " << config.input_desc[i].file_src << " contains " << h.kmer_len << "-mers\n";
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
	if (config.mode == CConfig::Mode::SIMPLE_SET)
	{
		uint32 min_cutoff_min = get_min_cutoff_min();
		uint32 max_cutoff_max = get_max_cutoff_max();
		uint32 max_counter_max = get_max_counter_max();
		
		for (auto& desc : config.simple_output_desc)
		{
			if (desc.cutoff_min == 0)
				desc.cutoff_min = min_cutoff_min;
			if (desc.cutoff_max == 0)
				desc.cutoff_max = max_cutoff_max;
			if (desc.counter_max == 0)
				desc.counter_max = max_counter_max;
		}
	}
	else if (config.mode == CConfig::Mode::TRANSFORM)
	{
		uint32 min_cutoff_min = get_min_cutoff_min();
		uint32 max_cutoff_max = get_max_cutoff_max();
		uint32 max_counter_max = get_max_counter_max();

		for (auto& desc : config.transform_output_desc)
		{
			if (desc.op_type == CTransformOutputDesc::OpType::SET_COUNTS)
				continue;

			if (desc.cutoff_min == 0)
				desc.cutoff_min = min_cutoff_min;
			if (desc.cutoff_max == 0)
			{
				if (desc.op_type == CTransformOutputDesc::OpType::HISTOGRAM)//for histogram default value differs				
				{
					desc.cutoff_max = MIN(config.headers.front().max_count, MIN(HISTOGRAM_MAX_COUNTER_DEFAULT, (uint32)((1ull << (8 * config.headers.front().counter_size)) - 1)));					
				}
				else
					desc.cutoff_max = max_cutoff_max;
			}
			if (desc.counter_max == 0)
				desc.counter_max = max_counter_max;
		}
	}
	else //old style operations
	{
		if (config.output_desc.cutoff_min == 0)
		{
			uint32 min_cutoff_min = get_min_cutoff_min();
			config.output_desc.cutoff_min = min_cutoff_min;
			if (config.verbose)
				cerr << "Warning: -ci was not specified for output. It will be set to " << min_cutoff_min << "\n";
		}

		if (config.output_desc.cutoff_max == 0)
		{
			if (config.mode == CConfig::Mode::HISTOGRAM) //for histogram default value differs
				config.output_desc.cutoff_max = MIN(config.headers.front().max_count, MIN(HISTOGRAM_MAX_COUNTER_DEFAULT, (uint32)((1ull << (8 * config.headers.front().counter_size)) - 1)));
			else
			{
				uint32 max_cutoff_max = get_max_cutoff_max();
				config.output_desc.cutoff_max = max_cutoff_max;
			}
			if (config.verbose)
				cerr << "Warning: -cx was not specified for output. It will be set to " << config.output_desc.cutoff_max << "\n";
		}
		if (config.output_desc.counter_max == 0)
		{
			uint32 max_counter_max = get_max_counter_max();			
			config.output_desc.counter_max = max_counter_max;
			if (config.verbose)
				cerr << "Warning: -cs was not specified for output. It will be set to " << max_counter_max << "\n";
		}
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
		else
		{
			if (threads_left > 2)
				threads_left -= 2; //2 threads for kmc1 input
			else
				threads_left = 1;
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