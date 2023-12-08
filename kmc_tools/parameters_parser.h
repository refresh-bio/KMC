/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.2.3
  Date   : 2023-12-08
*/

#ifndef _PARAMETERS_PARSER_H
#define _PARAMETERS_PARSER_H

#include "defs.h"
#include "parser.h"
#include <memory>
class CParametersParser
{	
	std::unique_ptr<CParser> complex_parser;
	int argc;
	char** argv;
	int pos;
	CConfig& config;

	uint32 replace_zero(uint32 val, const char* param_name, uint32 value_to_set_if_zero);
	void parse_int_or_float(bool& force_float, bool& force_int, float& float_value, uint32& int_val, const char* param_name);
	void parse_global_params();
	void read_input_fastq_desc();
	void read_output_fastq_desc();
	void read_input_desc();
	void read_check_params();
	void read_filter_params();
	bool read_output_desc_for_simple();
	bool read_output_for_transform();

	uint32 get_max_counter_max();
	uint64 get_max_cutoff_max();
	uint32 get_min_cutoff_min();
public:
	
	CParametersParser(int argc, char** argv);
	void Usage();
	
	template<unsigned SIZE>
	CExpressionNode<SIZE>* GetExpressionRoot();

	void Parse();
	bool validate_input_dbs();
	void SetThreads();
	
};

template<unsigned SIZE>
CExpressionNode<SIZE>* CParametersParser::GetExpressionRoot()
{
	return complex_parser->GetExpressionRoot<SIZE>();
}
#endif


// ***** EOF