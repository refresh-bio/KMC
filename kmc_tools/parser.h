/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _PARSER_H
#define _PARSER_H
#include "defs.h"
#include "expression_node.h"
#include "tokenizer.h"
#include "output_parser.h"

#include <iostream>
#include <fstream>
#include <regex>
#include <map>
#include <list>


//************************************************************************************************************
// CParser - parser for complex operations
//************************************************************************************************************
class CParser
{	
	std::ifstream file;
	uint32 line_no;
	std::regex input_line_pattern;
	std::regex output_line_pattern;
	std::regex empty_line_pattern;	
	std::map<std::string, uint32> input;	
	void parseInputLine(const std::string& line);

	template<unsigned SIZE>
	CExpressionNode<SIZE>* parseOutputLine(const std::string& line);
	void parseOtuputParamsLine();
	bool nextLine(std::string& line);
	CConfig& config;

public:
	CParser(const std::string& src);
	void ParseInputs();

	template<unsigned SIZE>
	CExpressionNode<SIZE>* ParseOutput();	

};

//************************************************************************************************************
template<unsigned SIZE> CExpressionNode<SIZE>* CParser::ParseOutput()
{
	std::string line;
	if (!nextLine(line) || line.find("OUTPUT_PARAMS:") != std::string::npos)
	{
		std::cout << "Error: None output was defined\n";
		exit(1);
	}

	auto result = parseOutputLine<SIZE>(line);

	while (nextLine(line))
	{
		if (line.find("OUTPUT_PARAMS:") != std::string::npos)
		{
			parseOtuputParamsLine();
			break;
		}
	}

	return result;
}

//************************************************************************************************************
template<unsigned SIZE> CExpressionNode<SIZE>* CParser::parseOutputLine(const std::string& line)
{
	std::smatch match;
	if (std::regex_search(line, match, output_line_pattern))
	{
#ifdef ENABLE_DEBUG
		std::cout << "out file name " << match[1] << "\n";
		std::cout << "rest of output " << match[2] << "\n";

		std::cout << "Tokenize resf of output\n";
#endif
		config.output_desc.file_src = match[1];

		//trim whitespaces at the end
		static const std::string whitespace = " \t\r\n\v\f";
		auto end = config.output_desc.file_src.find_last_not_of(whitespace);
		config.output_desc.file_src.erase(end + 1);
		if (config.output_desc.file_src == "")
		{
			std::cout << "Error: wrong line format, line: " << line_no << " (output file name is not specified)\n";
			exit(1);
		}

		CTokenizer tokenizer;

		std::list<Token> tokens;
		tokenizer.Tokenize(match[2], tokens);

		COutputParser<SIZE> out_parser(tokens, input);
		return out_parser.Parse();

	}
	else
	{
		std::cout << "Error: wrong line format, line: " << line_no << "\n";
		exit(1);
	}
}
#endif

// ***** EOF