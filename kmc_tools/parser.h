/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
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

	void parseOutputLine(const std::string& line);
	void parseOtuputParamsLine();
	bool nextLine(std::string& line);
	CConfig& config;
	
	CTokenizer tokenizer;
	std::list<Token> tokens;

public:
	CParser(const std::string& src);
	void ParseInputs();
	void ParseOutput();

	template<unsigned SIZE>
	CExpressionNode<SIZE>* GetExpressionRoot();

};

//************************************************************************************************************
template<unsigned SIZE> CExpressionNode<SIZE>* CParser::GetExpressionRoot()
{
	COutputParser<SIZE> out_parser(tokens, input);
	return out_parser.Parse();
}

#endif

// ***** EOF