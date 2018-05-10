/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include "parser.h"
#include "tokenizer.h"
#include "output_parser.h"
#include "config.h"

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

CParser::CParser(const std::string& src):
	config(CConfig::GetInstance())
{
	line_no = 0;
	file.open(src);
	if (!file.is_open())
	{
		std::cerr << "Error: cannot open file: " << src << "\n";
		exit(1);
	}
	//input_line_pattern = "\\s*(\\w*)\\s*=\\s*(.*)$";
	input_line_pattern = "^\\s*([\\w-+]*)\\s*=\\s*(.*)$";
	output_line_pattern = "^\\s*(.*)\\s*=\\s*(.*)$"; //TODO: consider valid file name	
	empty_line_pattern = "^\\s*$";
}

/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

void CParser::ParseInputs()
{
	std::string line;
	while (true)
	{
		if (!nextLine(line))
		{
			std::cerr << "Error: 'INPUT:' missing\n";
			exit(1);
		}
		if (line.find("INPUT:") != std::string::npos)
			break;
	}

	if (!nextLine(line) || line.find("OUTPUT:") != std::string::npos)
	{
		std::cerr << "Error: None input was defined\n";
		exit(1);
	}

	while (true)
	{
		parseInputLine(line);
		if (!nextLine(line))
		{
			std::cerr << "Error: 'OUTPUT:' missing\n";
			exit(1);
		}
		if (line.find("OUTPUT:") != std::string::npos)
			break;
	}
}
//************************************************************************************************************
 void CParser::ParseOutput()
 {
	std::string line;
	if (!nextLine(line) || line.find("OUTPUT_PARAMS:") != std::string::npos)
	{
 		std::cerr << "Error: None output was defined\n";
 		exit(1);
	}
 
	parseOutputLine(line);
 
	while (nextLine(line))
	{
 		if (line.find("OUTPUT_PARAMS:") != std::string::npos)
 		{
 			parseOtuputParamsLine();
 			break;
 		}
	}
 }

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
void CParser::parseInputLine(const std::string& line)
{
	std::smatch match;
	if (std::regex_search(line, match, input_line_pattern))
	{
#ifdef ENABLE_DEBUG
		std::cout << "\ninput name: " << match[1];
		std::cout << "\nafter = " << match[2];
#endif
		if (input.find(match[1]) != input.end())
		{
			std::cerr << "Error: Name redefinition(" << match[1] << ")" << " line: " << line_no << "\n";
			exit(1);
		}
		if (CTokenizer::GetKeywords().find(match[1]) != CTokenizer::GetKeywords().end())
		{
			std::cerr << "Error: `" << match[1] << "` is not valid name, line: " << line_no << "\n";
			exit(1);
		}
		else
		{
			std::string file_name;
			std::istringstream stream(match[2]);
			
			CInputDesc desc;

			if (!(stream >> desc.file_src))
			{
				std::cerr << "Error: file name for " << match[1] << " was not specified, line: "<< line_no <<"\n";
				exit(1);
			}
			std::string tmp;
			while (stream >> tmp)
			{
				if (strncmp(tmp.c_str(), "-ci", 3) == 0)
				{
					desc.cutoff_min = atoi(tmp.c_str() + 3);
					continue;
				}
				else if (strncmp(tmp.c_str(), "-cx", 3) == 0)
				{
					desc.cutoff_max = atoi(tmp.c_str() + 3);
					continue;
				}
				std::cerr << "Error: Unknow parameter " << tmp << " for variable " << match[1] << ", line: "<< line_no <<"\n";
				exit(1);
			}

			config.input_desc.push_back(desc);
			input[match[1]] = (uint32)(config.input_desc.size() - 1);				
		}
	}
	else
	{
		std::cerr << "Error: wrong line format, line: " << line_no << "\n";
		exit(1);
	}
}


//************************************************************************************************************
void CParser::parseOutputLine(const std::string& line)
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
			std::cerr << "Error: wrong line format, line: " << line_no << " (output file name is not specified)\n";
			exit(1);
		}

		tokenizer.Tokenize(match[2], tokens);
	}
	else
	{
		std::cerr << "Error: wrong line format, line: " << line_no << "\n";
		exit(1);
	}
}

/*****************************************************************************************************************************/
void CParser::parseOtuputParamsLine()
{
	std::string line;

	if (!nextLine(line))
	{
		std::cerr << "Warning: OUTPUT_PARAMS exists, but no parameters are defined\n";
	}
	else
	{
		std::istringstream stream(line);
		std::string tmp;
		while (stream >> tmp)
		{
			if (strncmp(tmp.c_str(), "-ci", 3) == 0)
			{				
				config.output_desc.cutoff_min = atoi(tmp.c_str() + 3);
				continue;
			}
			else if (strncmp(tmp.c_str(), "-cx", 3) == 0)
			{			
				config.output_desc.cutoff_max = atoi(tmp.c_str() + 3);
				continue;
			}
			else if ((strncmp(tmp.c_str(), "-cs", 3) == 0))
			{
				config.output_desc.counter_max = atoi(tmp.c_str() + 3);
				continue;
			}
			std::cerr << "Error: Unknow parameter " << tmp << " for variable " << tmp << ", line: " << line_no << "\n";
			exit(1);
		}
	}
}

/*****************************************************************************************************************************/
bool CParser::nextLine(std::string& line)
{
	while (true)
	{
		if (file.eof())
			return false;
		std::getline(file, line);
		++line_no;
		std::smatch match;
		if (!std::regex_search(line, match, empty_line_pattern))
			return true;
	}
}

// ***** EOF