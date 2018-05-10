/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2017-01-28
*/

#include "stdafx.h"
#include "tokenizer.h"

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

CTokenizer::CTokenizer()
{
	token_patterns.resize(13);
	token_patterns[0] = std::make_pair("^(\\()", TokenType::PARENTHESIS_OPEN);
	token_patterns[1] = std::make_pair("^(\\))", TokenType::PARENTHESIS_CLOSE);
	token_patterns[2] = std::make_pair("^(\\-)", TokenType::STRICT_MINUS_OPER);
	token_patterns[3] = std::make_pair("^(\\~)", TokenType::COUNTER_MINUS_OPER);
	token_patterns[4] = std::make_pair("^(\\+)", TokenType::PLUS_OPER);
	token_patterns[5] = std::make_pair("^(\\*)", TokenType::MUL_OPER);
	
	//those are keywords
	token_patterns[6] = std::make_pair("^(min)", TokenType::MIN_MODIFIER);
	token_patterns[7] = std::make_pair("^(max)", TokenType::MAX_MODIFIER);
	token_patterns[8] = std::make_pair("^(diff)", TokenType::DIFF_MODIFIER);
	token_patterns[9] = std::make_pair("^(sum)", TokenType::SUM_MODIFIER);
	token_patterns[10] = std::make_pair("^(left)", TokenType::LEFT_MODIFIER);
	token_patterns[11] = std::make_pair("^(right)", TokenType::RIGHT_MODIFIER);

	token_patterns[12] = std::make_pair("^(\\w*)", TokenType::VARIABLE);
}


/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/


const std::set<std::string>& CTokenizer::GetKeywords()
{
	static std::set<std::string> keywords = {"min", "max", "sum", "diff", "left", "right"}; //related to tokens created in CTokenizer ctor
	return keywords;
}

void CTokenizer::Tokenize(const std::string& _expression, std::list<Token>& tokens)
{
	std::string expression = _expression;
	std::smatch match;
	leftTrimString(expression, 0);
	while (!expression.empty())
	{
		bool valid_token = false;
		for (const auto& pattern : token_patterns)
		{
			if (std::regex_search(expression, match, pattern.first))
			{
#ifdef ENABLE_DEBUG
				std::cout << match[1];
#endif
				tokens.push_back(std::make_pair(match[1], pattern.second));
				leftTrimString(expression, (int)match[1].length());
				valid_token = true;
				break;
			}
		}
		if (!valid_token)
		{
			std::cerr << "Error: wrong output format near : " << expression << "\n";
			exit(1);
		}
	}
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

void CTokenizer::leftTrimString(std::string& str, int start_pos)
{
	static const std::string whitespace = " \t\r\n\v\f";
	auto next_pos = str.find_first_not_of(whitespace, start_pos);
	str.erase(0, next_pos);
}

// ***** EOF