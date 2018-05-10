/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _TOKENIZER_H
#define _TOKENIZER_H

#include "defs.h"
#include <vector>
#include <regex>
#include <list>
#include <set>
#include <iostream>

enum class TokenType{ VARIABLE, PLUS_OPER, STRICT_MINUS_OPER, COUNTER_MINUS_OPER, MUL_OPER, PARENTHESIS_OPEN, PARENTHESIS_CLOSE, TERMINATOR, DIFF_MODIFIER, SUM_MODIFIER, MIN_MODIFIER, MAX_MODIFIER, LEFT_MODIFIER, RIGHT_MODIFIER };
using Token = std::pair<std::string, TokenType>;

//************************************************************************************************************
// CTokenizer - Tokenizer for k-mers set operations
//************************************************************************************************************
class CTokenizer
{
public:
	static const std::set<std::string>& GetKeywords();
	CTokenizer();
	void Tokenize(const std::string& _expression, std::list<Token>& tokens);

private:
	std::vector<std::pair<std::regex, TokenType>> token_patterns;
	void leftTrimString(std::string& str, int start_pos);	
};

#endif

// ***** EOF