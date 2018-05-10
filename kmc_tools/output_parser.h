/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _OUTPUT_PARSER_H
#define _OUTPUT_PARSER_H

#include "defs.h"
#include <list>
#include <map>
#include "tokenizer.h"
#include "expression_node.h"


/*****************************************************************************************************************************/
// This parser validate below grammar:
// expr -> term sum_op
// sum_op -> PLUSMINUS term sum_op
// sum_op -> TERMINATOR
// 
// term -> argument term_op
// term_op -> MUL argument term_op
// term_op -> TERMINATOR
// argument -> VARIABLE
// argument -> OPEN_BRACKET expr CLOSE_BRACKET
// This code is based on: https://github.com/mikailsheikh/cogitolearning-examples/tree/master/CogPar
/*****************************************************************************************************************************/

template<unsigned SIZE> class COutputParser
{
	std::list<Token> tokens;
	const std::map<std::string, uint32>& input;
	Token curr_token;
	void nextToken();
	CExpressionNode<SIZE>* argument();
	CExpressionNode<SIZE>* term_op(CExpressionNode<SIZE>* left);
	CExpressionNode<SIZE>* term();
	CExpressionNode<SIZE>* sum_op(CExpressionNode<SIZE>* left);
	CExpressionNode<SIZE>* expr();

	void modifier(COperNode<SIZE>* exp);
public:
	COutputParser(std::list<Token>& tokens, const std::map<std::string, uint32>& input) :
		tokens(tokens), input(input)
	{
		curr_token = tokens.front();
	}

	CExpressionNode<SIZE>* Parse();
};


/*****************************************************************************************************************************/
/********************************************************** PUBLIC ***********************************************************/
/*****************************************************************************************************************************/

template<unsigned SIZE>
CExpressionNode<SIZE>* COutputParser<SIZE>::Parse()
{
	CExpressionNode<SIZE>* res = expr();
	if (curr_token.second != TokenType::TERMINATOR)
	{
		std::cerr << "Error: wrong symbol :" << curr_token.first <<"\n";
		exit(1);
	}
#ifdef ENABLE_DEBUG
	std::cout << "\n";
	res->Display();
#endif
	return res;
}

/*****************************************************************************************************************************/
/********************************************************** PRIVATE **********************************************************/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
template<unsigned SIZE> void COutputParser<SIZE>::nextToken()
{
	tokens.pop_front();
	if (tokens.empty())
		curr_token.second = TokenType::TERMINATOR;
	else
		curr_token = tokens.front();
}

/*****************************************************************************************************************************/
template<unsigned SIZE> CExpressionNode<SIZE>* COutputParser<SIZE>::argument()
{
	if (curr_token.second == TokenType::VARIABLE)
	{
		//check if this variable was defined
		auto elem = input.find(curr_token.first);
		if (elem == input.end())
		{
			std::cerr << "Error: variable " << curr_token.first << " was not defined\n";
			exit(1);
		}
		CExpressionNode<SIZE>* res = new CInputNode<SIZE>(elem->second);
		nextToken();
		return res;
	}
	else if (curr_token.second == TokenType::PARENTHESIS_OPEN)
	{
		nextToken();
		CExpressionNode<SIZE>* res = expr();
		if (curr_token.second != TokenType::PARENTHESIS_CLOSE)
		{
			std::cerr << "Error: close  parenthesis expected, but " << curr_token.first << " found\n";
			exit(1);
		}
		nextToken();
		return res;
	}
	return nullptr;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> CExpressionNode<SIZE>* COutputParser<SIZE>::term_op(CExpressionNode<SIZE>* left)
{
	if (curr_token.second == TokenType::MUL_OPER)
	{
		COperNode<SIZE>* res = new CIntersectionNode<SIZE>;
		res->AddLeftChild(left);
		nextToken();
		modifier(res);
		auto right = argument();
		res->AddRightChild(right);
		return term_op(res);
	}
	return left;
}
template<unsigned SIZE> CExpressionNode<SIZE>* COutputParser<SIZE>::term()
{
	auto left = argument();
	return term_op(left);
}

/*****************************************************************************************************************************/
template<unsigned SIZE> void COutputParser<SIZE>::modifier(COperNode<SIZE>* exp)
{	
	if (curr_token.second == TokenType::DIFF_MODIFIER)
	{
		exp->SetCounterOpType(CounterOpType::DIFF);
		nextToken();
	}
	else if (curr_token.second == TokenType::LEFT_MODIFIER)
	{
		exp->SetCounterOpType(CounterOpType::FROM_DB1);
		nextToken();
	}
	else if (curr_token.second == TokenType::MAX_MODIFIER)
	{
		exp->SetCounterOpType(CounterOpType::MAX);
		nextToken();
	}
	else if (curr_token.second == TokenType::MIN_MODIFIER)
	{
		exp->SetCounterOpType(CounterOpType::MIN);
		nextToken();
	}
	else if (curr_token.second == TokenType::RIGHT_MODIFIER)
	{
		exp->SetCounterOpType(CounterOpType::FROM_DB2);
		nextToken();
	}
	else if (curr_token.second == TokenType::SUM_MODIFIER)
	{
		exp->SetCounterOpType(CounterOpType::SUM);
		nextToken();
	}
}
/*****************************************************************************************************************************/
template<unsigned SIZE> CExpressionNode<SIZE>* COutputParser<SIZE>::sum_op(CExpressionNode<SIZE>* left)
{
	if (curr_token.second == TokenType::PLUS_OPER || curr_token.second == TokenType::STRICT_MINUS_OPER || curr_token.second == TokenType::COUNTER_MINUS_OPER)
	{
		COperNode<SIZE>* res = nullptr;
		if (curr_token.second == TokenType::PLUS_OPER)
			res = new CUnionNode<SIZE>;
		else if (curr_token.second == TokenType::STRICT_MINUS_OPER)
			res = new CKmersSubtractionNode<SIZE>;
		else
			res = new CCountersSubtractionNode<SIZE>;
		res->AddLeftChild(left);		
		bool modifier_allowed = !(curr_token.second == TokenType::STRICT_MINUS_OPER);
		nextToken();		
		if(modifier_allowed)
			modifier(res);
		auto right = term();
		res->AddRightChild(right);
		return sum_op(res);
	}
	return left;
}

/*****************************************************************************************************************************/
template<unsigned SIZE> CExpressionNode<SIZE>* COutputParser<SIZE>::expr()
{
	auto left = term();
	return sum_op(left);
}



#endif

// ***** EOF