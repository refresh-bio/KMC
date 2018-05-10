/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _EXPRESSION_NODE_H
#define _EXPRESSION_NODE_H
#include "defs.h"
#include "operations.h"
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include "kmc1_db_reader.h"
#include "kmc2_db_reader.h"

//************************************************************************************************************
// CExpressionNode - Base abstract class representing expression node. In first stage of algorithm from
// user input there is created binary tree. Node type represents operation. This tree is only for generating
// another tree (check out CInput and CBundle)
//************************************************************************************************************
template<unsigned SIZE> class CExpressionNode
{
public:
	CExpressionNode() :left(nullptr), right(nullptr)
	{

	}
	CExpressionNode* GetLeftChild() const
	{
		return left;
	}
	CExpressionNode* GetRightChild() const
	{
		return right;
	}

	virtual CBundle<SIZE>* GetExecutionRoot() = 0;

	void AddLeftChild(CExpressionNode* child)
	{
#ifdef ENABLE_DEBUG
		if (left)
		{
			std::cout << "This child node already exists\n";
			exit(1);
		}
#endif
		left = child;
	}

	void AddRightChild(CExpressionNode* child)
	{
#ifdef ENABLE_DEBUG
		if (right)
		{
			std::cout << "This child node already exists\n";
			exit(1);
		}
#endif
		right = child;
	}

#ifdef ENABLE_DEBUG	
	virtual void Info() = 0;
	void Display(int adient = 0)
	{
		if (right)
			right->Display(adient + 5);

		for (int i = 0; i < adient; ++i)
			std::cout << " ";
		Info();
		std::cout << "\n";
		if (left)
			left->Display(adient + 5);
	}
#endif

	virtual ~CExpressionNode()
	{
		delete left;
		delete right;
	}

protected:
	CExpressionNode* left, *right;
};

template<unsigned SIZE> class COperNode : public CExpressionNode<SIZE> 
{
protected:
	CounterOpType counter_op_type;
public:
	COperNode(CounterOpType counter_op_type) :counter_op_type(counter_op_type)
	{

	}
	void SetCounterOpType(CounterOpType _counter_op_type)
	{
		this->counter_op_type = _counter_op_type;
	}
};

//************************************************************************************************************
// CUnionNode - represents node for union operation
//************************************************************************************************************
template<unsigned SIZE> class CUnionNode : public COperNode<SIZE>
{
public:
	CUnionNode():COperNode<SIZE>(CounterOpType::SUM)
	{
	}
	CBundle<SIZE>* GetExecutionRoot() override
	{
		return new CBundle<SIZE>(new CUnion<SIZE>(this->left->GetExecutionRoot(), this->right->GetExecutionRoot(), this->counter_op_type));
	}
#ifdef ENABLE_DEBUG
	void Info() override
	{
		std::cout << "+";
	}
#endif
};

//************************************************************************************************************
// CKmersSubtractionNode - represents node for subtraction of k-mers (if k-mer exists in both input,
//	it is absent in result) operation
//************************************************************************************************************
template<unsigned SIZE> class CKmersSubtractionNode : public COperNode<SIZE>
{
public:
	CKmersSubtractionNode() :COperNode<SIZE>(CounterOpType::NONE)
	{
	}
	CBundle<SIZE>* GetExecutionRoot() override
	{
		return new CBundle<SIZE>(new CKmersSubtract<SIZE>(this->left->GetExecutionRoot(), this->right->GetExecutionRoot()));
	}
#ifdef ENABLE_DEBUG
	void Info() override
	{
		std::cout << "-";
	}
#endif
};


template<unsigned SIZE> class CCountersSubtractionNode : public COperNode<SIZE>
{
public:
	CCountersSubtractionNode():
		COperNode<SIZE>(CounterOpType::DIFF)
	{
	}
	CBundle<SIZE>* GetExecutionRoot() override
	{
		return new CBundle<SIZE>(new CCountersSubtract<SIZE>(this->left->GetExecutionRoot(), this->right->GetExecutionRoot(), this->counter_op_type));
	}
#ifdef ENABLE_DEBUG
	void Info() override
	{
		std::cout << "~";
	}
#endif
};
//************************************************************************************************************
// CIntersectionNode - represents node for intersection operation
//************************************************************************************************************
template<unsigned SIZE> class CIntersectionNode : public COperNode<SIZE>
{
public:	
	CIntersectionNode():
		COperNode<SIZE>(CounterOpType::MIN)
	{
	}
	CBundle<SIZE>* GetExecutionRoot() override
	{		
		return new CBundle<SIZE>(new CIntersection<SIZE>(this->left->GetExecutionRoot(), this->right->GetExecutionRoot(), this->counter_op_type));
	}
#ifdef ENABLE_DEBUG
	void Info() override
	{
		std::cout << "*";
	}
#endif
};

//************************************************************************************************************
// CInputNode - represents node (leaf) - KMC1 or KMC2 database
//************************************************************************************************************
template<unsigned SIZE> class CInputNode : public CExpressionNode<SIZE>
{	
	uint32 desc_pos;
public:
	CInputNode(uint32 desc_pos) : desc_pos(desc_pos)
	{
	}
	CBundle<SIZE>* GetExecutionRoot() override
	{		
		CConfig& config = CConfig::GetInstance();
		CInput<SIZE>* db = nullptr;
		if (!config.headers[desc_pos].IsKMC2())
			db = new CKMC1DbReader<SIZE>(config.headers[desc_pos], config.input_desc[desc_pos], CConfig::GetInstance().percent_progress, KMCDBOpenMode::sorted);
		else		
			db = new CKMC2DbReader<SIZE>(config.headers[desc_pos], config.input_desc[desc_pos], CConfig::GetInstance().percent_progress, KMCDBOpenMode::sorted);
		return new CBundle<SIZE>(db);
	}

#ifdef ENABLE_DEBUG
	void Info() override
	{
		std::cout << "In: " << CConfig::GetInstance().input_desc[desc_pos].file_src;
	}
#endif
};


#endif


// ***** EOF