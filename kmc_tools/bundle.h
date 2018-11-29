/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#ifndef _BUNDLE_H
#define _BUNDLE_H
#include "config.h"
#include "defs.h"
#include "kmer.h"

//************************************************************************************************************
// CBundle and CInput are CORE classes of this application. CInputs are nodes of binary tree which 
// represent operations. Leafs of this tree are kmc database (1 or 2) inputs (sets of k-mers).
// Each node represents an operation like intersection, subtraction, etc. Because this class is abstract 
// calling virtual method to get each single k-mer may be costly. To prevent high const, between tree nodes there
//are instances of CBundle which contains buffer of k-mers and its counters.
//
// The algorithm works as follow (conceptually):
// Build a tree with CBundles and CInputs, as a root take a some output writer (kmc database). 
// Root has its bundle and get from it k-mers, but at the beginning there is nothing in bundle. Each bundle 
// contains pointer to CInput below in tree. The CBundle is getting k-mers from its CInput. 
// This is repeated from top of tree to leafs
//************************************************************************************************************

//Forward declaration
template<unsigned SIZE> class CBundle;

//************************************************************************************************************
// CInput - Base abstract class representing data source for CBundle class
//************************************************************************************************************
template<unsigned SIZE> class CInput
{
public:
	virtual void NextBundle(CBundle<SIZE>& bundle) = 0;
	virtual void IgnoreRest() = 0;
	bool Finished(){ return finished; }
	virtual ~CInput(){}
protected:
	bool finished = false;

};

template<unsigned SIZE> class CSimpleOperation; 
template<unsigned SIZE> class CKMC1DbReader; 
template<unsigned SIZE> class CMergerParent;
template<unsigned SIZE> class CMergerParentSubthread;

//************************************************************************************************************
// CBundleData - class containing a buffer of k-mers and its counters. 
//************************************************************************************************************
template<unsigned SIZE> class CBundleData
{
	struct CKmerWithCounter
	{
		CKmer<SIZE> kmer;
		uint32 counter;
	};

public:
	CBundleData(uint32 size) : insert_pos(0), get_pos(0), size(size)
	{
		kmers_with_counters = new CKmerWithCounter[size];
	}
	CBundleData() : insert_pos(0), get_pos(0), size(BUNDLE_CAPACITY)
	{
		kmers_with_counters = new CKmerWithCounter[size];
	}
	~CBundleData()
	{
		delete[] kmers_with_counters;
	}
	CBundleData(CBundleData<SIZE>&& rhs) :
		insert_pos(rhs.insert_pos), get_pos(rhs.get_pos), size(rhs.size), kmers_with_counters(rhs.kmers_with_counters)
	{
		rhs.kmers_with_counters = nullptr;
		rhs.get_pos = rhs.size = rhs.insert_pos = 0;
	}

	CBundleData<SIZE>& operator=(CBundleData<SIZE>&& rhs)
	{
		if (this != &rhs)
		{
			delete[] kmers_with_counters;

			kmers_with_counters = rhs.kmers_with_counters;
			get_pos = rhs.get_pos;
			size = rhs.size;
			insert_pos = rhs.insert_pos;

			rhs.kmers_with_counters = nullptr;
			rhs.get_pos = rhs.size = rhs.insert_pos = 0;
		}
		return *this;
	}

	//deprecated
	//void CopyFrom(CBundleData<SIZE>& rhs) //similar to assign operator but I want assign operator deleted, this method should be used carefully. this->size and rhs.size must quals
	//{
	//	memcpy(kmers_with_counters, rhs.kmers_with_counters, rhs.insert_pos * sizeof(CKmerWithCounter));
	//	insert_pos = rhs.insert_pos;
	//	get_pos = 0;
	//}


	CBundleData(const CBundleData<SIZE>&) = delete;
	CBundle<SIZE>& operator=(const CBundleData<SIZE>&) = delete;

	CKmer<SIZE>& TopKmer() const
	{
		return kmers_with_counters[get_pos].kmer;
	}

	uint32& TopCounter() const
	{
		return kmers_with_counters[get_pos].counter;
	}

	bool Full()
	{
		return insert_pos >= size;
	}

	bool Empty()
	{
		return get_pos >= insert_pos;
	}

	void Insert(CKmer<SIZE>& kmer, uint32 counter)
	{
		kmers_with_counters[insert_pos].kmer = kmer;
		kmers_with_counters[insert_pos++].counter = counter;
	}
	void Pop()
	{
		++get_pos;
	}

	void Clear()
	{
		insert_pos = get_pos = 0;
	}

	uint32 NRecLeft()
	{
		return insert_pos - get_pos;
	}


private:
	friend class CKMC1DbReader<SIZE>; //improve performance, but CKMC1DbReader takes responsibility for CBundleData state!
	friend class CSimpleOperation<SIZE>; //as above
	friend class CMergerParent<SIZE>;
	friend class CMergerParentSubthread<SIZE>;
	friend class CBundle<SIZE>;
	uint32 insert_pos, get_pos, size;
	CKmerWithCounter* kmers_with_counters;
};


//************************************************************************************************************
// CBundle - connector between CBundleData and CInput
//************************************************************************************************************
template<unsigned SIZE> class CBundle
{
public:
	CBundle(CInput<SIZE>* input) : input(input)
	{
		
	}

	//deprecated
	//void CopyFrom(CBundle<SIZE>& rhs) //Use carefully. Look at comment in called function
	//{
	//	data.CopyFrom(rhs.data);
	//}

	CKmer<SIZE>& TopKmer() const
	{
		return data.TopKmer();
	}

	uint32& TopCounter() const
	{
		return data.TopCounter();
	}


	bool Full()
	{		
		return data.Full();
	}
	void Insert(CKmer<SIZE>& kmer, uint32 counter)
	{		
		data.Insert(kmer, counter);
	}
	void Pop()
	{		
		data.Pop();
	}
	~CBundle()
	{
		delete input;
	}

	bool Empty()
	{
		return data.Empty();
	}
	
	CBundleData<SIZE>& Data() {
		return data;
	}

	inline bool Finished();

	void IgnoreRest()
	{		
		input->IgnoreRest();
	}
	uint32 Size()
	{
		return data.insert_pos;
	}

	uint32 NRecLeft()
	{
		return data.insert_pos - data.get_pos;
	}
	
protected:
	friend class CSimpleOperation<SIZE>; //improve performance
	CBundleData<SIZE> data;
	CInput<SIZE>* input;
	bool finished = false;
};

//forward declaration
template <unsigned SIZE> class CKMC1DbWriter;
template<unsigned SIZE>
class COutputBundle : public CBundle<SIZE>
{
private:
	CSimpleOutputDesc::OpType op_type;
	CounterOpType counter_op;
	CKMC1DbWriter<SIZE>& db_writer;
	
public:

	uint32 GetCounter(uint32 counter1, uint32 counter2)
	{
		switch (counter_op)
		{
		case CounterOpType::DIFF:
			return counter1 > counter2 ? counter1 - counter2 : 0;
		case CounterOpType::MAX:
			return MAX(counter1, counter2);
		case CounterOpType::MIN:
			return MIN(counter1, counter2);
		case CounterOpType::SUM:
			return counter1 + counter2;
		case CounterOpType::FROM_DB1:
			return counter1;
		case CounterOpType::FROM_DB2:
			return counter2;
		case CounterOpType::NONE://should never be here
			std::cerr << "Error: trying to use undefined counter calculation mode!\n";
			exit(1);
		}
		return 0;
	}
	
	CSimpleOutputDesc::OpType GetOpType()
	{
		return op_type;
	}


	COutputBundle(CSimpleOutputDesc::OpType op_type, CounterOpType counter_op, CKMC1DbWriter<SIZE>& db_writer) :
		CBundle<SIZE>(nullptr),
		op_type(op_type),
		counter_op(counter_op),
		db_writer(db_writer)
	{

	}

	void InsertAndSendIfFull(CKmer<SIZE>& kmer, uint32 counter)
	{
		this->data.Insert(kmer, counter);
		if (this->Full())
		{			
			db_writer.MultiOptputAddResultPart(*this);
		}
	}

	void NotifyFinish()
	{
		if (!this->Empty())
			db_writer.MultiOptputAddResultPart(*this);
	}
};


//************************************************************************************************************
template<unsigned SIZE> inline bool CBundle<SIZE>::Finished()
{	
	if (finished)
		return true;
	if (data.get_pos >= data.insert_pos)
	{
		if (input->Finished())
		{
			finished = true;
			return true;
		}
		data.get_pos = data.insert_pos = 0;
		input->NextBundle(*this); 
		if (data.insert_pos == 0)//Because maybe NextBundle did not add anything, which means there is nothing to take
		{
			finished = true;
			return true;
		}
	}
	return false;
}

#endif


// ***** EOF