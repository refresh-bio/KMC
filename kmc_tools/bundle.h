/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _BUNDLE_H
#define _BUNDLE_H
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



//************************************************************************************************************
// CBundleData - class containing a buffer of k-mers and its counters. 
//************************************************************************************************************
template<unsigned SIZE> class CBundleData
{
public:
	CBundleData() : insert_pos(0), get_pos(0), size(BUNDLE_CAPACITY)
	{
		kmers = new CKmer<SIZE>[size]; 
		counters = new uint32[size];
	}
	~CBundleData()
	{
		delete[] kmers;
		delete[] counters;
	}
	CBundleData(CBundleData<SIZE>&& rhs):
		insert_pos(rhs.insert_pos), get_pos(rhs.get_pos), size(rhs.size), kmers(rhs.kmers), counters(rhs.counters)
	{
		rhs.counters = nullptr;
		rhs.kmers = nullptr;
		rhs.get_pos = rhs.size = rhs.insert_pos = 0;
	}

	CBundleData<SIZE>& operator=(CBundleData<SIZE>&& rhs)
	{
		if (this != &rhs)
		{
			delete[] kmers;
			delete[] counters;

			kmers = rhs.kmers;
			counters = rhs.counters;
			get_pos = rhs.get_pos;
			size = rhs.size;
			insert_pos = rhs.insert_pos;

			rhs.counters = nullptr;
			rhs.kmers = nullptr;
			rhs.get_pos = rhs.size = rhs.insert_pos = 0;
		}
		return *this;
	}

	CBundleData(const CBundleData<SIZE>&) = delete;
	CBundle<SIZE>& operator=(const CBundleData<SIZE>&) = delete;

	CKmer<SIZE>& TopKmer() const
	{
		return kmers[get_pos];
	}

	uint32& TopCounter() const
	{
		return counters[get_pos];
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
		kmers[insert_pos] = kmer;
		counters[insert_pos++] = counter;
	}
	void Pop()
	{
		++get_pos;
	}

	void Clear()
	{
		insert_pos = get_pos = 0;
	}

private:
	friend class CBundle<SIZE>;	
	uint32 insert_pos, get_pos, size;
	CKmer<SIZE>* kmers;
	uint32* counters;
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

private:
	CBundleData<SIZE> data;
	CInput<SIZE>* input;
	bool finished = false;
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