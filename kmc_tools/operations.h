/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#ifndef _OPERATIONS_H
#define _OPERATIONS_H


#ifdef ENABLE_DEBUG
#include "config.h"
#endif
#include <iostream>
#include "bundle.h"

//************************************************************************************************************
// C2ArgOper - abstract class representing 2 argument's operation
//************************************************************************************************************
template<unsigned SIZE> class C2ArgOper : public CInput<SIZE>
{
protected:
	CBundle<SIZE>* input1, *input2;
public:
	C2ArgOper(CBundle<SIZE>* input1, CBundle<SIZE>* input2) :
		input1(input1), input2(input2)
	{
	}

	void IgnoreRest() override
	{		
		input1->IgnoreRest();
		input2->IgnoreRest();
	}

	~C2ArgOper() override
	{
		delete input1;
		delete input2;
	}
};

//************************************************************************************************************
// CUnion - implementation of union operation on 2 k-mer's sets.
//************************************************************************************************************
template <unsigned SIZE> class CUnion : public C2ArgOper<SIZE>
{		
public:
	CUnion(CBundle<SIZE>* input1, CBundle<SIZE>* input2) : C2ArgOper<SIZE>(input1, input2)
	{
	}
	void NextBundle(CBundle<SIZE>& bundle) override
	{
		while (!this->input1->Finished() && !this->input2->Finished())
		{
			if (bundle.Full())
			{
				return;
			}	
			if (this->input1->TopKmer() == this->input2->TopKmer())
			{
				bundle.Insert(this->input1->TopKmer(), this->input1->TopCounter() + this->input2->TopCounter());
				this->input1->Pop();
				this->input2->Pop();
			}
			else if (this->input1->TopKmer() < this->input2->TopKmer())
			{
				bundle.Insert(this->input1->TopKmer(), this->input1->TopCounter());
				this->input1->Pop();
			}
			else
			{
				bundle.Insert(this->input2->TopKmer(), this->input2->TopCounter());
				this->input2->Pop();
			}
		}
		CBundle<SIZE>* non_empty_bundle = this->input1->Finished() ? this->input2 : this->input1;
		while (!non_empty_bundle->Finished())
		{
			if (bundle.Full())
			{
				return;
			}
			bundle.Insert(non_empty_bundle->TopKmer(), non_empty_bundle->TopCounter());
			non_empty_bundle->Pop();
		}
		this->finished = true;
	}
};

//************************************************************************************************************
// CIntersection - implementation of intersection operation on 2 k-mer's sets.
//************************************************************************************************************
template<unsigned SIZE> class CIntersection : public C2ArgOper<SIZE>
{		
public:
	CIntersection(CBundle<SIZE>* input1, CBundle<SIZE>* input2) : C2ArgOper<SIZE>(input1, input2)
	{
	}
	void NextBundle(CBundle<SIZE>& bundle) override
	{
		while (!this->input1->Finished() && !this->input2->Finished())
		{
			/*this->input1->Top(kmer1, counter1);
			this->input2->Top(kmer2, counter2);*/

			if (this->input1->TopKmer() == this->input2->TopKmer())
			{
				bundle.Insert(this->input1->TopKmer(), MIN(this->input1->TopCounter(), this->input2->TopCounter()));
				this->input1->Pop();
				this->input2->Pop();
				if (bundle.Full())
					return;
			}
			else if (this->input1->TopKmer() < this->input2->TopKmer())
				this->input1->Pop();
			else
				this->input2->Pop();
		}
		if (!this->input1->Finished())
			this->input1->IgnoreRest();
		if (!this->input2->Finished())
			this->input2->IgnoreRest();
		this->finished = true;
	}
};

//************************************************************************************************************
// CKmersSubtract - implementation of subtraction operation of 2 k-mer's sets.
// If k-mer exists in both input it is absent in result (counters does not matter).
//************************************************************************************************************
template<unsigned SIZE> class CKmersSubtract : public C2ArgOper<SIZE>
{	
	//CKmer<SIZE> kmer1, kmer2;
	//uint32 counter1, counter2;
public:
	CKmersSubtract(CBundle<SIZE>* input1, CBundle<SIZE>* input2) : C2ArgOper<SIZE>(input1, input2)
	{
	}
	void NextBundle(CBundle<SIZE>& bundle) override
	{
		while (!this->input1->Finished() && !this->input2->Finished())
		{
			//this->input1->Top(kmer1, counter1);
			//this->input2->Top(kmer2, counter2);
			if (this->input2->TopKmer() < this->input1->TopKmer())
				this->input2->Pop();
			else if (this->input2->TopKmer() == this->input1->TopKmer())
			{
				this->input1->Pop();
				this->input2->Pop();
			}
			else
			{
				bundle.Insert(this->input1->TopKmer(), this->input1->TopCounter());
				this->input1->Pop();
				if (bundle.Full())
					return;
			}
		}
		
		if(!this->input2->Finished())
			this->input2->IgnoreRest(); 

		while (!this->input1->Finished())
		{
			if (bundle.Full())
				return;
			bundle.Insert(this->input1->TopKmer(), this->input1->TopCounter());
			this->input1->Pop();
		}
		this->finished = true;
	}
};




//************************************************************************************************************
// CCountersSubtract - implementation of subtraction operation of 2 k-mer's sets.
// If k-mer exists in both input their counters are subtracted.
//************************************************************************************************************
template<unsigned SIZE> class CCountersSubtract : public C2ArgOper<SIZE>
{
	//CKmer<SIZE> kmer1, kmer2;
	//uint32 counter1, counter2;
public:
	CCountersSubtract(CBundle<SIZE>* input1, CBundle<SIZE>* input2) : C2ArgOper<SIZE>(input1, input2)
	{
	}
	void NextBundle(CBundle<SIZE>& bundle) override
	{
		while (!this->input1->Finished() && !this->input2->Finished())
		{
			//this->input1->Top(kmer1, counter1);
			//this->input2->Top(kmer2, counter2);
			if (this->input2->TopKmer() < this->input1->TopKmer())
				this->input2->Pop();
			else if (this->input2->TopKmer() == this->input1->TopKmer())
			{
				if (this->input1->TopCounter() > this->input2->TopCounter())
				{
					bundle.Insert(this->input1->TopKmer(), this->input1->TopCounter() - this->input2->TopCounter());
					this->input1->Pop();
					this->input2->Pop();
					if (bundle.Full())
						return;
				}
				else
				{
					this->input1->Pop();
					this->input2->Pop();
				}				
			}
			else
			{
				bundle.Insert(this->input1->TopKmer(), this->input1->TopCounter());
				this->input1->Pop();
				if (bundle.Full())
					return;
			}
		}

		if (!this->input2->Finished())
			this->input2->IgnoreRest();

		while (!this->input1->Finished())
		{
			if (bundle.Full())
				return;			
			bundle.Insert(this->input1->TopKmer(), this->input1->TopCounter());
			this->input1->Pop();
		}
		this->finished = true;
	}
};

template<unsigned SIZE> class CComparer 
{
	CBundle<SIZE>* input1, *input2;
public:
	CComparer(CBundle<SIZE>* input1, CBundle<SIZE>* input2) : input1(input1), input2(input2)
	{
	}

	bool Equals()
	{		
		while (!this->input1->Finished() && !this->input2->Finished())
		{
			if (this->input1->TopCounter() != this->input2->TopCounter())
			{
				this->input1->IgnoreRest();
				this->input2->IgnoreRest();
				return false;
			}
			if (!(this->input1->TopKmer() == this->input2->TopKmer()))
			{
				this->input1->IgnoreRest();
				this->input2->IgnoreRest();
				return false;
			}

			this->input1->Pop();
			this->input2->Pop();
		}
		if (!this->input1->Finished() || !this->input2->Finished())
		{
			std::cout << "one of input is not finished\n";
			this->input1->IgnoreRest();
			this->input2->IgnoreRest();
			return false;
		}
		return true;
	}
};


#endif

// ***** EOF