/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.1.1
Date   : 2019-05-19
*/

#ifndef _DB_WRITER_H
#define _DB_WRITER_H
#include "bundle.h"

template<unsigned SIZE>
class CDbWriter
{
public:
	virtual bool Process() = 0;
	virtual void MultiOptputInit() = 0;
	virtual void MultiOptputAddResultPart(COutputBundle<SIZE>& bundle) = 0;
	virtual void MultiOptputAddResultPart(CBundle<SIZE>& bundle) = 0;
	virtual void MultiOptputFinish() = 0;
	virtual ~CDbWriter() = default;
};
#endif
