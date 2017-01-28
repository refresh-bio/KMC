/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.0.0
Date   : 2017-01-28
*/

#ifndef _PROB_QUAL_H
#define _PROB_QUAL_H

struct CProbQual
{
	static double prob_qual[94];
	static double inv_prob_qual[94];
	static double MIN_PROB_QUAL_VALUE;
};
#endif