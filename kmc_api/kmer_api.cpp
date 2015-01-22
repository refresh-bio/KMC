/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz and Agnieszka Debudaj-Grabysz

  Version: 2.1.1
  Date   : 2015-01-22
*/


#include "stdafx.h"
#include "kmer_api.h"
#include <vector>
#include <math.h>

using namespace std;

const char CKmerAPI::char_codes[] = {'A','C', 'G', 'T'};	
char CKmerAPI::num_codes[];
CKmerAPI::_si CKmerAPI::_init; 

// ***** EOF
