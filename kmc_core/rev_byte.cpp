/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.2
  Date   : 2023-03-09
*/

#include "rev_byte.h"

uchar CRev_byte::lut[256];
CRev_byte::_si CRev_byte::_init;

// ***** EOF