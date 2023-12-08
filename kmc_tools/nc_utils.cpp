/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.2.3
  Date   : 2023-12-08
*/

#include "nc_utils.h"


uchar CNumericConversions::digits[100000*5];
int CNumericConversions::powOf10[30];
CNumericConversions::_si CNumericConversions::_init;