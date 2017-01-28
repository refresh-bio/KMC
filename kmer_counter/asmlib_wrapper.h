/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#ifndef _ASMLIB_WRAPPER_H
#define _ASMLIB_WRAPPER_H

#include "defs.h"
#ifdef DISABLE_ASMLIB
#define A_memcpy memcpy
#define A_memmove memmove
#define SetMemcpyCacheLimit(X) 
#else
#include "libs/asmlib.h"
#endif

#endif