/*
    This file is a part of KMC software distributed under GNU GPL 3 licence.
    The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

    Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

    Version: 2.2.0
    Date   : 2015-04-15
*/

#ifndef _ASMLIB_WRAPPER_H
#define _ASMLIB_WRAPPER_H

#include "../kmc/definitions.h"

#ifdef DISABLE_ASMLIB
#define A_memcpy memcpy
#define SetMemcpyCacheLimit(X)
#else
#include "../external/asmlib.h"
#endif

#endif
