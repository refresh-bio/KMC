/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _META_OPER_H
#define _META_OPER_H

//#include <functional>


template <size_t N> struct uint_{ };

// For loop (forward)
template <size_t N, typename Lambda>
inline void IterFwd(const Lambda &oper, uint_<N>) {
	IterFwd(oper, uint_<N-1>());
	oper(N);
}

template <typename Lambda>
inline void IterFwd(const Lambda &oper, uint_<0>) {
	oper(0);
}

// For loop (backward)
template <size_t N, typename Lambda>
inline void IterRev(const Lambda &oper, uint_<N>) {
	oper(N);
	IterRev(oper, uint_<N-1>());
}

template <typename Lambda>
inline void IterRev(const Lambda &oper, uint_<0>) {
	oper(0);
}

#endif

// ***** EOF
