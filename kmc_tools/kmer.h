/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.0
Date   : 2018-05-10
*/

#ifndef _KMER_H
#define _KMER_H

// Important remark: there is no inheritance here to guarantee that all classes defined here are POD according to C++11

#include "defs.h"
#include "meta_oper.h"
#include <string>
#include <cstring>

// *************************************************************************
// Ckmer class for k > 32 with classic kmer counting
template<unsigned SIZE> struct CKmer {
	unsigned long long data[SIZE];

	typedef unsigned long long data_t;
	inline void set(const CKmer<SIZE> &x);

	inline void mask(const CKmer<SIZE> &x);
	inline uint32 end_mask(const uint32 mask);
	inline void set_2bits(const uint64 x, const uint32 p);
	inline uchar get_2bits(const uint32 p);
	inline uchar get_byte(const uint32 p);
	inline void set_byte(const uint32 p, uchar x);
	inline void set_bytes(const uint32 p, const uint32 n, uint32 x);
	inline void set_bits(const uint32 p, const uint32 n, uint64 x);

	inline void increment_at(uint32 suffix_bytes);

	inline void SHL_insert_2bits(const uint64 x);
	inline void SHR_insert_2bits(const uint64 x, const uint32 p);

	inline void SHR(const uint32 p);
	inline void SHL(const uint32 p);

	inline uint64 remove_suffix(const uint32 n) const;
	inline void set_n_1(const uint32 n);
	inline void set_n_01(const uint32 n);

	inline void store(uchar *&buffer, int32 n);
	inline void store(uchar *buffer, int32 p, int32 n);

	inline void load(uchar *&buffer, int32 n);

	inline void load_fast(uchar *&buffer, int32 n, bool little_endian);


	inline bool operator==(const CKmer<SIZE> &x);
	inline bool operator<(const CKmer<SIZE> &x)const;

	inline void clear(void);

	inline char get_symbol(int p);


	inline void set_prefix(CKmer<SIZE>& rhs, uint32 suffix_bytes);
};


// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set(const CKmer<SIZE> &x)
{
#ifdef USE_META_PROG
	IterFwd([&](const int &i){
		data[i] = x.data[i];
	}, uint_<SIZE - 1>());
#else
	for (uint32 i = 0; i < SIZE; ++i)
		data[i] = x.data[i];
#endif
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::mask(const CKmer<SIZE> &x)
{
#ifdef USE_META_PROG
	IterFwd([&](const int &i){
		data[i] &= x.data[i];
	}, uint_<SIZE - 1>());
#else
	for (uint32 i = 0; i < SIZE; ++i)
		data[i] &= x.data[i];
#endif
}

// *********************************************************************
template<unsigned SIZE> inline uint32 CKmer<SIZE>::end_mask(const uint32 mask)
{
	return data[0] & mask;
}
// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_2bits(const uint64 x, const uint32 p)
{
	//	data[p >> 6] |= x << (p & 63);
	data[p >> 6] += x << (p & 63);
}

template<unsigned SIZE> inline uchar CKmer<SIZE>::get_2bits(const uint32 p)
{
	return (data[p >> 6] >> (p & 63)) & 3;
}
// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHR_insert_2bits(const uint64 x, const uint32 p)
{
#ifdef USE_META_PROG
	IterFwd([&](const int &i){
		data[i] >>= 2;
		//		data[i] |= data[i+1] << (64-2);
		data[i] += data[i + 1] << (64 - 2);
	}, uint_<SIZE - 2>());
#else
	for (uint32 i = 0; i < SIZE - 1; ++i)
	{
		data[i] >>= 2;
		//		data[i] |= data[i+1] << (64-2);
		data[i] += data[i + 1] << (64 - 2);
	}
#endif
	data[SIZE - 1] >>= 2;

	//	data[p >> 6] |= x << (p & 63);
	data[p >> 6] += x << (p & 63);
}



// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHR(const uint32 p)
{
#ifdef USE_META_PROG
	IterFwd([&](const int &i){
		data[i] >>= 2 * p;
		//		data[i] |= data[i+1] << (64-2*p);
		data[i] += data[i + 1] << (64 - 2 * p);
	}, uint_<SIZE - 2>());
#else
	for (uint32 i = 0; i < SIZE - 1; ++i)
	{
		data[i] >>= 2 * p;
		//		data[i] |= data[i+1] << (64-2*p);
		data[i] += data[i + 1] << (64 - 2 * p);
	}
#endif
	data[SIZE - 1] >>= 2 * p;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHL(const uint32 p)
{
#ifdef USE_META_PROG
	IterRev([&](const int &i){
		data[i + 1] <<= p * 2;
		//		data[i+1] |= data[i] >> (64-p*2);
		data[i + 1] += data[i] >> (64 - p * 2);
	}, uint_<SIZE - 2>());
#else
	for (uint32 i = SIZE - 1; i > 0; --i)
	{
		data[i] <<= p * 2;
		//		data[i] |= data[i-1] >> (64-p*2);
		data[i] += data[i - 1] >> (64 - p * 2);
	}
#endif
	data[0] <<= p * 2;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHL_insert_2bits(const uint64 x)
{
#ifdef USE_META_PROG
	IterRev([&](const int &i){
		data[i + 1] <<= 2;
		//		data[i+1] |= data[i] >> (64-2);
		data[i + 1] += data[i] >> (64 - 2);
	}, uint_<SIZE - 2>());
#else
	for (uint32 i = SIZE - 1; i > 0; --i)
	{
		data[i] <<= 2;
		//		data[i] |= data[i-1] >> (64-2);
		data[i] += data[i - 1] >> (64 - 2);
	}
#endif
	data[0] <<= 2;
	//	data[0] |= x;
	data[0] += x;
}

// *********************************************************************
template<unsigned SIZE> inline uchar CKmer<SIZE>::get_byte(const uint32 p)
{
	return (data[p >> 3] >> ((p << 3) & 63)) & 0xFF;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_byte(const uint32 p, uchar x)
{
	//	data[p >> 3] |= ((uint64) x) << ((p & 7) << 3);
	data[p >> 3] += ((uint64)x) << ((p & 7) << 3);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_bytes(const uint32 p, const uint32 n, uint32 x)
{
	data[p >> 3] += ((uint64)x) << ((p & 7) << 3);
	//if (8 - (p & 7) < n)
	if((p>>3) != ((p + n - 1) >> 3))
		data[(p >> 3) + 1] += ((uint64)x) >> ((8 - (p & 7)) << 3);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_bits(const uint32 p, const uint32 n, uint64 x)
{
	//	data[p >> 6] |= x << (p & 63);
	data[p >> 6] += x << (p & 63);
	if ((p >> 6) != ((p + n - 1) >> 6))
		//		data[(p >> 6) + 1] |= x >> (64 - (p & 63));
		data[(p >> 6) + 1] += x >> (64 - (p & 63));
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::increment_at(uint32 suffix_bytes)
{
	data[(suffix_bytes) >> 3] += 1ull << ((suffix_bytes & 7) << 3);
	if (((suffix_bytes) >> 3) < (SIZE - 1) && (data[(suffix_bytes) >> 3]) == 0) //overflow
		++data[SIZE - 1];

}

// *********************************************************************
template<unsigned SIZE> inline bool CKmer<SIZE>::operator==(const CKmer<SIZE> &x) {
	for (uint32 i = 0; i < SIZE; ++i)
	if (data[i] != x.data[i])
		return false;

	return true;
}

// *********************************************************************
template<unsigned SIZE> inline bool CKmer<SIZE>::operator<(const CKmer<SIZE> &x)const {
	for (int32 i = SIZE - 1; i >= 0; --i)
	if (data[i] < x.data[i])
		return true;
	else if (data[i] > x.data[i])
		return false;
	return false;
}



// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::clear(void)
{
#ifdef USE_META_PROG
	IterFwd([&](const int &i){
		data[i] = 0;
	}, uint_<SIZE - 1>());
#else	
	for (uint32 i = 0; i < SIZE; ++i)
		data[i] = 0;
#endif
}

// *********************************************************************
template<unsigned SIZE> inline uint64 CKmer<SIZE>::remove_suffix(const uint32 n) const
{
	uint32 p = n >> 6; // / 64;
	uint32 r = n & 63;	// % 64;

	if (p == SIZE - 1)
		return data[p] >> r;
	else
		//		return (data[p+1] << (64-r)) | (data[p] >> r);
		return (data[p + 1] << (64 - r)) + (data[p] >> r);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_n_1(const uint32 n)
{
	clear();

	for (uint32 i = 0; i < (n >> 6); ++i)
		data[i] = ~((uint64)0);

	uint32 r = n & 63;

	if (r)
		data[n >> 6] = (1ull << r) - 1;
}


// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_n_01(const uint32 n)
{
	clear();

	for (uint32 i = 0; i < n; ++i)
	if (!(i & 1))
		//			data[i >> 6] |= (1ull << (i & 63));
		data[i >> 6] += (1ull << (i & 63));
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::store(uchar *&buffer, int32 n)
{
	for (int32 i = n - 1; i >= 0; --i)
		*buffer++ = get_byte(i);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::store(uchar *buffer, int32 p, int32 n)
{
	for (int32 i = n - 1; i >= 0; --i)
		buffer[p++] = get_byte(i);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::load(uchar *&buffer, int32 n)
{
	clear();
	for (int32 i = n - 1; i >= 0; --i)
		set_byte(i, *buffer++);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::load_fast(uchar *&buffer, int32 n, bool little_endian)
{	
	uint32 p = ((n + 7) >> 3) - 1;
	for (uint32 i = SIZE - 1; i > p; --i)
		data[i] = 0;
	
	if (!(n & 7))
		++p;
	else
	{
		memcpy(&data[p], buffer, sizeof(uint64));
		if (little_endian)
			data[p] = _bswap_uint64(data[p]);
		data[p] >>= (sizeof(uint64)-(n & 7)) << 3;
		buffer += n & 7;
	}
	for (int i = p - 1; i >= 0; --i)
	{
		memcpy(&data[i], buffer, sizeof(uint64));
		if (little_endian)
			data[i] = _bswap_uint64(data[i]);
		buffer += 8;
	}
}


// *********************************************************************
template<unsigned SIZE> inline char CKmer<SIZE>::get_symbol(int p)
{
	uint32 x = (data[p >> 5] >> (2 * (p & 31))) & 0x03;

	switch (x)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	default: return 'T';
	}
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_prefix(CKmer<SIZE>& rhs, uint32 suffix_bytes)
{
	data[(suffix_bytes) >> 3] += rhs.data[(suffix_bytes) >> 3];
	if (((suffix_bytes) >> 3) < SIZE - 1)
		data[SIZE - 1] += rhs.data[SIZE - 1];
}
// *********************************************************************
// *********************************************************************
// *********************************************************************
// *********************************************************************
// Ckmer class for k <= 32 with classic kmer counting
template<> struct CKmer<1> {
	unsigned long long data;

	typedef unsigned long long data_t;
	static uint32 QUALITY_SIZE;
	static uint32 KMER_SIZE;

	void set(const CKmer<1> &x);	

	void mask(const CKmer<1> &x);
	uint32 end_mask(const uint32 mask);
	void set_2bits(const uint64 x, const uint32 p);
	uchar get_2bits(const uint32 p);
	uchar get_byte(const uint32 p);
	void set_byte(const uint32 p, uchar x);
	void set_bytes(const uint32 p, const uint32 n, uint32 x);
	void set_bits(const uint32 p, const uint32 n, uint64 x);

	void increment_at(uint32 suffix_bytes);

	void SHL_insert_2bits(const uint64 x);
	void SHR_insert_2bits(const uint64 x, const uint32 p);

	void SHR(const uint32 p);
	void SHL(const uint32 p);

	uint64 remove_suffix(const uint32 n) const;
	void set_n_1(const uint32 n);
	void set_n_01(const uint32 n);

	void store(uchar *&buffer, int32 n);
	void store(uchar *buffer, int32 p, int32 n);


	void load(uchar *&buffer, int32 n);

	void load_fast(uchar *&buffer, int32 n, bool little_endian);

	bool operator==(const CKmer<1> &x);
	bool operator<(const CKmer<1> &x)const;

	void clear(void);

	inline char get_symbol(int p);

	inline void set_prefix(CKmer<1>& rhs, uint32 suffix_bytes);
};


// *********************************************************************
inline void CKmer<1>::mask(const CKmer<1> &x)
{
	data &= x.data;
}


// *********************************************************************
inline uint32 CKmer<1>::end_mask(const uint32 mask)
{
	return data & mask;
}
// *********************************************************************
inline void CKmer<1>::set(const CKmer<1> &x)
{
	data = x.data;
}

// *********************************************************************
inline void CKmer<1>::set_2bits(const uint64 x, const uint32 p)
{
	//	data |= x << p;
	data += x << p;
}

inline uchar CKmer<1>::get_2bits(const uint32 p)
{
	return (data >> p) & 3;
}
// *********************************************************************
inline void CKmer<1>::SHR_insert_2bits(const uint64 x, const uint32 p)
{
	data >>= 2;
	//	data |= x << p;
	data += x << p;
}

// *********************************************************************
inline void CKmer<1>::SHR(const uint32 p)
{
	data >>= 2 * p;
}

// *********************************************************************
inline void CKmer<1>::SHL(const uint32 p)
{
	data <<= p * 2;
}
// *********************************************************************
inline void CKmer<1>::SHL_insert_2bits(const uint64 x)
{
	//	data = (data << 2) | x;
	data = (data << 2) + x;
}

// *********************************************************************
inline uchar CKmer<1>::get_byte(const uint32 p)
{
	return (data >> (p << 3)) & 0xFF;
}

// *********************************************************************
inline void CKmer<1>::set_byte(const uint32 p, uchar x)
{
	//	data |= ((uint64) x) << (p << 3);
	data += ((uint64)x) << (p << 3);
}

// *********************************************************************
inline void CKmer<1>::set_bytes(const uint32 p, const uint32 /*n*/, uint32 x)
{
	data += ((uint64) x) << (p << 3);
}

// *********************************************************************
inline void CKmer<1>::set_bits(const uint32 p, const uint32 /*n*/, uint64 x)
{
	//	data |= x << p;
	data += x << p;
}

// *********************************************************************
inline void CKmer<1>::increment_at(uint32 suffix_bytes)
{
	data += (1ull << (suffix_bytes << 3));
}
// *********************************************************************
inline bool CKmer<1>::operator==(const CKmer<1> &x) {
	return data == x.data;
}

// *********************************************************************
inline bool CKmer<1>::operator<(const CKmer<1> &x)const {
	return data < x.data;
}

// *********************************************************************
inline void CKmer<1>::clear(void)
{
	data = 0ull;
}

// *********************************************************************
inline uint64 CKmer<1>::remove_suffix(const uint32 n) const
{
	return data >> n;
}

// *********************************************************************
inline void CKmer<1>::set_n_1(const uint32 n)
{
	if (n == 64)
		data = ~(0ull);
	else
		data = (1ull << n) - 1;
}

// *********************************************************************
inline void CKmer<1>::set_n_01(const uint32 n)
{
	data = 0ull;

	for (uint32 i = 0; i < n; ++i)
	if (!(i & 1))
		data += (1ull << i);
}

// *********************************************************************
inline void CKmer<1>::store(uchar *&buffer, int32 n)
{
	for (int32 i = n - 1; i >= 0; --i)
		*buffer++ = get_byte(i);
}

// *********************************************************************
inline void CKmer<1>::store(uchar *buffer, int32 p, int32 n)
{
	for (int32 i = n - 1; i >= 0; --i)
		buffer[p++] = get_byte(i);
}

// *********************************************************************
inline void CKmer<1>::load_fast(uchar *&buffer, int32 n, bool little_endian)
{
	//for short k-mers n may be 0
	if (!n)
	{
		data = 0;
		return;
	}
	//It compiles to the same as data = *(uint64*)buffer; ->  mov	rax, QWORD PTR [rcx]
	//i am not sure about other platforms than x86
	memcpy(&data, buffer, sizeof(data)); 
	if (little_endian)
		data = _bswap_uint64(data);
	data >>= (sizeof(data)-n) << 3;
	buffer += n;
}

// *********************************************************************
inline void CKmer<1>::load(uchar *&buffer, int32 n)
{
	clear();

	if (!n)
		return;
	unsigned long long tmp;
	tmp = 0;

	switch (n-1)
	{
	case 7:
		data = ((uint64) *buffer++) << 56;
	case 6:
		tmp = ((uint64) *buffer++) << 48;
	case 5:
		data += ((uint64) *buffer++) << 40;
	case 4:
		tmp += ((uint64) *buffer++) << 32;
	case 3:
		data += ((uint64) *buffer++) << 24;
	case 2:
		tmp += ((uint64) *buffer++) << 16;
	case 1:
		data += ((uint64) *buffer++) << 8;
	}

	tmp += *buffer++;

	data += tmp;
}

// *********************************************************************
char CKmer<1>::get_symbol(int p)
{
	uint32 x = (data >> (2 * p)) & 0x03;

	switch (x)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	default: return 'T';
	}
}

// *********************************************************************
void CKmer<1>::set_prefix(CKmer<1>& rhs, uint32 /*suffix_bytes*/)
{
	data += rhs.data;
}
#endif

// ***** EOF


