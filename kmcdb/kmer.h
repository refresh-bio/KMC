#ifndef KMER_H_
#define KMER_H_
#include <cstdint>
#include <cassert>
#include <emmintrin.h>

#ifdef _WIN32
#define _bswap64(x) _byteswap_uint64(x)
#else
#define _bswap64(x) __builtin_bswap64(x)
#endif

namespace kmcdb
{
	inline uint32_t get_rev_compl_shift(uint32_t len)
	{
		return (len + 31) / 32 * 32 - len;
	}

	//instead of x << 2*p there is (x<<p) << p, because
	//for p=32 we would have x << 64 which is UB
	inline uint64_t shl_2p(uint64_t x, uint64_t p) {
		return (x << p) << p;
	}
	//instead of x >> 2*p there is (x>>p) >> p, because
	//for p=32 we would have x >> 64 which is UB
	inline uint64_t shr_2p(uint64_t x, uint64_t p) {
		return (x >> p) >> p;
	}
// *************************************************************************
// Ckmer class for k > 32 with classic kmer counting
template<unsigned SIZE> struct CKmer {
	unsigned long long data[SIZE];

	typedef unsigned long long data_t;
	static uint32_t KMER_SIZE;

	inline void set(const CKmer<SIZE>& x);

	inline void from_kxmer(const CKmer<SIZE>& x, uint32_t _shr, const CKmer<SIZE>& _mask);

	template<unsigned X_SIZE> inline void to_kxmer(CKmer<X_SIZE>& x);

	inline void mask(const CKmer<SIZE>& x);
	inline uint32_t end_mask(const uint32_t mask);
	inline void set_2bits(const uint64_t x, const uint32_t p);
	inline uint8_t get_2bits(const uint32_t p) const;
	inline uint8_t get_byte(const uint32_t p) const;
	inline void set_byte(const uint32_t p, uint8_t x);
	inline void set_bits(const uint32_t p, const uint32_t n, uint64_t x);

	inline void SHL_insert_2bits(const uint64_t x);
	inline void SHR_insert_2bits(const uint64_t x, const uint32_t p);

	inline void SHR(const uint32_t p);
	inline void SHL(const uint32_t p);

	inline uint64_t remove_suffix(const uint32_t n) const;

	inline void set_prefix(uint64_t prefix, const uint32_t n);
	inline void set_n_1(const uint32_t n);
	inline void set_n_01(const uint32_t n);

	//mkokot_TODO: read comment in comment in CKmer<1>
	void store_left_aligned(uint8_t*& buffer, uint64_t k) const;

	//mkokot_TODO: read comment in comment in CKmer<1>
	void load_from_left_aligned(const uint8_t*& buffer, uint64_t k);

	inline void store(uint8_t*& buffer, int32_t n);
	inline void store(uint8_t* buffer, int32_t p, int32_t n);

	inline void load(uint8_t*& buffer, int32_t n);

	inline bool operator==(const CKmer<SIZE>& x) const;
	inline bool operator<(const CKmer<SIZE>& x) const;

	inline void clear(void);

	inline char get_symbol(int p);

	inline void fill_T();

	inline void random_init(uint32_t pos, uint64_t value);

	inline void to_string(uint64_t kmer_len, char* out) const;

	inline std::string to_string(uint64_t kmer_len) const;

	inline void from_string(const std::string& str);

	inline unsigned long long& operator[](int idx)
	{
		return data[idx];
	}

	//similar tu substring
	//n must be <= 32
	//p = 0 takes the least significant symbols of k-mer
	//p is in symbols (p = 1 means skip two least significant bits)
	inline uint64_t subkmer(uint64_t n, uint64_t p) const
	{
		assert(n + p <= SIZE * 32);

		uint64_t result = data[p >> 5] >> ((2 * p) & 63);

		if ((p >> 5) != ((p + n - 1) >> 5))
			result += data[(p >> 5) + 1] << (64 - ((2 * p) & 63));

		result &= ((1ull << n) << n) - 1ull;

		return result;
	}


	CKmer<SIZE> GetRevCompl(const uint32_t shift) const
	{
		CKmer<SIZE> res;

		const uint64_t mask4 = 0b0000111100001111000011110000111100001111000011110000111100001111ull;
		const uint64_t mask2 = 0b0011001100110011001100110011001100110011001100110011001100110011ull;

		const __m128i v_mask4 = _mm_set1_epi64x(mask4);
		const __m128i v_mask2 = _mm_set1_epi64x(mask2);
		const __m128i v_mask_xor = _mm_set1_epi64x(-1);

		for (uint32_t i = 0; i + 1 < SIZE; i += 2)
		{
			__m128i d = _mm_set_epi64x(_bswap64(data[i]), _bswap64(data[i + 1]));
			__m128i sf4 = _mm_add_epi64(_mm_slli_epi64(_mm_and_si128(d, v_mask4), 4), _mm_and_si128(_mm_srli_epi64(d, 4), v_mask4));
			__m128i sf2 = _mm_add_epi64(_mm_slli_epi64(_mm_and_si128(sf4, v_mask2), 2), _mm_and_si128(_mm_srli_epi64(sf4, 2), v_mask2));

			sf2 = _mm_xor_si128(sf2, v_mask_xor);

			_mm_storeu_si128((__m128i*) (res.data + SIZE - 2 - i), sf2);
		}

		if constexpr (SIZE % 2 == 1)
		{
			uint64_t sf4 = ((data[SIZE - 1] & mask4) << 4) + ((data[SIZE - 1] >> 4) & mask4);
			uint64_t sf2 = ((sf4 & mask2) << 2) + ((sf4 >> 2) & mask2);

			res.data[0] = ~_bswap64(sf2);
		}

		res.SHR(shift);

		return res;
	}

};

template <unsigned SIZE> uint32_t CKmer<SIZE>::KMER_SIZE = SIZE;

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set(const CKmer<SIZE>& x)
{
	for (uint32_t i = 0; i < SIZE; ++i)
		data[i] = x.data[i];
}


// *********************************************************************
template<unsigned SIZE>
template<unsigned X_SIZE> inline void CKmer<SIZE>::to_kxmer(CKmer<X_SIZE>& x)
{
	x.data[X_SIZE - 1] = 0;
	for (uint32_t i = 0; i < SIZE; ++i)
		x.data[i] = data[i];
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::from_kxmer(const CKmer<SIZE>& x, uint32_t _shr, const CKmer<SIZE>& _mask)
{
	if (_shr)
	{
		for (uint32_t i = 0; i < SIZE - 1; ++i)
		{
			data[i] = x.data[i] >> (2 * _shr);
			data[i] += x.data[i + 1] << (64 - 2 * _shr);
			data[i] &= _mask.data[i];
		}
		data[SIZE - 1] = x.data[SIZE - 1] >> (2 * _shr);
		data[SIZE - 1] &= _mask.data[SIZE - 1];
	}
	else
	{
		for (uint32_t i = 0; i < SIZE; ++i)
			data[i] = x.data[i] & _mask.data[i];
	}
	//	mask(_mask);
}


// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::mask(const CKmer<SIZE>& x)
{
	for (uint32_t i = 0; i < SIZE; ++i)
		data[i] &= x.data[i];
}

// *********************************************************************
template<unsigned SIZE> inline uint32_t CKmer<SIZE>::end_mask(const uint32_t mask)
{
	return data[0] & mask;
}
// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_2bits(const uint64_t x, const uint32_t p)
{
	data[p >> 6] += x << (p & 63);
}

template<unsigned SIZE> inline uint8_t CKmer<SIZE>::get_2bits(const uint32_t p) const
{
	return (data[p >> 6] >> (p & 63)) & 3;
}
// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHR_insert_2bits(const uint64_t x, const uint32_t p)
{
	for (uint32_t i = 0; i < SIZE - 1; ++i)
	{
		data[i] >>= 2;
		data[i] += data[i + 1] << (64 - 2);
	}

	data[SIZE - 1] >>= 2;
	data[p >> 6] += x << (p & 63);
}



// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHR(const uint32_t p)
{
	for (uint32_t i = 0; i < SIZE - 1; ++i)
	{
		data[i] >>= 2 * p;
		data[i] += data[i + 1] << (64 - 2 * p);
	}
	data[SIZE - 1] >>= 2 * p;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHL(const uint32_t p)
{
	for (uint32_t i = SIZE - 1; i > 0; --i)
	{
		data[i] <<= p * 2;
		data[i] += data[i - 1] >> (64 - p * 2);
	}
	data[0] <<= p * 2;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::SHL_insert_2bits(const uint64_t x)
{
	for (uint32_t i = SIZE - 1; i > 0; --i)
	{
		data[i] <<= 2;
		data[i] += data[i - 1] >> (64 - 2);
	}
	data[0] <<= 2;
	data[0] += x;
}

// *********************************************************************
template<unsigned SIZE> inline uint8_t CKmer<SIZE>::get_byte(const uint32_t p) const
{
	return (data[p >> 3] >> ((p << 3) & 63)) & 0xFF;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_byte(const uint32_t p, uint8_t x)
{
	data[p >> 3] += ((uint64_t)x) << ((p & 7) << 3);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_bits(const uint32_t p, const uint32_t n, uint64_t x)
{
	data[p >> 6] += x << (p & 63);
	if ((p >> 6) != ((p + n - 1) >> 6))
		data[(p >> 6) + 1] += x >> (64 - (p & 63));
}

// *********************************************************************
template<unsigned SIZE> inline bool CKmer<SIZE>::operator==(const CKmer<SIZE>& x) const {
	for (uint32_t i = 0; i < SIZE; ++i)
		if (data[i] != x.data[i])
			return false;

	return true;
}

// *********************************************************************
template<unsigned SIZE> inline bool CKmer<SIZE>::operator<(const CKmer<SIZE>& x) const {
	for (int32_t i = SIZE - 1; i >= 0; --i)
		if (data[i] < x.data[i])
			return true;
		else if (data[i] > x.data[i])
			return false;
	return false;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::clear(void)
{
	for (uint32_t i = 0; i < SIZE; ++i)
		data[i] = 0;
}

// *********************************************************************
template<unsigned SIZE> inline uint64_t CKmer<SIZE>::remove_suffix(const uint32_t n) const
{
	uint32_t p = n >> 6; // / 64;
	uint32_t r = n & 63;	// % 64;

	if (p == SIZE - 1)
		return data[p] >> r;
	else
		return (data[p + 1] << (64 - r)) + (data[p] >> r);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_prefix(uint64_t prefix, const uint32_t n)
{
	uint32_t p = n >> 6; // / 64;
	uint32_t r = n & 63;	// % 64;

	if (p == SIZE - 1)
		data[p] += prefix << r;
	else
	{
		data[p + 1] += prefix >> (64 - r);
		data[p] += prefix << r;
	}
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_n_1(const uint32_t n)
{
	clear();

	for (uint32_t i = 0; i < (n >> 6); ++i)
		data[i] = ~((uint64_t)0);

	uint32_t r = n & 63;

	if (r)
		data[n >> 6] = (1ull << r) - 1;
}


// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::set_n_01(const uint32_t n)
{
	clear();

	for (uint32_t i = 0; i < n; ++i)
		if (!(i & 1))
			data[i >> 6] += (1ull << (i & 63));
}

// *********************************************************************
template<unsigned SIZE> void CKmer<SIZE>::store_left_aligned(uint8_t*& buffer, uint64_t k) const
{
	//first make a copy to work on, after changing the representation the copy will not be needed
	CKmer<SIZE> copy = *this;

	//this could be precomputed
	//if k % 32 != 0 we need to shl first
	if (k & 31)
		copy.SHL(32 - (k & 31)); //align toward most significant
	if ((k + 31) / 32 != SIZE) //if we store prefix a suffix separately we may pass k that is lower than it may seem from SIZE
	{
		//mkokot_TODO: this is quite ugly, but lets make if work first
		copy.SHL(16);
		copy.SHL(16);
	}

	uint32_t bytes = static_cast<uint32_t>((k + 3) / 4);

	for (uint32_t pos = 0 ; pos < bytes ; ++pos)
		*buffer++ = copy.get_byte(SIZE * 8 - 1 - pos);
}


template<unsigned SIZE> void CKmer<SIZE>::load_from_left_aligned(const uint8_t*& buffer, uint64_t k)
{
	//mkokot_TODO: is clear needed, probably yes
	clear();

	uint32_t bytes = static_cast<uint32_t>((k + 3) / 4);

	for (uint32_t pos = 0 ; pos < bytes ; ++pos)
		set_byte(SIZE * 8 - 1 - pos, *buffer++);

	if (k & 31)
		SHR(32 - (k & 31)); //alignt toward least significant
	if((k+31)/32 != SIZE)
	{
		//mkokot_TODO: this is quite ugly, but lets make if work first
		SHR(16);
		SHR(16);
	}
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::store(uint8_t*& buffer, int32_t n)
{
	for (int32_t i = n - 1; i >= 0; --i)
		*buffer++ = get_byte(i);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::store(uint8_t* buffer, int32_t p, int32_t n)
{
	for (int32_t i = n - 1; i >= 0; --i)
		buffer[p++] = get_byte(i);
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::load(uint8_t*& buffer, int32_t n)
{
	clear();
	for (int32_t i = n - 1; i >= 0; --i)
		set_byte(i, *buffer++);
}

// *********************************************************************
template<unsigned SIZE> inline char CKmer<SIZE>::get_symbol(int p)
{
	uint32_t x = (data[p >> 5] >> (2 * (p & 31))) & 0x03;

	switch (x)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	default: return 'T';
	}
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::fill_T()
{
	for (uint32_t i = 0; i < SIZE; ++i)
		data[i] = ~0ull;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::random_init(uint32_t pos, uint64_t value)
{
	data[pos] = value;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::to_string(uint64_t kmer_len, char* out) const
{
	auto pos = static_cast<uint32_t>(2 * kmer_len - 2);
	for (uint32_t i = 0; i < kmer_len; ++i, pos -= 2)
		out[i] = "ACGT"[get_2bits(pos)];
}

// *********************************************************************
template<unsigned SIZE> inline std::string CKmer<SIZE>::to_string(uint64_t kmer_len) const
{
	std::string res(kmer_len, ' ');
	to_string(kmer_len, res.data());
	return res;
}

// *********************************************************************
template<unsigned SIZE> inline void CKmer<SIZE>::from_string(const std::string& str)
{
	//mkokot_TODO: should probably be optimized
	clear();
	if (str.empty())
		return;

	auto pos = static_cast<uint32_t>(2 * str.length() - 2);
	for (uint32_t i = 0; i < str.length(); ++i, pos -= 2)
	{
		switch (str[i])
		{
		case 'A': set_2bits(0, pos); break;
		case 'C': set_2bits(1, pos); break;
		case 'G': set_2bits(2, pos); break;
		case 'T': set_2bits(3, pos); break;
		}
	}
}

// *********************************************************************
// *********************************************************************
// *********************************************************************
// *********************************************************************
// Ckmer class for k <= 32 with classic kmer counting
template<> struct CKmer<1> {
	unsigned long long data;

	typedef unsigned long long data_t;
	static uint32_t KMER_SIZE;

	void set(const CKmer<1>& x);

	void from_kxmer(const CKmer<1>& x, uint32_t _shr, const CKmer<1>& _mask);

	template <unsigned X_SIZE> void to_kxmer(CKmer<X_SIZE>& x);

	void mask(const CKmer<1>& x);
	uint32_t end_mask(const uint32_t mask);
	void set_2bits(const uint64_t x, const uint32_t p);
	uint8_t get_2bits(const uint32_t p) const;
	uint8_t get_byte(const uint32_t p) const;
	void set_byte(const uint32_t p, uint8_t x);
	void set_bits(const uint32_t p, const uint32_t n, uint64_t x);

	void SHL_insert_2bits(const uint64_t x);
	void SHR_insert_2bits(const uint64_t x, const uint32_t p);

	void SHR(const uint32_t p);
	void SHL(const uint32_t p);

	uint64_t remove_suffix(const uint32_t n) const;

	inline void set_prefix(uint64_t prefix, const uint32_t n);

	void set_n_1(const uint32_t n);
	void set_n_01(const uint32_t n);

	void store(uint8_t*& buffer, int32_t n);

	//mkokot_TODO: this representation is in general right aligned
	//but we will switch to left aligned
	//for this reason I am adding this method to be used by the archive based storage
	//maybe its not the best way to pass k here, but compute what is needed outsize and use, like the number of bytes, etc.
	void store_left_aligned(uint8_t*& buffer, uint64_t k) const;

	//mkokot_TODO: similarly to store_left_aligned, the idea is that in this file we currently have
	//data right aligned, but we store as left aligned, so we need to make additional shr at the deserialization end
	void load_from_left_aligned(const uint8_t*& buffer, uint64_t k);

	void store(uint8_t* buffer, int32_t p, int32_t n);
	void load(uint8_t*& buffer, int32_t n);

	bool operator==(const CKmer<1>& x) const;
	bool operator<(const CKmer<1>& x)const;

	void clear(void);

	inline char get_symbol(int p);

	inline void fill_T();

	inline void random_init(uint32_t pos, uint64_t value);

	inline void to_string(uint64_t kmer_len, char* out) const;

	inline std::string to_string(uint64_t kmer_len) const;

	inline void from_string(const std::string& str);

	inline unsigned long long& operator[](int idx)
	{
		(void)idx;
		assert(idx == 0);
		return data;
	}

	//similar tu substring
	//n must be <= 32
	//p = 0 takes the least significant symbols of k-mer
	//p is in symbols (p = 1 means skip two least significant bits)
	inline uint64_t subkmer(uint64_t n, uint64_t p) const
	{
		assert(p + n <= 32);
		return (data >> (2ull * p)) & (((1ull << n) << n) - 1ull);
	}

	CKmer<1> GetRevCompl(const uint32_t shift) const
	{
		CKmer<1> res;

		const uint64_t mask4 = 0b0000111100001111000011110000111100001111000011110000111100001111ull;
		const uint64_t mask2 = 0b0011001100110011001100110011001100110011001100110011001100110011ull;

		uint64_t sf4 = ((data & mask4) << 4) + ((data >> 4) & mask4);
		uint64_t sf2 = ((sf4 & mask2) << 2) + ((sf4 >> 2) & mask2);

		res.data = shr_2p((~_bswap64(sf2)), shift);

		return res;
	}
};


// *********************************************************************
template <unsigned X_SIZE> inline void CKmer<1>::to_kxmer(CKmer<X_SIZE>& x)
{
	x.data[X_SIZE - 1] = 0;
	x.data[0] = data;
}

// *********************************************************************
template<> inline void CKmer<1>::to_kxmer(CKmer<1>& x)
{
	x.data = data;
}


// *********************************************************************
inline void CKmer<1>::mask(const CKmer<1>& x)
{
	data &= x.data;
}


// *********************************************************************
inline uint32_t CKmer<1>::end_mask(const uint32_t mask)
{
	return data & mask;
}
// *********************************************************************
inline void CKmer<1>::set(const CKmer<1>& x)
{
	data = x.data;
}

// *********************************************************************
inline void CKmer<1>::from_kxmer(const CKmer<1>& x, uint32_t _shr, const CKmer<1>& _mask)
{
	data = (x.data >> (2 * _shr)) & _mask.data;
}


// *********************************************************************
inline void CKmer<1>::set_2bits(const uint64_t x, const uint32_t p)
{
	data += x << p;
}

inline uint8_t CKmer<1>::get_2bits(const uint32_t p) const
{
	return (data >> p) & 3;
}
// *********************************************************************
inline void CKmer<1>::SHR_insert_2bits(const uint64_t x, const uint32_t p)
{
	data >>= 2;
	data += x << p;
}

// *********************************************************************
inline void CKmer<1>::SHR(const uint32_t p)
{
	data >>= 2 * p;
}

// *********************************************************************
inline void CKmer<1>::SHL(const uint32_t p)
{
	data <<= p * 2;
}
// *********************************************************************
inline void CKmer<1>::SHL_insert_2bits(const uint64_t x)
{
	data = (data << 2) + x;
}

// *********************************************************************
inline uint8_t CKmer<1>::get_byte(const uint32_t p) const
{
	return (data >> (p << 3)) & 0xFF;
}

// *********************************************************************
inline void CKmer<1>::set_byte(const uint32_t p, uint8_t x)
{
	data += ((uint64_t)x) << (p << 3);
}

// *********************************************************************
inline void CKmer<1>::set_bits(const uint32_t p, const uint32_t /*n*/, uint64_t x)
{
	data += x << p;
}

// *********************************************************************
inline bool CKmer<1>::operator==(const CKmer<1>& x) const {
	return data == x.data;
}

// *********************************************************************
inline bool CKmer<1>::operator<(const CKmer<1>& x) const {
	return data < x.data;
}

// *********************************************************************
inline void CKmer<1>::clear(void)
{
	data = 0ull;
}

// *********************************************************************
inline uint64_t CKmer<1>::remove_suffix(const uint32_t n) const
{
	return data >> n;
}

// *********************************************************************
inline void CKmer<1>::set_prefix(uint64_t prefix, const uint32_t n)
{
	data += prefix << n;
}
// *********************************************************************
inline void CKmer<1>::set_n_1(const uint32_t n)
{
	if (n == 64)
		data = ~(0ull);
	else
		data = (1ull << n) - 1;
}

// *********************************************************************
inline void CKmer<1>::set_n_01(const uint32_t n)
{
	data = 0ull;

	for (uint32_t i = 0; i < n; ++i)
		if (!(i & 1))
			data += (1ull << i);
}

// *********************************************************************
inline void CKmer<1>::store(uint8_t*& buffer, int32_t n)
{
	for (int32_t i = n - 1; i >= 0; --i)
		*buffer++ = get_byte(i);
}
// *********************************************************************
inline void CKmer<1>::store_left_aligned(uint8_t*& buffer, uint64_t k) const
{
	//mkokot_TODO: first make a copy to work on, after changing the representation the copy will not be needed
	CKmer<1> copy(*this);

	//this could be precomputed
	//if k % 32 != 0 we need to shl first
	if (k & 31)
		copy.SHL(32 - (k & 31)); //align toward most significant
	uint32_t bytes = static_cast<uint32_t>((k + 3) / 4);

	for (uint32_t pos = 0 ; pos < bytes ; ++pos)
		*buffer++ = copy.get_byte(7 - pos);
}

// *********************************************************************
inline void CKmer<1>::load_from_left_aligned(const uint8_t*& buffer, uint64_t k)
{
	//mkokot_TODO: is clear needed, probably yes
	clear();

	uint32_t bytes = static_cast<uint32_t>((k + 3) / 4);

	for (uint32_t pos = 0 ; pos < bytes ; ++pos)
		set_byte(7 - pos, *buffer++);

	if (k & 31)
		SHR(32 - (k & 31)); //alignt toward least significant
}
// *********************************************************************
inline void CKmer<1>::store(uint8_t* buffer, int32_t p, int32_t n)
{
	for (int32_t i = n - 1; i >= 0; --i)
		buffer[p++] = get_byte(i);
}

// *********************************************************************
inline void CKmer<1>::load(uint8_t*& buffer, int32_t n)
{
	clear();
	for (int32_t i = n - 1; i >= 0; --i)
		set_byte(i, *buffer++);
}


// *********************************************************************
char CKmer<1>::get_symbol(int p)
{
	uint32_t x = (data >> (2 * p)) & 0x03;

	switch (x)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	default: return 'T';
	}
}

// *********************************************************************
inline void CKmer<1>::fill_T()
{
	data = ~0ull;
}

// *********************************************************************
inline void CKmer<1>::random_init(uint32_t /*pos*/, uint64_t value)
{
	data = value;
}

// *********************************************************************
inline void CKmer<1>::to_string(uint64_t kmer_len, char* out) const
{
	auto pos = static_cast<uint32_t>(2 * kmer_len - 2);
	for (uint32_t i = 0; i < kmer_len; ++i, pos -= 2)
		out[i] = "ACGT"[get_2bits(pos)];
}

// *********************************************************************
inline std::string CKmer<1>::to_string(uint64_t kmer_len) const
{
	std::string res(kmer_len, ' ');
	to_string(kmer_len, res.data());
	return res;
}

// *********************************************************************
inline void CKmer<1>::from_string(const std::string& str)
{
	//mkokot_TODO: should probably be optimized
	clear();
	if (str.empty())
		return;

	auto pos = static_cast<uint32_t>(2 * str.length() - 2);
	for (uint32_t i = 0; i < str.length(); ++i, pos -= 2)
	{
		switch (str[i])
		{
		case 'A': set_2bits(0, pos); break;
		case 'C': set_2bits(1, pos); break;
		case 'G': set_2bits(2, pos); break;
		case 'T': set_2bits(3, pos); break;
		}
	}
}

}

#undef _bswap64

#endif


