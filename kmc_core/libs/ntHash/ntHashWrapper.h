#pragma once
#include "ntHash.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifdef __GNUC__
//#include <intrin.h>
#endif

template<typename T> class cyclic_buffer
{
	std::vector<T> buf;
	size_t size;
	size_t size_m1;
	size_t no_elements;
	size_t idx_start;
	size_t idx_end;

	void _resize(size_t _size)
	{
		size = _size;
		size_m1 = size - 1u;
		buf.resize(size);
		clear();
	}

public:
	cyclic_buffer(size_t _size)
	{
		_resize(_size);
	}

	void resize(size_t _size)
	{
		_resize(_size);
	}

	void push(T x)
	{
//		if (no_elements == size)
//			return false;

		buf[idx_end] = x;
		if (++idx_end == size)
			idx_end = 0;

		++no_elements;

//		return true;
	}

	void pop(T& x)
	{
//		if (no_elements == 0)
//			return false;

		x = buf[idx_start];

		if (++idx_start == size)
			idx_start = 0;

		--no_elements;

//		return true;
	}

	void clear()
	{
		idx_start = 0;
		idx_end = 0;
		no_elements = 0;
	}

	size_t count()
	{
		return no_elements;
	}

	bool full()
	{
		return no_elements == size;
	}

	bool almost_full()
	{
		return no_elements == size_m1;
	}
};

class CntHashEstimator
{
	uint32_t k;
	uint32_t s, r;

	uint64_t insert_mask;
	const uint64_t accept_mask0 = 1ull;
	uint64_t accept_mask1;

	uint32_t* counters[2];
	uint32_t counters_size;
	const uint32_t max_occ_value = 65536;

	char convert[256] = {
		'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
	};

	const uint64_t* M_msTab33r[256] = {
	A33r, C33r, G33r, T33r, N33r, N33r, N33r, N33r, // 0..7
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 8..15
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 16..23
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 24..31
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 32..39
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 40..47
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 48..55
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 56..63
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 64..71
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 72..79
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 80..87
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 88..95
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 96..103
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 104..111
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 112..119
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 120..127
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 128..135
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 136..143
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 144..151
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 152..159
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 160..167
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 168..175
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 176..183
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 184..191
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 192..199
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 200..207
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 208..215
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 216..223
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 224..231
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 232..239
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 240..247
	N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r  // 248..255
	};

	const uint64_t* M_msTab31l[256] = {
		A31l, C31l, G31l, T31l, N31l, N31l, N31l, N31l, // 0..7
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 8..15
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 16..23
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 24..31
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 32..39
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 40..47
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 48..55
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 56..63
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 64..71
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 72..79
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 80..87
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 88..95
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 96..103
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 104..111
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 112..119
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 120..127
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 128..135
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 136..143
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 144..151
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 152..159
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 160..167
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 168..175
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 176..183
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 184..191
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 192..199
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 200..207
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 208..215
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 216..223
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 224..231
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 232..239
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 240..247
		N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l  // 248..255
	};

	const uint64_t M_seedTab[256] = {
	seedA, seedC, seedG, seedT, seedN, seedN, seedN, seedN, // 0..7
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 40..47
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 48..55
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 64..71
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 72..79
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 80..87
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 96..103
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 112..119
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 120..127
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 128..135
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 136..143
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 144..151
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 152..159
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 160..167
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 168..175
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 176..183
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 184..191
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 192..199
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 200..207
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 208..215
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 216..223
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 224..231
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 232..239
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 240..247
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN  // 248..255
	};

	uint64_t P_msTab[256];
	uint64_t P_msTab_seedTab[4][4];

	// **************************************************
	// rotate "v" to the right by 1 position
	inline uint64_t M_ror1(const uint64_t v) {
		return (v >> 1) | (v << 63);
//		return _rotr64(v, 1);
	}

	// **************************************************
	// rotate "v" to the left 1 position
	inline uint64_t M_rol1(const uint64_t v) {
		return (v << 1) | (v >> 63);
//		return _rotl64(v, 1);
	}

	// **************************************************
	// forward-strand ntHash for sliding k-mers
	inline uint64_t M_NTF64(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn) {
		uint64_t hVal = M_rol1(fhVal);
		hVal = swapbits033(hVal);
		hVal ^= P_msTab_seedTab[charOut][charIn];
//		hVal ^= M_seedTab[charIn];
//		hVal ^= P_msTab[charOut];
		return hVal;
	}

	// **************************************************
	// forward-strand ntHash for sliding k-mers
	inline uint64_t M_NTF64_N(const uint64_t fhVal, const unsigned char charIn) {
		uint64_t hVal = M_rol1(fhVal);
		hVal = swapbits033(hVal);
		hVal ^= M_seedTab[charIn];
		return hVal;
	}

	// **************************************************
	// reverse-complement ntHash for sliding k-mers
	inline uint64_t M_NTR64(const uint64_t rhVal, const unsigned char charOut, const unsigned char charIn) {
//		uint64_t hVal = rhVal ^ P_msTab[3 - charIn];
//		hVal ^= M_seedTab[3 - charOut];
		uint64_t hVal = rhVal ^ P_msTab_seedTab[3 - charIn][3 - charOut];
		hVal = ror1(hVal);
		hVal = swapbits3263(hVal);
		return hVal;
	}

	// **************************************************
	// reverse-complement ntHash for sliding k-mers
	inline uint64_t M_NTR64_N(const uint64_t rhVal, const unsigned char charIn) {
		uint64_t hVal = rhVal ^ P_msTab[3 - charIn];
		hVal = M_ror1(hVal);
		hVal = swapbits3263(hVal);
		return hVal;
	}

	// **************************************************
	inline uint64_t M_NTC64(const unsigned char charOut, const unsigned char charIn, uint64_t& fhVal, uint64_t& rhVal) {
		fhVal = M_NTF64(fhVal, charOut, charIn);
		rhVal = M_NTR64(rhVal, charOut, charIn);
		return (rhVal < fhVal) ? rhVal : fhVal;
	}

	// **************************************************
	inline uint64_t M_NTC64_N(const unsigned char charIn, uint64_t& fhVal, uint64_t& rhVal) {
		fhVal = M_NTF64_N(fhVal, charIn);
		rhVal = M_NTR64_N(rhVal, charIn);
		return (rhVal < fhVal) ? rhVal : fhVal;
	}

	// **************************************************
	void try_insert(uint64_t x)
	{
		uint64_t pref = x >> (63 - s);

		if (pref == accept_mask0)
			insert(0, x & insert_mask);
		
		pref >>= 1;
		if (pref == accept_mask1)
			insert(1, x & insert_mask);
	}

	// **************************************************
	void insert(uint32_t type, uint64_t x)
	{
#ifdef _MSC_VER
		_InterlockedIncrement((long*) &counters[type][x]);
#endif
#ifdef __GNUC__
		__atomic_add_fetch(&counters[type][x], 1, __ATOMIC_RELAXED);
#endif
	}

public:
	// **************************************************
	CntHashEstimator(uint32_t _k, uint32_t _s = 11, uint32_t _r = 27) : k(_k), s(_s), r(_r)
	{
		counters_size = 1ull << r;

		counters[0] = new uint32_t[counters_size];
		counters[1] = new uint32_t[counters_size];

		std::fill_n(counters[0], counters_size, 0ull);
		std::fill_n(counters[1], counters_size, 0ull);

		accept_mask1 = (1ull << (s - 1)) - 1ull;

		insert_mask = ~0ull >> (64u - r);

		uint32_t k31 = k % 31;
		uint32_t k33 = k % 33;

		for(int i = 0; i < 256; ++i)
			P_msTab[i] = M_msTab31l[i][k31] | M_msTab33r[i][k33];

		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j)
				P_msTab_seedTab[i][j] = P_msTab[i] ^ M_seedTab[j];
	}

	// **************************************************
	~CntHashEstimator()
	{
		delete[] counters[0];
		delete[] counters[1];
	}

	// **************************************************
	void ProcessCanonical(const char* str, uint32_t len)
	{
		if (len < k)
			return;

		uint64_t hVal;
		uint64_t fhVal = 0;
		uint64_t rhVal = 0;

		uint8_t c_out;
		uint8_t c_in;

		cyclic_buffer<uint8_t> c_buf(k);

		for (uint32_t i = 0; i < len; ++i)
		{
			c_in = (uint8_t) str[i];

			if (c_in == 255)		// N
			{
				c_buf.clear();
				fhVal = rhVal = 0;
			}
			else
			{
				if (c_buf.full())
				{
					c_buf.pop(c_out);

					hVal = M_NTC64(c_out, c_in, fhVal, rhVal);
					try_insert(hVal);
				}
				else
				{
					hVal = M_NTC64_N(c_in, fhVal, rhVal);

					if (c_buf.almost_full())
						try_insert(hVal);
				}

				c_buf.push(c_in);
			}
		}
	}

	void ProcessDirect(const char* str, uint32_t len)
	{
		if (len < k)
			return;

		uint64_t fhVal = 0;

		uint8_t c_out;
		uint8_t c_in;

		cyclic_buffer<uint8_t> c_buf(k);

		for (uint32_t i = 0; i < len; ++i)
		{
			c_in = (uint8_t) str[i];

			if (c_in == 255)		// N
			{
				c_buf.clear();
				fhVal = 0;
			}
			else
			{
				if (c_buf.full())
				{
					c_buf.pop(c_out);

					fhVal = M_NTF64(fhVal, c_out, c_in);
					try_insert(fhVal);
				}
				else
				{
					fhVal = M_NTF64_N(fhVal, c_in);

					if (c_buf.almost_full())
						try_insert(fhVal);
				}

				c_buf.push(c_in);
			}
		}
	}

	void EstimateHistogram(std::vector<uint64_t>& v_histogram)
	{
		uint32_t max_occ0 = *std::max_element(counters[0], counters[0] + counters_size);
		uint32_t max_occ1 = *std::max_element(counters[1], counters[1] + counters_size);
		uint32_t max_occ = std::max(max_occ0, max_occ1);
		max_occ = std::min(std::max(max_occ, 2u), max_occ_value);

		std::vector<uint32_t> v_hist0(max_occ + 1u);
		std::vector<uint32_t> v_hist1(max_occ + 1u);

		for (uint32_t i = 0; i < counters_size; ++i)
			++v_hist0[std::min(max_occ_value, counters[0][i])];

		for (uint32_t i = 0; i < counters_size; ++i)
			++v_hist1[std::min(max_occ_value, counters[1][i])];

		std::vector<double> v_mean(max_occ + 1u);
		for (uint32_t i = 0; i <= max_occ; ++i)
			v_mean[i] = (v_hist0[i] + v_hist1[i]) / 2.0;

		double c_log2 = log(2);
		double log_v_mean0 = log(v_mean[0]);
		double F0Mean = (int64_t) (((double) r * c_log2 - log_v_mean0) * (1ull << (s + r)));

		v_histogram.clear();
		v_histogram.resize(max_occ + 1, 0);

		if (v_mean[0] * (log_v_mean0 - r * c_log2) == 0)
			return;

		std::vector<double> v_hist;
		v_hist.resize(max_occ + 1, 0.0);

		v_hist[1] = -1.0 * v_mean[1] / (v_mean[0] * (log_v_mean0 - r * c_log2));
		for (uint32_t i = 2; i <= max_occ; ++i)
		{
			double sum = 0.0;
			for (uint32_t j = 1; j < i; ++j)
				sum += j * v_mean[i - j] * v_hist[j];
			v_hist[i] = -1.0 * v_mean[i] / (v_mean[0] * (log_v_mean0 - r * c_log2)) - sum / (i * v_mean[0]);
		}

		for (uint32_t i = 1; i <= max_occ; ++i)
			v_histogram[i] = abs((int64_t) (v_hist[i] * F0Mean));
	}
};

// EOF
