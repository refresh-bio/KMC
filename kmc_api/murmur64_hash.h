/*******************************************************************************

 CoLoRd
 Copyright (C) 2021, M. Kokot, S. Deorowicz, and A. Gudys
 https://github.com/refresh-bio/CoLoRd

 This program is free software: you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this
 program. If not, see https://www.gnu.org/licenses/.

******************************************************************************/
#pragma once
#include <cstddef>
#include <cinttypes>

struct MurMur32Hash
{
	std::size_t operator()(uint32_t h) const noexcept
	{
		h ^= h >> 16;
		h *= 0x85ebca6b;
		h ^= h >> 13;
		h *= 0xc2b2ae35;
		h ^= h >> 16;

		return h;
	}

	/*	static inline uint32_t murmur_32_scramble(uint32_t k) {
			k *= 0xcc9e2d51;
			k = (k << 15) | (k >> 17);
			k *= 0x1b873593;
			return k;
		}
		uint32_t operator()(uint32_t key)
		{
			uint32_t h = 0;
			uint32_t k = key;

			key += sizeof(uint32_t);
			h ^= murmur_32_scramble(k);
			h = (h << 13) | (h >> 19);
			h = h * 5 + 0xe6546b64;

			h ^= murmur_32_scramble(k);

			h ^= 4;
			h ^= h >> 16;
			h *= 0x85ebca6b;
			h ^= h >> 13;
			h *= 0xc2b2ae35;
			h ^= h >> 16;

			return h;
		}*/
};

// MurMurHash3
struct MurMur64Hash
{
	std::size_t operator()(size_t h) const noexcept
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;

		return h;

		/*		h = (~h) + (h << 21); // key = (key << 21) - key - 1;
				h = h ^ (h >> 24);
				h = (h + (h << 3)) + (h << 8); // key * 265
				h = h ^ (h >> 14);
				h = (h + (h << 2)) + (h << 4); // key * 21
				h = h ^ (h >> 28);
				h = h + (h << 31);
				return h;*/
	}
};

struct CityHash64
	//struct MurMur64Hash
{
	static const uint64_t k2 = 0x9ae16a3b2f90404fULL;

	static uint64_t Rotate(uint64_t val, int shift) {
		// Avoid shifting by 64: doing so yields an undefined result.
		return shift == 0 ? val : ((val >> shift) | (val << (64 - shift)));
	}

	static uint64_t HashLen16(uint64_t u, uint64_t v, uint64_t mul) {
		// Murmur-inspired hashing.
		uint64_t a = (u ^ v) * mul;
		a ^= (a >> 47);
		uint64_t b = (v ^ a) * mul;
		b ^= (b >> 47);
		b *= mul;
		return b;
	}

	std::size_t operator()(size_t x) const noexcept
	{
		uint64_t mul = k2 + 16;
		uint64_t a = x + k2;
		uint64_t b = x;
		uint64_t c = Rotate(b, 37) * mul + a;
		uint64_t d = (Rotate(a, 25) + b) * mul;

		return HashLen16(c, d, mul);
	}

};