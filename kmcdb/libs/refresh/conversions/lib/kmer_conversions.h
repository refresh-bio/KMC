#pragma once

#include <cinttypes>
#include <string>
#include <cassert>

namespace refresh
{
	namespace conversions
	{
		// ************************************************************************************
		class kmer_conversions
		{
			static const size_t bases_per_chunk = 6;										// max is 10
			static const size_t dict_size = (1ull << (2 * bases_per_chunk));
			static char inline bases[dict_size * bases_per_chunk];
			static const uint64_t chunk_shift_to_right = 64 - 2 * bases_per_chunk;
			static const uint64_t chunk_shift_to_left = 2 * bases_per_chunk;

			struct _si {
				_si()
				{
					size_t idx = 0;
					for (size_t i = 0; i < dict_size; ++i)
					{
						uint64_t x = (uint64_t)i;

						for (size_t j = bases_per_chunk; j > 0; --j)
						{
							bases[idx + j - 1] = "ACGT"[x & 3];
							x >>= 2;
						}

						idx += bases_per_chunk;
					}
				}
			} static inline _init;

			static void copy_chunk_bases(char* dest, char* src)
			{
				static_assert(bases_per_chunk <= 10);

				if constexpr (bases_per_chunk > 0)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 1)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 2)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 3)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 4)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 5)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 6)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 7)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 8)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 9)
					*dest++ = *src++;
				if constexpr (bases_per_chunk > 10)
					*dest++ = *src++;
			}

			static void copy_bases(char* dest, char* src, size_t no_bases)
			{
				assert(bases_per_chunk <= 10);

				if (no_bases > 0)	*dest++ = *src++;		else return;
				if (no_bases > 1)	*dest++ = *src++;		else return;
				if (no_bases > 2)	*dest++ = *src++;		else return;
				if (no_bases > 3)	*dest++ = *src++;		else return;
				if (no_bases > 4)	*dest++ = *src++;		else return;
				if (no_bases > 5)	*dest++ = *src++;		else return;
				if (no_bases > 6)	*dest++ = *src++;		else return;
				if (no_bases > 7)	*dest++ = *src++;		else return;
				if (no_bases > 8)	*dest++ = *src++;		else return;
				if (no_bases > 9)	*dest++ = *src++;		else return;
			}

		public:
			static size_t kmer_to_pchar(uint64_t kmer, char* out, size_t no_bases, bool left_align = true)
			{
				if (!left_align)
					kmer <<= (64 - 2 * no_bases);

				auto left_bases = no_bases;
				for (; left_bases >= bases_per_chunk; left_bases -= bases_per_chunk)
				{
					uint64_t chunk = kmer >> chunk_shift_to_right;
					copy_chunk_bases(out, &bases[bases_per_chunk * chunk]);
					kmer <<= chunk_shift_to_left;

					out += bases_per_chunk;
				}

				uint64_t chunk = kmer >> chunk_shift_to_right;
				copy_bases(out, &bases[bases_per_chunk * chunk], left_bases);

				return no_bases;
			}

			static size_t kmer_to_pchar(uint64_t* kmer, char* out, size_t no_bases, bool left_align = true)
			{
				size_t x = no_bases % 32;
				if (x == 0)
					x = 32;

				if (!left_align)
				{
					kmer_to_pchar(kmer[(no_bases + 31) / 32 - 1], out, x, false);
					out += x;

					for (std::ptrdiff_t i = (no_bases + 31) / 32 - 2; i >= 0; --i, out += 32)
						kmer_to_pchar(kmer[i], out, 32);
				}
				else
				{
					for (size_t i = (no_bases + 31) / 32 - 1; i > 0 ; --i, out += 32)
						kmer_to_pchar(kmer[i], out, 32);

					kmer_to_pchar(kmer[0], out, x);
				}

				return no_bases;
			}
		};
	}

	// ************************************************************************************
	inline size_t kmer_to_pchar(uint64_t kmer, char* out, size_t no_bases, bool left_align = true, char term = '\0')
	{
		assert(no_bases <= 32);

		auto r = conversions::kmer_conversions::kmer_to_pchar(kmer, out, no_bases, left_align);
		out[r] = term;

		return r + 1;
	}

	// ************************************************************************************
	inline size_t kmer_to_pchar(uint64_t* kmer, char* out, size_t no_bases, bool left_align = true, char term = '\0')
	{
		auto r = conversions::kmer_conversions::kmer_to_pchar(kmer, out, no_bases, left_align);
		out[r] = term;

		return r + 1;
	}

	// ************************************************************************************
	inline std::string kmer_to_string(uint64_t kmer, size_t no_bases, bool left_align = true)
	{
		assert(no_bases <= 32);

		std::string out;
		out.resize(no_bases);

		conversions::kmer_conversions::kmer_to_pchar(kmer, out.data(), no_bases, left_align);

		return out;
	}

	// ************************************************************************************
	inline std::string kmer_to_string(uint64_t* kmer, size_t no_bases, bool left_align = true)
	{
		std::string out;
		out.resize(no_bases);

		conversions::kmer_conversions::kmer_to_pchar(kmer, out.data(), no_bases, left_align);

		return out;
	}
};
