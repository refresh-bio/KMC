#ifndef MMER_H_
#define MMER_H_

#include <cstdint>
#include <algorithm>
#include <cassert>
namespace kmcdb
{
	//mkokot_TODO: consider instead of looking at lookup tables
	//to just perform normal comparison and only if the value is lower validate if its allowed
	//maybe a little complex so probably must be done carefully and tested well
	class MmerSignature
	{
		uint32_t str{};
		uint32_t mask{};
		uint32_t current_val{};
		uint32_t* norm{};
		uint32_t len{};
		static inline uint32_t norm4[1 << 8];
		static inline uint32_t norm5[1 << 10];
		static inline uint32_t norm6[1 << 12];
		static inline uint32_t norm7[1 << 14];
		static inline uint32_t norm8[1 << 16];
		static inline uint32_t norm9[1 << 18];
		static inline uint32_t norm10[1 << 20];
		static inline uint32_t norm11[1 << 22];

		friend class CSignatureMapper;
		struct _si
		{
			static uint32_t get_rev(uint32_t mmer, uint32_t len)
			{
				uint32_t rev = 0;
				uint32_t shift = len * 2 - 2;
				for (uint32_t i = 0; i < len; ++i)
				{
					rev += (3 - (mmer & 3)) << shift;
					mmer >>= 2;
					shift -= 2;
				}
				return rev;
			}

			static void init_norm(uint32_t* norm, uint32_t len)
			{
				uint32_t special = 1 << len * 2;
				for (uint32_t i = 0; i < special; ++i)
				{
					uint32_t rev = get_rev(i, len);
					uint32_t str_val = is_allowed(i, len) ? i : special;
					uint32_t rev_val = is_allowed(rev, len) ? rev : special;
					norm[i] = (std::min)(str_val, rev_val);
				}
			}

			_si()
			{
				init_norm(norm4, 4);
				init_norm(norm5, 5);
				init_norm(norm6, 6);
				init_norm(norm7, 7);
				init_norm(norm8, 8);
				init_norm(norm9, 9);
				init_norm(norm10, 10);
				init_norm(norm11, 11);
			}

		} static inline _init;
	public:
		static bool is_allowed(uint32_t mmer, uint32_t len)
		{
			if ((mmer & 0x3f) == 0x3f)            // TTT suffix
				return false;
			if ((mmer & 0x3f) == 0x3b)            // TGT suffix
				return false;
			if ((mmer & 0x3c) == 0x3c)            // TG* suffix !!!! consider issue #152
				return false;

			for (uint32_t j = 0; j < len - 3; ++j)
				if ((mmer & 0xf) == 0)                // AA inside
					return false;
				else
					mmer >>= 2;

			if (mmer == 0)            // AAA prefix
				return false;
			if (mmer == 0x04)        // ACA prefix
				return false;
			if ((mmer & 0xf) == 0)    // *AA prefix
				return false;

			return true;
		}
		MmerSignature(uint32_t _len)
		{
			switch (_len)
			{
			case 4:
				norm = norm4;
				break;
			case 5:
				norm = norm5;
				break;
			case 6:
				norm = norm6;
				break;
			case 7:
				norm = norm7;
				break;
			case 8:
				norm = norm8;
				break;
			case 9:
				norm = norm9;
				break;
			case 10:
				norm = norm10;
				break;
			case 11:
				norm = norm11;
				break;
			default:
				assert(false); //this len is not supported
				break;
			}
			len = _len;
			mask = (1 << _len * 2) - 1;
			str = 0;
		}
		inline void insert(uint8_t symb)
		{
			str <<= 2;
			str += symb;
			str &= mask;

			current_val = norm[str];
		}
		inline uint32_t get() const
		{
			return current_val;
		}
		inline bool operator==(const MmerSignature& x) const
		{
			return current_val == x.current_val;
		}
		inline bool operator<(const MmerSignature& x) const
		{
			return current_val < x.current_val;
		}
		inline void clear()
		{
			str = 0;
		}
		inline bool operator<=(const MmerSignature& x) const
		{
			return current_val <= x.current_val;
		}
		inline void set(const MmerSignature& x)
		{
			str = x.str;
			current_val = x.current_val;
		}

		void set_from_binary(uint32_t value)
		{
			str = value;
			current_val = norm[str];
		}

		inline void insert(const char* seq)
		{
			switch (len)
			{
			case 4:
				str = (seq[0] << 6) + (seq[1] << 4) + (seq[2] << 2) + seq[3];
				break;
			case 5:
				str = (seq[0] << 8) + (seq[1] << 6) + (seq[2] << 4) + (seq[3] << 2) + (seq[4]);
				break;
			case 6:
				str = (seq[0] << 10) + (seq[1] << 8) + (seq[2] << 6) + (seq[3] << 4) + (seq[4] << 2) + (seq[5]);
				break;
			case 7:
				str = (seq[0] << 12) + (seq[1] << 10) + (seq[2] << 8) + (seq[3] << 6) + (seq[4] << 4) + (seq[5] << 2) + (seq[6]);
				break;
			case 8:
				str = (seq[0] << 14) + (seq[1] << 12) + (seq[2] << 10) + (seq[3] << 8) + (seq[4] << 6) + (seq[5] << 4) + (seq[6] << 2) + (seq[7]);
				break;
			case 9:
				str = (seq[0] << 16) + (seq[1] << 14) + (seq[2] << 12) + (seq[3] << 10) + (seq[4] << 8) + (seq[5] << 6) + (seq[6] << 4) + (seq[7] << 2) + (seq[8]);
				break;
			case 10:
				str = (seq[0] << 18) + (seq[1] << 16) + (seq[2] << 14) + (seq[3] << 12) + (seq[4] << 10) + (seq[5] << 8) + (seq[6] << 6) + (seq[7] << 4) + (seq[8] << 2) + (seq[9]);
				break;
			case 11:
				str = (seq[0] << 20) + (seq[1] << 18) + (seq[2] << 16) + (seq[3] << 14) + (seq[4] << 12) + (seq[5] << 10) + (seq[6] << 8) + (seq[7] << 6) + (seq[8] << 4) + (seq[9] << 2) + (seq[10]);
				break;
			default:
				break;
			}

			current_val = norm[str];
		}
	};

	template<typename HASH_T>
	class MmerMinHash
	{
		uint32_t len;
		uint64_t str;
		uint64_t rev;
		uint64_t mask;
		uint64_t current_hash;

		uint32_t rc_shift;

		static uint32_t get_rev_compl_shift(uint32_t len)
		{
			return (len + 31) / 32 * 32 - len;
		}

		//instead of x >> 2*p there is (x>>p) >> p, because
		//for p=32 we would have x >> 64 which is UB
		static uint64_t shr_2p(uint64_t x, uint64_t p)
		{
			return (x >> p) >> p;
		}

		static uint64_t bswap_uint64(uint64_t val)
		{
#ifdef _MSC_VER
			return _byteswap_uint64(val);
#elif defined(__GNUC__) || defined(__clang__)
			return __builtin_bswap64(val);
#else //unknown. Use the fastest "standard" way I've found
			val = ((val << 8) & 0xFF00FF00FF00FF00ULL) + ((val >> 8) & 0x00FF00FF00FF00FFULL);
			val = ((val << 16) & 0xFFFF0000FFFF0000ULL) + ((val >> 16) & 0x0000FFFF0000FFFFULL);
			val = (val << 32) + (val >> 32);
			return val;
#endif
		}

		static uint64_t get_rev_compl(uint64_t str, uint32_t shift)
		{
			const uint64_t mask4 = 0b0000111100001111000011110000111100001111000011110000111100001111ull;
			const uint64_t mask2 = 0b0011001100110011001100110011001100110011001100110011001100110011ull;

			uint64_t sf4 = ((str & mask4) << 4) + ((str >> 4) & mask4);
			uint64_t sf2 = ((sf4 & mask2) << 2) + ((sf4 >> 2) & mask2);

			return shr_2p((~bswap_uint64(sf2)), shift);

		}
	public:
		explicit MmerMinHash(uint32_t len) :
			len(len),
			mask((1ull << len * 2) - 1),
			rc_shift(get_rev_compl_shift(len))
		{
			clear();
		}

		inline void insert(uint8_t symb)
		{
			str <<= 2;
			str += symb;
			str &= mask;
			rev >>= 2;
			rev += (3 - (uint64_t)symb) << (len * 2 - 2);

			//current_hash = MurMur64Hash{}(MIN(str, rev)) & mask; //mkokot_TODO: should I "& mask" ?
			//current_hash = MurMur64Hash{}(MIN(str, rev)); //mkokot_TODO: should I "& mask" ?
			current_hash = HASH_T{}((std::min)(str, rev));
		}
		inline uint64_t get() const {
			return current_hash;
		}
		inline bool operator==(const MmerMinHash& x) const {
			return current_hash < x.current_hash;
		}
		inline bool operator<(const MmerMinHash& x) const {
			return current_hash < x.current_hash;
		}
		inline void clear() {
			str = rev = 0;
			current_hash = 0;
		}
		inline bool operator<=(const MmerMinHash& x) const {
			return current_hash <= x.current_hash;
		}
		inline void set(const MmerMinHash& x) {
			str = x.str;
			rev = x.rev;
			current_hash = x.current_hash;
		}

		void set_from_binary(uint64_t value)
		{
			str = value;
			rev = get_rev_compl(str, rc_shift);
			current_hash = HASH_T{}((std::min)(str, rev));
		}

		void set_from_binary(uint64_t new_str, uint64_t new_rev)
		{
			str = new_str;
			rev = new_rev;
			current_hash = HASH_T{}((std::min)(str, rev));
		}

		inline void insert(const char* seq) {
			//mkokot_TODO: rewrite more efficiently...
			for (uint32_t i = 0; i < len; ++i)
				insert(seq[i]);
		}
	};
}

#endif // ! MMER_H_
