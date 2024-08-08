#pragma once

#include <memory>
#include <algorithm>
#include <cassert>
#include <string>
#include <cstdint>
#include <cinttypes>
#include <cstring>
#include <type_traits>
#include <iterator>
#include <cmath>

#include "dragonbox.h"

namespace refresh
{
	namespace conversions
	{
		// ************************************************************************************
		class numeric_conversions
		{
			static const uint64_t e1 = 10ull;
			static const uint64_t e2 = 100ull;
			static const uint64_t e3 = 1000ull;
			static const uint64_t e4 = 10000ull;
			static const uint64_t e5 = 100000ull;
			static const uint64_t e6 = e3 * e3;
			static const uint64_t e7 = 10ull * e6;
			static const uint64_t e8 = 100ull * e6;
			static const uint64_t e9 = e3 * e6;
			static const uint64_t e10 = 10ull * e9;
			static const uint64_t e11 = 100ull * e9;
			static const uint64_t e12 = e6 * e6;
			static const uint64_t e13 = 10ull * e12;
			static const uint64_t e14 = 100ull * e12;
			static const uint64_t e15 = e3 * e12;
			static const uint64_t e16 = 10ull * e15;
			static const uint64_t e17 = 100ull * e15;
			static const uint64_t e18 = e9 * e9;
			static const uint64_t e19 = 10ull * e18;

			static char inline digits[4 * 10000];
			static uint64_t inline int_powers10[20] = { 1, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19 };
			struct _si {
				_si()
				{
					for (int i = 0; i < 10000; ++i)
					{
						int dig = i;

						digits[i * 4 + 3] = '0' + (dig % 10);
						dig /= 10;
						digits[i * 4 + 2] = '0' + (dig % 10);
						dig /= 10;
						digits[i * 4 + 1] = '0' + (dig % 10);
						dig /= 10;
						digits[i * 4 + 0] = '0' + static_cast<char>(dig);
					}
				}
			} static inline _init;
					
			// ************************************************************************************
			static inline int num_digits(uint64_t v)
			{
				if (v < e9)
				{
					if (v < e5)
						return (v < e1) ? 1 : (v < e2) ? 2 : (v < e3) ? 3 : (v < e4) ? 4 : 5;
					else
						return (v < e6) ? 6 : (v < e7) ? 7 : (v < e8) ? 8 : 9;
				}
				else
				{
					if (v < e14)
						return (v < e10) ? 10 : (v < e11) ? 11 : (v < e12) ? 12 : (v < e13) ? 13 : 14;
					else
						return (v < e15) ? 15 : (v < e16) ? 16 : (v < e17) ? 17 : (v < e18) ? 18 : (v < e19) ? 19 : 20;
				}
			}

			// ************************************************************************************
			static inline int num_digits_small(uint64_t v)
			{
				return (v < e1) ? 1 : (v < e2) ? 2 : (v < e3) ? 3 : 4;
			}

			// ************************************************************************************
			// Works only for len <= 4
			static inline void short_str_cpy_upto_4(char* dest, char* src, size_t len)
			{
//				assert(len <= 4);

				if (len == 1)
					dest[0] = src[0];
				else if (len == 2)
				{
					dest[0] = src[0];
					dest[1] = src[1];
				}
				else if (len == 3)
				{
					dest[0] = src[0];
					dest[1] = src[1];
					dest[2] = src[2];
				}
				else 
				{
					dest[0] = src[0];
					dest[1] = src[1];
					dest[2] = src[2];
					dest[3] = src[3];
				}
			}

			// ************************************************************************************
			static inline void short_str_cpy_4(char* dest, char* src)
			{
				dest[0] = src[0];
				dest[1] = src[1];
				dest[2] = src[2];
				dest[3] = src[3];
			}

			// ************************************************************************************
			static inline size_t append_exponent(int exponent, char* ptr)
			{
				*ptr++ = 'e';

				if (exponent < 0)
				{
					*ptr++ = '-';
					exponent = -exponent;
				}
				else
					*ptr++ = '+';

				int n_dig = (exponent < 100) ? 2 : (exponent < 1000) ? 3 : 4;

				short_str_cpy_upto_4(ptr, digits + exponent * 4 + (4 - n_dig), n_dig);

				return n_dig + 2;
			}

		public:
			// ************************************************************************************
			static inline size_t uint_to_pchar(uint64_t val, char* str)
			{
//				return jeaiii::to_text_from_integer(str, val) - str;

				if(val < e4)
				{
					size_t ndig = num_digits_small(val);

					short_str_cpy_upto_4(str, digits + val * 4 + (4 - ndig), ndig);

					return ndig;
				}
				else if (val < e8)
				{
					uint64_t dig1 = val / e4;
					uint64_t dig2 = val - dig1 * e4;

					int ndig = num_digits_small(dig1);

					short_str_cpy_upto_4(str, digits + dig1 * 4 + (4 - ndig), ndig);
					short_str_cpy_4(str + ndig, digits + dig2 * 4);

					return ndig + 4;
				}
				else if (val < e12)
				{
					uint64_t dig1 = val / e8;
					val -= dig1 * e8;
					uint64_t dig2 = val / e4;
					uint64_t dig3 = val - dig2 * e4;

					int ndig = num_digits_small(dig1);

					short_str_cpy_upto_4(str, digits + dig1 * 4 + (4 - ndig), ndig);
					short_str_cpy_4(str + ndig, digits + dig2 * 4);
					short_str_cpy_4(str + ndig + 4, digits + dig3 * 4);

					return ndig + 8;
				}
				else if (val < e16)
				{
					uint64_t dig1 = val / e12;
					val -= dig1 * e12;
					uint64_t dig2 = val / e8;
					val -= dig2 * e8;
					uint64_t dig3 = val / e4;
					uint64_t dig4 = val - dig3 * e4;

					int ndig = num_digits_small(dig1);

					short_str_cpy_upto_4(str, digits + dig1 * 4 + (4 - ndig), ndig);
					short_str_cpy_4(str + ndig, digits + dig2 * 4);
					short_str_cpy_4(str + ndig + 4, digits + dig3 * 4);
					short_str_cpy_4(str + ndig + 8, digits + dig4 * 4);

					return ndig + 12;
				}
				else
				{
					uint64_t dig1 = val / e16;
					val -= dig1 * e16;
					uint64_t dig2 = val / e12;
					val -= dig2 * e12;
					uint64_t dig3 = val / e8;
					val -= dig3 * e8;
					uint64_t dig4 = val / e4;
					uint64_t dig5 = val - dig4 * e4;

					int ndig = num_digits_small(dig1);

					short_str_cpy_upto_4(str, digits + dig1 * 4 + (4 - ndig), ndig);
					short_str_cpy_4(str + ndig, digits + dig2 * 4);
					short_str_cpy_4(str + ndig + 4, digits + dig3 * 4);
					short_str_cpy_4(str + ndig + 8, digits + dig4 * 4);
					short_str_cpy_4(str + ndig + 12, digits + dig5 * 4);

					return ndig + 16;
				}
			}

			// ************************************************************************************
			template<typename Floating>
			static size_t real_to_pchar(Floating val, char* out, size_t precision)
			{
				char* ptr = out;

				auto v = jkj::dragonbox::to_decimal(val);

				*ptr = '-';
				ptr += static_cast<int>(v.is_negative);

				uint64_t significand = static_cast<uint64_t>(v.significand);
				int exponent = static_cast<int>(v.exponent);

				int n_dig = num_digits(significand);

				if (static_cast<size_t>(n_dig) > precision)
				{
					significand += int_powers10[n_dig - precision] / 2;
					significand /= int_powers10[n_dig - precision];
					exponent += n_dig - static_cast<int>(precision);
					n_dig = static_cast<int>(precision);

					if (significand >= int_powers10[precision])
					{
						significand /= 10;
						exponent++;
					}
				}

				if (exponent == 0)
				{
					ptr += uint_to_pchar(significand, ptr);
				}
				else if (exponent > 0 || -exponent >= n_dig + 4)
				{
					if (n_dig == 1)
						*ptr++ = '0' + (char)significand;
					else
					{
						uint_to_pchar(significand, ptr + 1);
						ptr[0] = ptr[1];
						ptr[1] = '.';
						ptr += n_dig + 1;
						exponent += n_dig - 1;
					}

					ptr += append_exponent(exponent, ptr);
				}
				else if (-exponent < n_dig)
				{
					uint_to_pchar(significand, ptr);

					char* dot_ptr = ptr + n_dig + exponent;
					std::memmove(dot_ptr + 1, dot_ptr, -exponent);
					*dot_ptr = '.';
					ptr += n_dig + 1;
				}
				else// if (-exponent < n_dig + 4)
				{
					*ptr++ = '0';
					*ptr++ = '.';
					ptr[0] = '0';
					ptr[1] = '0';
					ptr[2] = '0';
					ptr[3] = '0';
					ptr[4] = '0';

					ptr += -exponent - n_dig;			// no. of zeros
					ptr += uint_to_pchar(significand, ptr);
				}

				return ptr - out;
			}
		};
	}

	// ************************************************************************************
	// Public functions
	// ************************************************************************************

	// ************************************************************************************
	// integral specialization
	template <typename Integer, typename std::enable_if<std::is_integral<Integer>::value, int>::type* = nullptr>
	size_t int_to_pchar(const Integer val, char* out, char term = '\0')
	{
		if constexpr (std::is_unsigned<Integer>::value)
		{
			auto r = conversions::numeric_conversions::uint_to_pchar((uint64_t)val, out);
			out[r] = term;
			return r + 1;
		}
		else
		{
			if (val < 0)
			{
				*out++ = '-';
				auto r = conversions::numeric_conversions::uint_to_pchar((uint64_t)-val, out);
				out[r] = term;

				return r + 2;
			}
			else
			{
				auto r = conversions::numeric_conversions::uint_to_pchar((uint64_t)val, out);
				out[r] = term;

				return r + 1;
			}
		}
	}

	// ************************************************************************************
	// floating point specialization
	template <typename Floating, typename std::enable_if<std::is_floating_point<Floating>::value, int>::type* = nullptr>
	size_t real_to_pchar(const Floating val, char* out, const size_t prec, char term = '\0')
	{
		if (val == 0) {
			*out++ = '0';
			*out = term;
			return 1 + 1;
		}
		if (std::isnan(val))
		{
			*out++ = 'n';
			*out++ = 'a';
			*out++ = 'n';
			*out = term;
			return 4;
		}
		if (std::isinf(val))
		{
			if (val < static_cast<Floating>(0)) {
				*out++ = '-';
				*out++ = 'i';
				*out++ = 'n';
				*out++ = 'f';
				*out = term;
				return 5;
			}
			else
			{
				*out++ = 'i';
				*out++ = 'n';
				*out++ = 'f';
				*out = term;
				return 4;
			}
		}

		size_t prec_bounded;

		if constexpr (sizeof(Floating) == sizeof(float))
			prec_bounded = std::clamp(prec, (size_t)1, (size_t)7);
		else if constexpr (sizeof(Floating) == sizeof(double))
			prec_bounded = std::clamp(prec, (size_t)1, (size_t)15);
		else if constexpr (sizeof(Floating) == sizeof(long double))
			prec_bounded = std::clamp(prec, (size_t)1, (size_t)19);

		auto r = conversions::numeric_conversions::real_to_pchar(val, out, prec_bounded);
		out[r] = term;

		return r + 1;
	}

	// ************************************************************************************
	template<typename Num1, typename Num2>
	size_t pair_to_pchar(const std::pair<Num1, Num2> val, char* out, const char sep, const size_t prec1, const size_t prec2, char term = '\0')
	{
		size_t n = 0;

		if constexpr (std::is_floating_point<decltype(val.first)>::value)
			n = real_to_pchar(val.first, out, prec1, sep);
		else if constexpr (std::is_integral<decltype(val.first)>::value)
			n = int_to_pchar(val.first, out, sep);
		else
			assert(0);

		if constexpr (std::is_floating_point<decltype(val.second)>::value)
			n += real_to_pchar(val.second, out + n, prec2, term);
		else if constexpr (std::is_integral<decltype(val.second)>::value)
			n += int_to_pchar(val.second, out + n, term);
		else
			assert(0);

		return n;
	}

	// ************************************************************************************
	template<typename Iterator>
	size_t sequence_to_pchar(Iterator first, Iterator last, char* out, const char sep, const size_t prec, char term = '\0')
	{
		char* ptr = out;

		if constexpr (std::is_floating_point<typename std::iterator_traits<Iterator>::value_type>::value)
		{
			for (auto p = first; p != last; ++p)
				ptr += real_to_pchar(*p, ptr, prec, sep);
		}
		else if constexpr (std::is_integral<typename std::iterator_traits<Iterator>::value_type>::value)
		{
			for (auto p = first; p != last; ++p)
				ptr += int_to_pchar(*p, ptr, sep);
		}
		else
			assert(0);

		if (first != last)		// remove final separator if any
			--ptr;

		*ptr++ = term;

		return ptr - out;
	}

	// ************************************************************************************
	template<typename Iterator>
	size_t pair_sequence_to_pchar(Iterator first, Iterator last, char* out, const char in_sep, const char out_sep, const size_t prec1, const size_t prec2, char term = '\0')
	{
		char* ptr = out;

		for (auto p = first; p != last; ++p)
			ptr += pair_to_pchar(*p, ptr, in_sep, prec1, prec2, out_sep);

		if (first != last)		// remove final separator if any
			--ptr;

		*ptr++ = term;

		return ptr - out;
	}

	// ************************************************************************************
	template<typename T>
	static inline size_t numeric_conversion_max_length()
	{
		if constexpr (std::is_integral<T>::value)
		{
			if constexpr (sizeof(T) == 1)
				return 3 + 1;
			else if constexpr (sizeof(T) == 2)
				return 5 + 1;
			else if constexpr (sizeof(T) == 4)
				return 10 + 1;
			else if constexpr (sizeof(T) == 8)
				return 20;
			else
				return 3 * sizeof(T);			// should be unnecessary
		}
		else if constexpr(std::is_floating_point<T>::value)
		{
			if constexpr (sizeof(T) == 4)
				return 7 + 4 + 1 + 1;				// value + exponent + dot + sign
			else if constexpr (sizeof(T) == 8)
				return 15 + 5 + 1 + 1;
			else if constexpr (sizeof(T) == 10 || sizeof(T) == 12 || sizeof(T) == 16)
				return 19 + 6 + 1 + 1;
			else
				return 0;
		}
		else
			return 0;			// unknown type
	}

	// ************************************************************************************
	template <typename Integer, typename std::enable_if<std::is_integral<Integer>::value, int>::type* = nullptr>
	std::string int_to_string(const Integer val)
	{
		std::string out;

		out.resize(numeric_conversion_max_length<Integer>());
		out.resize(int_to_pchar(val, out.data()));
		out.pop_back();										// remove terminator
				
		return out;
	}

	// ************************************************************************************
	template <typename Floating, typename std::enable_if<std::is_floating_point<Floating>::value, int>::type* = nullptr>
	std::string real_to_string(const Floating val, const size_t prec)
	{
		std::string out;

		out.resize(numeric_conversion_max_length<Floating>());
		out.resize(real_to_pchar(val, out.data(), prec));
		out.pop_back();										// remove terminator

		return out;
	}

	// ************************************************************************************
	template<typename Pair>
	std::string pair_to_string(const Pair val, const char sep, const size_t prec1, const size_t prec2)
	{
		std::string out;

		out.resize(numeric_conversion_max_length<decltype(val.first)>() + 1 + numeric_conversion_max_length<decltype(val.second)>());
		out.resize(pair_to_pchar(val, out.data(), sep, prec1, prec2));
		out.pop_back();										// remove terminator

		return out;
	}

	// ************************************************************************************
	template<typename Iterator>
	std::string sequence_to_string(Iterator first, Iterator last, const char sep, const size_t prec)
	{
		std::string out;
		std::string tmp;

		if constexpr (std::is_floating_point<typename std::iterator_traits<Iterator>::value_type>::value)
		{
			for (auto p = first; p != last; ++p)
			{
				out.append(real_to_string(*p, prec));
				out.push_back(sep);
			}
		}
		else if constexpr (std::is_integral<typename std::iterator_traits<Iterator>::value_type>::value)
		{
			for (auto p = first; p != last; ++p)
			{
				out.append(int_to_string(*p));
				out.push_back(sep);
			}
		}
		else
			assert(0);

		if (first != last)		// remove final separator if any
			out.pop_back();

		return out;
	}

	// ************************************************************************************
	template<typename Iterator>
	std::string pair_sequence_to_string(Iterator first, Iterator last, const char in_sep, const char out_sep, const size_t prec1, const size_t prec2)
	{
		std::string out;
		std::string tmp;

		for (auto p = first; p != last; ++p)
		{
			out.append(pair_to_string(*p, in_sep, prec1, prec2));
			out.push_back(out_sep);
		}

		if (first != last)		// remove final separator if any
			out.pop_back();

		return out;
	}
}
