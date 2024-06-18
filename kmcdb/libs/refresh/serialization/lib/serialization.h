#ifndef _SERIALIZATION
#define _SERIALIZATION

#include <cassert>
#include <type_traits>
#include <ostream>
#include <istream>
#include <utility>

#include "serialization.h"
#if defined(__cpp_lib_endian) || defined(__cpp_lib_bit_cast)
#include <bit>
#else
#include <cstring>
#endif

namespace refresh
{
	const uint32_t REFRESH_BUILD_SERIALIZATION = 1;
	namespace serialization
	{
		namespace detail
		{
#ifdef __cpp_lib_endian
			using endian = std::endian;
#else
			enum class endian
			{
#if defined(_MSC_VER) && !defined(__clang__)
				little = 0,
				big = 1,
				native = little
#else
				little = __ORDER_LITTLE_ENDIAN__,
				big = __ORDER_BIG_ENDIAN__,
				native = __BYTE_ORDER__
#endif
		};
#endif

			static_assert(endian::native == endian::little || endian::native == endian::big, "refresh::serialization does not support mixed-endian");

#ifndef __cpp_lib_bit_cast
			template<class To, class From>
			std::enable_if_t<
				sizeof(To) == sizeof(From) &&
				std::is_trivially_copyable_v<From>&&
				std::is_trivially_copyable_v<To>,
				To>
				// constexpr support needs compiler magic
				constexpr bit_cast(const From& src) noexcept
			{
				static_assert(std::is_trivially_constructible_v<To>,
					"This implementation additionally requires "
					"destination type to be trivially constructible");

				To dst;
				std::memcpy(&dst, &src, sizeof(To));
				return dst;
			}
#endif

			template<size_t N_BYTES>
			struct UnsignedIntImpl {};

			template<>
			struct UnsignedIntImpl<1> {
				using type = uint8_t;
			};

			template<>
			struct UnsignedIntImpl<2> {
				using type = uint16_t;
			};

			template<>
			struct UnsignedIntImpl<4> {
				using type = uint32_t;
			};

			template<>
			struct UnsignedIntImpl<8> {
				using type = uint64_t;
			};

			template<size_t N_BYTES>
			using UnsignedInt = typename UnsignedIntImpl<N_BYTES>::type;

			//callback must be compatibile with void store_byte(uint8_t byte);
			template<typename T, typename STORE_BYTE_CALLBACK>
			void read_bytes_little_endian_impl(const T& value, const STORE_BYTE_CALLBACK& store_byte)
			{
				static_assert(std::is_unsigned<T>::value, "This only works for unsigned");
				for (uint32_t i = 0; i < sizeof(value); ++i)
				{
					uint8_t byte = static_cast<uint8_t>(value >> (8 * i));
					store_byte(byte);
				}
			}

			//callback must be compatibile with void store_byte(uint8_t byte);
			template<typename T, typename STORE_BYTE_CALLBACK>
			void read_bytes_little_endian_impl(const T& value, const STORE_BYTE_CALLBACK& store_byte, size_t only_num_bytes)
			{
				assert(only_num_bytes <= sizeof(value));
				static_assert(std::is_unsigned<T>::value, "This only works for unsigned");
				for (uint32_t i = 0; i < only_num_bytes; ++i)
				{
					uint8_t byte = static_cast<uint8_t>(value >> (8 * i));
					store_byte(byte);
				}
			}

			//callback must be compatibile with uint8_t read_byte();
			template<typename T, typename READ_BYTE_CALLBACK>
			void write_bytes_little_endian_impl(T& value, const READ_BYTE_CALLBACK& read_byte)
			{
				static_assert(std::is_unsigned<T>::value, "This only works for unsigned");
				value = T{};
				for (uint32_t i = 0; i < sizeof(value); ++i)
				{
					T byte = read_byte();
					value |= byte << (8 * i);
				}
			}

			//callback must be compatibile with uint8_t read_byte();
			template<typename T, typename READ_BYTE_CALLBACK>
			void write_bytes_little_endian_impl(T& value, const READ_BYTE_CALLBACK& read_byte, size_t serialized_size_bytes)
			{
				assert(serialized_size_bytes <= sizeof(T));
				static_assert(std::is_unsigned<T>::value, "This only works for unsigned");
				value = T{};
				for (uint32_t i = 0; i < serialized_size_bytes; ++i)
				{
					T byte = read_byte();
					value |= byte << (8 * i);
				}
			}

			template<typename UNSIGNED_T, typename BASE_T>
			UNSIGNED_T change_type(const BASE_T& val) {
#ifndef __cpp_lib_bit_cast
				return bit_cast<UNSIGNED_T, BASE_T>(val);
#else
				return std::bit_cast<UNSIGNED_T, BASE_T>(val);
#endif
			}

			//callback must be compatibile with void store_byte(uint8_t byte);
			template<typename T, typename STORE_BYTE_CALLBACK>
			void read_bytes_little_endian(const T& value, const STORE_BYTE_CALLBACK& store_byte)
			{
				if constexpr (std::is_unsigned_v<T>)
					detail::read_bytes_little_endian_impl(value, store_byte);
				else
				{
					using unsigned_type = detail::UnsignedInt<sizeof(T)>;
					auto converted_val = detail::change_type<unsigned_type>(value);
					detail::read_bytes_little_endian_impl(converted_val, store_byte);
				}
			}

			//callback must be compatibile with void store_byte(uint8_t byte);
			template<typename T, typename STORE_BYTE_CALLBACK>
			void read_bytes_little_endian(const T& value, const STORE_BYTE_CALLBACK& store_byte, size_t only_num_bytes)
			{
				if constexpr (std::is_unsigned_v<T>)
					detail::read_bytes_little_endian_impl(value, store_byte, only_num_bytes);
				else
				{
					using unsigned_type = detail::UnsignedInt<sizeof(T)>;
					auto converted_val = detail::change_type<unsigned_type>(value);
					detail::read_bytes_little_endian_impl(converted_val, store_byte, only_num_bytes);
				}
			}

			//callback must be compatibile with uint8_t read_byte();
			template<typename T, typename READ_BYTE_CALLBACK>
			void write_bytes_little_endian(T& value, const READ_BYTE_CALLBACK& read_byte)
			{
				if constexpr (std::is_unsigned_v<T>)
					detail::write_bytes_little_endian_impl(value, read_byte);
				else
				{
					using unsigned_type = detail::UnsignedInt<sizeof(T)>;
					unsigned_type res;
					detail::write_bytes_little_endian_impl(res, read_byte);
					value = detail::change_type<T>(res);
				}
			}

			//callback must be compatibile with uint8_t read_byte();
			template<typename T, typename READ_BYTE_CALLBACK>
			void write_bytes_little_endian(T& value, const READ_BYTE_CALLBACK& read_byte, size_t serialized_size_bytes)
			{
				if constexpr (std::is_unsigned_v<T>)
					detail::write_bytes_little_endian_impl(value, read_byte, serialized_size_bytes);
				else
				{
					using unsigned_type = detail::UnsignedInt<sizeof(T)>;
					unsigned_type res;
					detail::write_bytes_little_endian_impl(res, read_byte, serialized_size_bytes);
					value = detail::change_type<T>(res);
				}
			}
		} // namespace detail

		//single value
		template<typename T>
		void serialize_little_endian(const T& value, std::ostream& out)
		{
			detail::read_bytes_little_endian(value, [&out](uint8_t byte) {
				out.write(reinterpret_cast<const char*>(&byte), 1);
			});
		}

		template<typename T>
		void serialize_little_endian(const T& value, uint8_t*& out)
		{
			detail::read_bytes_little_endian(value, [&out](uint8_t byte) {
				*out++ = byte;
				});
		}

		template<typename T>
		void serialize_little_endian(const T& value, uint8_t*& out, size_t only_num_bytes)
		{
			assert(only_num_bytes <= sizeof(T));
			detail::read_bytes_little_endian(value, [&out](uint8_t byte) {
				*out++ = byte;
				}, only_num_bytes);
		}

		template<typename T>
		void serialize_little_endian(const T& value, std::vector<uint8_t>& out)
		{
			detail::read_bytes_little_endian(value, [&out](uint8_t byte) {
				out.push_back(byte);
			});
		}

		inline void serialize_string_no_len(const std::string& str, std::vector<uint8_t>& out) {
			out.insert(out.end(), str.begin(), str.end());
		}

		inline void serialize_string_with_len(const std::string& str, std::vector<uint8_t>& out) {
			size_t len = str.length();
			serialize_little_endian(len, out);
			serialize_string_no_len(str, out);
		}

		template<typename T>
		void load_little_endian(T& value, std::istream& in)
		{
			detail::write_bytes_little_endian(value, [&in]() {
				uint8_t byte;
				in.read(reinterpret_cast<char*>(&byte), 1);
				assert(in.gcount() == 1);
				return byte;
			});
		}

		template<typename T>
		void load_little_endian(T& value, const uint8_t* &in)
		{
			detail::write_bytes_little_endian(value, [&in]() {
				return *in++;
			});
		}

		template<typename T>
		void load_little_endian(T& value, const uint8_t*& in, size_t serialized_size_bytes)
		{
			detail::write_bytes_little_endian(value, [&in]() {
				return *in++;
				}, serialized_size_bytes);
		}

		template<typename T>
		void load_little_endian(T& value, const std::vector<uint8_t>& in, size_t& pos)
		{
			detail::write_bytes_little_endian(value, [&in, &pos]() {
				return in[pos++];
				});
		}

		template<typename T>
		void load_little_endian(T& value, const std::vector<uint8_t>& in, size_t& pos, size_t serialized_size_bytes)
		{
			detail::write_bytes_little_endian(value, [&in, &pos]() {
				return in[pos++];
				}, serialized_size_bytes);
		}

		inline void load_string(std::string& data, size_t len, const std::vector<uint8_t>& in, size_t& pos) {
			data.clear();
			data.insert(data.end(), in.begin() + pos, in.begin() + pos + len);
			pos += len;
		}

		inline void load_string(std::string& data, const std::vector<uint8_t>& in, size_t& pos) {
			size_t len;
			load_little_endian(len, in, pos);
			load_string(data, len, in, pos);
		}

		inline std::string load_string(const std::vector<uint8_t>& in, size_t& pos) {
			std::string res;
			load_string(res, in, pos);
			return res;
		}

		//std::pair
		template<typename T1, typename T2>
		void serialize_little_endian(const std::pair<T1, T2>& rec, std::vector<uint8_t>& out)
		{
			serialize_little_endian(rec.first, out);
			serialize_little_endian(rec.second, out);
		}

		template<typename T1, typename T2>
		void load_little_endian(std::pair<T1, T2>& rec, const std::vector<uint8_t>& in, size_t& pos)
		{
			load_little_endian(rec.first, in, pos);
			load_little_endian(rec.second, in, pos);
		}

		//std::tuple -> consider making this also more generic (variadic template?)
		template<typename T1, typename T2, typename T3>
		void serialize_little_endian(const std::tuple<T1, T2, T3>& rec, std::vector<uint8_t>& out)
		{
			serialize_little_endian(std::get<0>(rec), out);
			serialize_little_endian(std::get<1>(rec), out);
			serialize_little_endian(std::get<2>(rec), out);
		}

		template<typename T1, typename T2, typename T3>
		void load_little_endian(std::tuple<T1, T2, T3>& rec, const std::vector<uint8_t>& in, size_t& pos)
		{
			load_little_endian(std::get<0>(rec), in, pos);
			load_little_endian(std::get<1>(rec), in, pos);
			load_little_endian(std::get<2>(rec), in, pos);
		}

		//std::vector
		template<typename T>
		void serialize_little_endian(const std::vector<T>& recs, std::vector<uint8_t>& out)
		{
			uint64_t size = recs.size();
			serialize_little_endian(size, out);
			for (const auto& rec : recs)
				serialize_little_endian(rec, out);
		}

		template<typename T>
		void load_little_endian(std::vector<T>& recs, const std::vector<uint8_t>& in, size_t& pos)
		{
			uint64_t size;
			load_little_endian(size, in, pos);
			recs.resize(size);
			for (auto& rec : recs)
				load_little_endian(rec, in, pos);
		}

		template<typename T>
		void load_little_endian(std::vector<T>& recs, const uint8_t*& in)
		{
			uint64_t size;
			load_little_endian(size, in);
			recs.resize(size);
			for (auto& rec : recs)
				load_little_endian(rec, in);
		}

		inline uint64_t bswap_uint64(uint64_t val)
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

		inline uint32_t bswap_uint32(uint32_t val)
		{
#ifdef _MSC_VER
			return _byteswap_ulong(val);
#elif defined(__GNUC__) || defined(__clang__)
			return __builtin_bswap32(val);
#else //unknown. Use the fastest "standard" way I've found
			val = (val<<24) | ((val<<8) & 0x00ff0000) | ((val >> 8) & 0x0000ff00) | (val >> 24);
			return val;
#endif
		}
	} // namespace serialization
} // namespace refresh

#endif // _SERIALIZATION
