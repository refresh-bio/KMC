#ifndef UTILS_H_
#define UTILS_H_

#include "kmer.h"
#include <array>
#include "libs/refresh/archive/lib/archive_input.h"
#include "libs/refresh/archive/lib/archive_output.h"
#include "libs/refresh/serialization/lib/serialization.h"
#include <sstream>
namespace kmcdb
{
	struct Version
	{
		uint64_t major;
		uint64_t minor;
		uint64_t patch;
		auto operator<=>(const Version& rhs) const
		{
			return std::make_tuple(major, minor, patch) <=> std::make_tuple(rhs.major, rhs.minor, rhs.patch);
		}

		std::string to_string() const
		{
			std::ostringstream oss;
			oss << major << '.' << minor << '.' << patch;
			return oss.str();
		}

		bool from_string(const std::string& str_version)
		{
			std::istringstream iss(str_version);
			char dot1{}, dot2{};
			if (!(iss >> major >> dot1 >> minor >> dot2 >> patch))
				return false;

			if (dot1 != '.' || dot2 != '.')
				return false;

			return true;
		}
	};

	//software version
	//if the major version of kmcdb API is the same as the version of
	//kmcdb file it is compatible
	//mkokot_TODO: before public release we should change this to 4.0.0
	constexpr Version KMCDB_VERSION{ 0, 1, 0 };

	namespace detail
	{
		inline bool IsCompatible(const Version& ver)
		{
			//mkokot_TODO: may change in the future such that we will support older versions
			return KMCDB_VERSION.major == ver.major;
		}

		inline void serialize(const Version& version, std::vector<uint8_t>& serialized)
		{
			//serialize as string to be be able to easier find it in file using binary viewers
			using namespace refresh::serialization;
			serialize_string_with_len(version.to_string(), serialized);
		}
		inline void load(Version& version, const std::vector<uint8_t>& serialized, size_t& read_pos)
		{
			using namespace refresh::serialization;
			auto res = version.from_string(load_string(serialized, read_pos));
			assert(res); (void)res; //should throw exception if cannot load?
		}
	}

	//helpers to dispatch k-mer size
	template<unsigned SIZE>
	struct KmerSizeDispatcher
	{
		static void Dispatch(uint64_t kmer_len, const auto& callback)
		{
			auto min_k = 32 * (SIZE - 1);
			auto max_k = min_k + 32;
			if (kmer_len > min_k && kmer_len <= max_k)
				callback(std::integral_constant<unsigned, SIZE>{});
			else
				KmerSizeDispatcher<SIZE - 1>::Dispatch(kmer_len, callback);
		}
	};

	template<>
	struct KmerSizeDispatcher<0>
	{
		static void Dispatch(uint64_t /*kmer_len*/, const auto& /*callback*/)
		{
			throw std::runtime_error("k-mer size dispatcher failed!");
		}
	};

	template<unsigned MAX_KMER_LEN>
	void DispatchKmerSize(uint64_t kmer_length, const auto& callback)
	{
		if (kmer_length > MAX_KMER_LEN)
			throw std::runtime_error("k too large, use larger MAX_KMER_LEN");

		constexpr auto max_no_of_uint64_t_for_kmer = (MAX_KMER_LEN + 31) / 32;
		KmerSizeDispatcher<max_no_of_uint64_t_for_kmer>::Dispatch(kmer_length, callback);
	}

	//mkokot_TODO: move to details namespace
	using archive_output_t = refresh::archive_output;
	using archive_input_t = refresh::archive_input;

	namespace detail
	{
		template<typename T>
		struct is_tuple : std::false_type {};

		template<typename... Ts>
		struct is_tuple<std::tuple<Ts...>> : std::true_type {};

		template<typename T>
		inline constexpr bool is_tuple_v = is_tuple<T>::value;

		template<typename T>
		consteval size_t values_size()
		{
			if constexpr (is_tuple_v<T>)
				return std::tuple_size_v<T>;
			else
				return 1;
		}

		//not only for a tuple
		template<size_t Idx, typename T>
		auto& get(T& t)
		{
			if constexpr (is_tuple_v<T>)
				return std::get<Idx>(t);
			else {
				static_assert(Idx == 0);
				return t;
			}
		}

		// Implementations for tuple types
		template<typename Tuple, typename Func, size_t... Indices>
		void IterateImpl(Tuple&& tuple, Func&& func, std::index_sequence<Indices...>) {
			// Call the function with each tuple element and its index
			(func(std::integral_constant<size_t, Indices>{}, std::get<Indices>(tuple)), ...);
		}

		template<typename Tuple, typename Func>
		void IterateValues(Tuple&& tuple, Func&& func) {
			if constexpr (is_tuple_v<std::decay_t<Tuple>>) {
				constexpr size_t N = std::tuple_size_v<std::decay_t<Tuple>>;
				IterateImpl(std::forward<Tuple>(tuple), std::forward<Func>(func), std::make_index_sequence<N>{});
			}
			else {
				// Handle non-tuple types by treating them as a single-element tuple
				func(std::integral_constant<size_t, 0>{}, tuple);
			}
		}

		template<typename VALUE_T>
		void SetZeros(VALUE_T* values, size_t num_values)
		{
			for (uint64_t i = 0; i < num_values; ++i)
				IterateValues(values[i], [](auto /*idx*/, auto& val)
					{
						val = 0;
					});
		}

		template<typename VALUE_T>
		void LoadValues(VALUE_T* values, size_t num_values, const uint8_t* &ptr, const std::array<uint64_t, values_size<VALUE_T>()>& num_bytes_single_value)
		{
			for (uint64_t i = 0; i < num_values; ++i)
				IterateValues(values[i], [&]<typename T>(auto idx, T & val)
			{
				if constexpr (std::is_integral_v<T>)
					refresh::serialization::load_little_endian(val, ptr, num_bytes_single_value[idx]);
				else
					refresh::serialization::load_little_endian(val, ptr);
			});
		}

		template<typename VALUE_T>
		void LoadValues(VALUE_T* values, size_t num_values, std::vector<uint8_t>& data, size_t& read_pos, const std::array<uint64_t, values_size<VALUE_T>()>& num_bytes_single_value)
		{
			for (uint64_t i = 0; i < num_values; ++i)
				IterateValues(values[i], [&]<typename T>(auto idx, T & val)
			{
				if constexpr (std::is_integral_v<T>)
					refresh::serialization::load_little_endian(val, data, read_pos, num_bytes_single_value[idx]);
				else
					refresh::serialization::load_little_endian(val, data, read_pos);
			});
		}

		template<typename VALUE_T>
		void SerializeValues(const VALUE_T* values, size_t num_values, uint8_t* &ptr, const std::array<uint64_t, values_size<VALUE_T>()>& num_bytes_single_value)
		{
			for (size_t i = 0; i < num_values; ++i)
				detail::IterateValues(values[i], [&]<typename T>(auto idx, T & val)
			{
				if constexpr (std::is_integral_v<T>)
					refresh::serialization::serialize_little_endian(val, ptr, num_bytes_single_value[idx]);
				else
					refresh::serialization::serialize_little_endian(val, ptr);
			});
		}


		template<typename T, size_t Size>
		inline std::array<T, Size> vec_to_array(const std::vector<T>& vec)
		{
			assert(vec.size() == Size);
			std::array<T, Size> res;
			if constexpr (Size != 0)
				std::copy_n(vec.begin(), vec.size(), res.begin());
			return res;
		}

		inline std::string to_string_leading_zeros(uint64_t val, uint64_t len)
		{
			std::string res = std::to_string(val);
			while (res.length() < len)
				res = "0" + res;
			return res;
		}

		class serialization_buffer
		{
		public:
			serialization_buffer(const serialization_buffer&) = delete;
			serialization_buffer(serialization_buffer&&) = delete;
			serialization_buffer& operator=(serialization_buffer&&) = delete;
			serialization_buffer& operator=(const serialization_buffer&) = delete;

			serialization_buffer(size_t capacity) :
				capacity_(capacity),
				data_(new uint8_t[capacity]),
				cur_pos_(data_)
			{

			}

			size_t size() const
			{
				return cur_pos_ - data_;
			}

			size_t capacity() const
			{
				return capacity_;
			}

			void clear()
			{
				cur_pos_ = data_;
			}
			uint8_t*& get_end()
			{
				return cur_pos_;
			}
			const uint8_t* get()
			{
				return data_;
			}
			~serialization_buffer()
			{
				delete[] data_;
			}
		private:
			size_t capacity_{};
			uint8_t* data_{};
			uint8_t* cur_pos_{};
		};
	}
}

#endif // ! UTILS_H_
