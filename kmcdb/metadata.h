#ifndef METADATA_H_
#define METADATA_H_
#include "utils.h"

namespace kmcdb
{
	namespace detail
	{
		//mkokot_TODO: to be adjusted
		//maybe this should be exposed to the user, to be considered
		enum class KmersRepresentation
		{
			SortedPlain,
			SortedWithLUT //,
			//SortedCommonPrefixCompressed, //planned
			//Hashed // planned
		};

		//the idea is I will create classes
		//with naming convention
		//<Type><representation>
		//for example
		//ReaderSortedPlainForListing, BinReaderSortedPlainForListing
		//WriterSortedPlainForListing, BinWriterSortedPlainForListing
		//if there is further distinction if will be placed after, for example WriterSortedWithLUTRaw - for "raw" mode, i.e. separation to suffixes and lut done outsize (like in kmc)

		//for each representation I will have specific additional config
		//I will store the config as a std::variant

		inline std::string to_string(KmersRepresentation kmers_representation)
		{
			switch (kmers_representation)
			{
			case KmersRepresentation::SortedPlain: return "SortedPlain";
			case KmersRepresentation::SortedWithLUT: return "SortedWithLUT";
			//case KmersRepresentation::SortedCommonPrefixCompressed: return "SortedCommonPrefixCompressed"; // planned
			//case KmersRepresentation::Hashed: return "Hashed"; // planned
			}
			throw std::runtime_error("kmers_representation not covered by switch statement");
		}

		inline KmersRepresentation kmers_representation_from_string(const std::string& kmers_representation)
		{
			if (kmers_representation == "SortedPlain")
				return KmersRepresentation::SortedPlain;
			if (kmers_representation == "SortedWithLUT")
				return KmersRepresentation::SortedWithLUT;
			//if (kmers_representation == "SortedCommonPrefixCompressed") //planned
			//	return KmersRepresentation::SortedCommonPrefixCompressed;
			//if (kmers_representation == "Hashed") //planned
			//	return KmersRepresentation::Hashed;

			assert(false);
			throw std::runtime_error("Wrong kmers_representation");
		}

		enum class Value_Type : uint64_t { Double, Float, Uint64, Uint32, Uint16, Uint8 };

		inline std::string to_string(Value_Type value_type)
		{
			switch (value_type)
			{
			case Value_Type::Double: return "Double";
			case Value_Type::Float: return "Float";
			case Value_Type::Uint64: return "Uint64";
			case Value_Type::Uint32: return "Uint32";
			case Value_Type::Uint16: return "Uint16";
			case Value_Type::Uint8: return "Uint8";
			}
			throw std::runtime_error("value_type not covered by switch statement");
		}

		inline Value_Type value_type_from_string(const std::string& value_type)
		{
			if (value_type == "Double")
				return Value_Type::Double;
			if (value_type == "Float")
				return Value_Type::Float;
			if (value_type == "Uint64")
				return Value_Type::Uint64;
			if (value_type == "Uint32")
				return Value_Type::Uint32;
			if (value_type == "Uint16")
				return Value_Type::Uint16;
			if (value_type == "Uint8")
				return Value_Type::Uint8;

			assert(false);
			throw std::runtime_error("Wrong value_type");
		}

		template<typename VALUE_T>
		inline Value_Type ToValueType()
		{
			using PLAIN_VALUE_T = std::decay_t<VALUE_T>;
			if constexpr (std::is_same_v<PLAIN_VALUE_T, double>)
				return Value_Type::Double;
			else if constexpr (std::is_same_v<PLAIN_VALUE_T, float>)
				return Value_Type::Float;
			else if constexpr (std::is_same_v<PLAIN_VALUE_T, uint64_t>)
				return Value_Type::Uint64;
			else if constexpr (std::is_same_v<PLAIN_VALUE_T, uint32_t>)
				return Value_Type::Uint32;
			else if constexpr (std::is_same_v<PLAIN_VALUE_T, uint16_t>)
				return Value_Type::Uint16;
			else if constexpr (std::is_same_v<PLAIN_VALUE_T, uint8_t>)
				return Value_Type::Uint8;
			else
				static_assert(!sizeof(VALUE_T), "This type is not supported");
		}

		template<typename VALUE_T>
		inline std::vector<Value_Type> ToValueTypes()
		{
			if constexpr (detail::is_tuple_v<VALUE_T>)
			{
				std::vector<Value_Type> res(std::tuple_size_v<VALUE_T>);
				detail::IterateValues(VALUE_T{}, [&]<typename T>(auto idx, T& /*val*/)
				{
					res[idx] = detail::ToValueType<T>();
				});
				return res;
			}
			else
				return { detail::ToValueType<VALUE_T>() };
		}

		template<typename VALUE_T>
		inline bool AreSame(const std::vector<Value_Type>& value_types)
		{
			return value_types == ToValueTypes<VALUE_T>();
		}

		inline std::string to_string(const std::vector<Value_Type>& vec)
		{
			if (vec.size() == 1)
				return to_string(vec[0]);
			else
			{
				std::string res = "(";
				for (const auto& x : vec)
					res += to_string(x) + ", ";

				res.pop_back(); //remove ' '
				res.pop_back(); //remove ','
				res += ")";

				return res;
			}
		}
	}
	enum class SignatureToBinMapping
	{
		Modulo //just a value of signature (does not matter if it is KMC based value, or if its min-hash, we just somehow have a value, in this mode the bin id is this value % num_bins
		//ZigZag // planned, for example we are doing modulo 2 * n_bins and merge 0 with (2*n_bins-1), 1 with (2*n_bins -2), etc.
	};

	namespace detail
	{
		inline std::string to_string(SignatureToBinMapping signature_to_bin_mapping)
		{
			switch (signature_to_bin_mapping)
			{
			case SignatureToBinMapping::Modulo: return "Modulo";
			}
			throw std::runtime_error("signature_to_bin_mapping not covered by switch statement");
		}

		inline SignatureToBinMapping signature_to_bin_mapping_from_string(const std::string& signature_to_bin_mapping)
		{
			if (signature_to_bin_mapping == "Modulo")
				return SignatureToBinMapping::Modulo;

			assert(false);
			throw std::runtime_error("Wrong signature_to_bin_mapping");
		}
	}

	enum class SignatureSelectionScheme
	{
		MinHash
	};

	namespace detail
	{
		inline std::string to_string(SignatureSelectionScheme signature_selection_scheme)
		{
			switch (signature_selection_scheme)
			{
			case SignatureSelectionScheme::MinHash: return "MinHash";
			}
			throw std::runtime_error("signature_selection_scheme not covered by switch statement");
		}

		inline SignatureSelectionScheme signature_selection_scheme_from_string(const std::string& signature_selection_scheme)
		{
			if (signature_selection_scheme == "MinHash")
				return SignatureSelectionScheme::MinHash;

			assert(false);
			throw std::runtime_error("Wrong signature_selection_scheme: " + signature_selection_scheme);
		}
	}

	struct HistoryItem
	{
		//mkokot_TODO: consider adding hardware concurrency, some hardware info, like cpu, ram
		uint64_t open_time{}; //time when the writer object was created
		uint64_t close_time{}; //time when the writer object was closed
		uint64_t mem_peak_bytes{};
		std::string command_line;
		std::string info; //some additional info
		std::string system_info;

		//currently only capture std::cerr, and std::cout - capture other ways to stdout stderr seams to be tricky and I am not sure needed here
		std::string std_cout;
		std::string std_cerr;

	private:

		static void serialize_single(const std::string& name, uint64_t value, std::vector<uint8_t>& serialized)
		{
			using namespace refresh::serialization;
			serialize_string_with_len(name, serialized);
			serialize_little_endian(value, serialized);
		}
		static void serialize_single(const std::string& name, const std::string& value, std::vector<uint8_t>& serialized)
		{
			using namespace refresh::serialization;
			serialize_string_with_len(name, serialized);
			serialize_string_with_len(value, serialized);
		}

		static void load_single(uint64_t& value, const  std::vector<uint8_t>& serialized, size_t& read_pos)
		{
			using namespace refresh::serialization;
			load_little_endian(value, serialized, read_pos);
		}
		static void load_single(std::string& value, const std::vector<uint8_t>& serialized, size_t& read_pos)
		{
			using namespace refresh::serialization;
			load_string(value, serialized, read_pos);
		}

	public:
		void serialize(std::vector<uint8_t>& serialized) const
		{
#define SERIALIZE_WRAPPER(x) serialize_single(#x, x, serialized);
			SERIALIZE_WRAPPER(open_time)
			SERIALIZE_WRAPPER(close_time)
			SERIALIZE_WRAPPER(mem_peak_bytes)
			SERIALIZE_WRAPPER(command_line)
			SERIALIZE_WRAPPER(system_info)
			SERIALIZE_WRAPPER(info)
			SERIALIZE_WRAPPER(std_cout)
			SERIALIZE_WRAPPER(std_cerr)
#undef SERIALIZE_WRAPPER
		}

		void load(const std::vector<uint8_t>& serialized)
		{
			//unknown fields will be just skipped
			using namespace refresh::serialization;
			size_t read_pos{};
			std::string name;

#define LOAD_WRAPPER(x) if (#x == name) { load_single(x, serialized, read_pos); continue; }
			while (read_pos < serialized.size())
			{
				load_string(name, serialized, read_pos);

				LOAD_WRAPPER(open_time)
				LOAD_WRAPPER(close_time)
				LOAD_WRAPPER(mem_peak_bytes)
				LOAD_WRAPPER(command_line)
				LOAD_WRAPPER(system_info)
				LOAD_WRAPPER(info)
				LOAD_WRAPPER(std_cout)
				LOAD_WRAPPER(std_cerr)

				//if we are here we could signalize unknown field somehow
			}
#undef LOAD_WRAPPER
		}
	};

	struct Config
	{
		Version kmcdb_version = KMCDB_VERSION; // mkokot_TODO: hide this ?
		uint64_t kmer_len = 27;
		uint64_t num_samples = 1; //mkokot_TODO: consider rename to num_columns
		uint64_t num_bins = 512;
		uint64_t signature_len = 9; //mkokot_TODO: reconsider if this should be here, the problem is we don't use signatures for small k, but maybe setting zero is such case is just fine...
		std::vector<uint64_t> num_bytes_single_value;
		SignatureSelectionScheme signature_selection_scheme = SignatureSelectionScheme::MinHash;
		SignatureToBinMapping signature_to_bin_mapping = SignatureToBinMapping::Modulo;
	};

	struct ConfigSortedPlain
	{
		void serialize(std::vector<uint8_t>& /*serialized*/) const
		{
		}
		void load(const std::vector<uint8_t>& /*serialized*/, size_t& /*read_pos*/)
		{
		}
	};

	struct ConfigSortedWithLUT
	{
		uint64_t lut_prefix_len;
		void serialize(std::vector<uint8_t>& serialized) const
		{
			refresh::serialization::serialize_little_endian(lut_prefix_len, serialized);
		}
		void load(const std::vector<uint8_t>& serialized, size_t& read_pos)
		{
			refresh::serialization::load_little_endian(lut_prefix_len, serialized, read_pos);
		}
	};

	//planned
	//struct ConfigSortedCommonPrefixCompressed
	//{
	//	void serialize(std::vector<uint8_t>& /*serialized*/) const
	//	{
	//	}
	//	void load(const std::vector<uint8_t>& /*serialized*/, size_t& /*read_pos*/)
	//	{
	//	}
	//};

	//planned
	//struct ConfigHashed
	//{
	//	void serialize(std::vector<uint8_t>& /*serialized*/) const
	//	{
	//	}
	//	void load(const std::vector<uint8_t>& /*serialized*/, size_t& /*read_pos*/)
	//	{
	//	}
	//};

	namespace detail
	{
		using RepresentationConfigVariant = std::variant<
			ConfigSortedPlain,
			ConfigSortedWithLUT//,
			//ConfigSortedCommonPrefixCompressed, // planned
			//ConfigHashed // planned
		>;

		inline RepresentationConfigVariant to_variant(KmersRepresentation kmers_representation)
		{
			switch (kmers_representation)
			{
			case KmersRepresentation::SortedPlain: return ConfigSortedPlain{};
			case KmersRepresentation::SortedWithLUT: return ConfigSortedWithLUT{};
			//case KmersRepresentation::SortedCommonPrefixCompressed: return ConfigSortedCommonPrefixCompressed{}; //planned
			//case KmersRepresentation::Hashed: return ConfigHashed{}; //planned
			}
			throw std::runtime_error("kmers_representation not covered by switch statement");
		}

		struct Metadata
		{
			Config config;
			std::vector<Value_Type> value_types;
			KmersRepresentation kmers_representation = KmersRepresentation::SortedPlain;
			detail::RepresentationConfigVariant representation_config;

			void serialize(std::vector<uint8_t>& serialized) const
			{
				assert(value_types.size() == config.num_bytes_single_value.size());
				using namespace refresh::serialization;
				detail::serialize(config.kmcdb_version, serialized);
				serialize_little_endian(config.kmer_len, serialized);
				serialize_little_endian(config.num_samples, serialized);
				serialize_little_endian(config.num_bins, serialized);
				serialize_little_endian(config.signature_len, serialized);
				serialize_string_with_len(to_string(config.signature_selection_scheme), serialized);
				serialize_string_with_len(to_string(config.signature_to_bin_mapping), serialized);
				serialize_little_endian(value_types.size(), serialized);

				for (size_t i = 0; i < value_types.size(); ++i)
				{
					serialize_string_with_len(to_string(value_types[i]), serialized);
					serialize_little_endian(config.num_bytes_single_value[i], serialized);
				}

				serialize_string_with_len(to_string(kmers_representation), serialized);
				std::visit([&](auto&& arg)
					{
						arg.serialize(serialized);
					}, representation_config);
			}

			void load(const std::vector<uint8_t>& serialized)
			{
				using namespace refresh::serialization;
				size_t read_pos = 0;
				detail::load(config.kmcdb_version, serialized, read_pos);
				load_little_endian(config.kmer_len, serialized, read_pos);
				load_little_endian(config.num_samples, serialized, read_pos);
				load_little_endian(config.num_bins, serialized, read_pos);
				load_little_endian(config.signature_len, serialized, read_pos);
				config.signature_selection_scheme = signature_selection_scheme_from_string(
					load_string(serialized, read_pos));

				config.signature_to_bin_mapping = signature_to_bin_mapping_from_string(
					load_string(serialized, read_pos));

				size_t num_value_types;
				load_little_endian(num_value_types, serialized, read_pos);
				value_types.resize(num_value_types);
				config.num_bytes_single_value.resize(num_value_types);

				for (size_t i = 0; i < value_types.size(); ++i)
				{
					value_types[i] = value_type_from_string(
						load_string(serialized, read_pos));

					load_little_endian(config.num_bytes_single_value[i], serialized, read_pos);
				}

				kmers_representation = kmers_representation_from_string(
					load_string(serialized, read_pos));

				representation_config = detail::to_variant(kmers_representation);

				std::visit([&]<typename T>(T && arg)
				{
					arg.load(serialized, read_pos);
				}, representation_config);
			}
		};
	}

}

#endif // ! METADATA_H_
