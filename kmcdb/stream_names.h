#ifndef STREAM_NAMES_H_
#define STREAM_NAMES_H_
#include <string>

namespace kmcdb::stream_names
{
	const static std::string METADATA = "metadata";

	const static std::string HISTORY = "history";

	const static std::string SAMPLE_NAMES = "samples_names";

	static constexpr uint64_t BIN_ID_STR_LEN = 5;

	inline std::string BinMetadata(uint64_t bin_id)
	{
		return "bin_metadata_" + detail::to_string_leading_zeros(bin_id, BIN_ID_STR_LEN);
	}

	inline std::string Bin(uint64_t bin_id)
	{
		return "bin_" + detail::to_string_leading_zeros(bin_id, BIN_ID_STR_LEN);
	}

	inline std::string BinSufData(uint64_t bin_id)
	{
		return Bin(bin_id) + "_suf+data";
	}

	inline std::string BinLut(uint64_t bin_id)
	{
		return Bin(bin_id) + "_lut";
	}
}

#endif // ! STREAM_NAMES_H_
