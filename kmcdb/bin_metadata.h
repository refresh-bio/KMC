#ifndef BIN_METADATA_H_
#define BIN_METADATA_H_
#include "utils.h"

namespace kmcdb
{
	//mkokot_TODO: this is not in the detail namespace, but the Metadata struct is, wired
	//mkokot_TODO: move to different file
	struct BinMetadata
	{
		uint64_t total_kmers{};

		void serialize(std::vector<uint8_t>& serialized) const
		{
			using namespace refresh::serialization;
			serialize_little_endian(total_kmers, serialized);
		}

		void load(const std::vector<uint8_t>& serialized)
		{
			using namespace refresh::serialization;
			size_t read_pos = 0;
			load_little_endian(total_kmers, serialized, read_pos);
		}
	};
}
#endif // ! BIN_METADATA_H_
