#pragma once

#include <string>
#include <vector>
#include <cstdint>

template<typename T>
void StoreBigEndian(uint8_t* buff, const T& data)
{
	using unsigned_T = typename std::make_unsigned<T>::type;
	unsigned_T _data = static_cast<unsigned_T>(data);

	for (int32_t b = sizeof(_data) - 1; b >= 0; --b)
		*buff++ = _data >> (8 * b);
}

class CKFFWriter
{
	static constexpr uint8_t VER_MAJOR = 1;
	static constexpr uint8_t VER_MINOR = 0;
	FILE* file;

	uint64_t k;

	uint64_t counter_size;
	uint64_t min_count;
	uint64_t max_count;

	std::vector<int64_t> index;
	uint64_t cur_file_size = 0;

	void storeIndexPair(const char* type, int64_t val, std::vector<uint8_t>& tmp);



	struct SectionPartState
	{
		uint64_t where_to_store_nb_recs;
		uint64_t nb_recs;
	};
	SectionPartState section_part_state;

public:
	explicit CKFFWriter(const std::string& path, uint8_t canonical, uint64_t k, uint64_t counter_size, uint64_t min_count, uint64_t max_count, uint8_t encoding = 0b00011011);

	void StoreWholeSection(uint8_t* data, uint64_t n_kmers);

	void InitSection();
	void StoreSectionPart(uint8_t* data, uint64_t n_recs);
	void FinishSection();

	~CKFFWriter();
};