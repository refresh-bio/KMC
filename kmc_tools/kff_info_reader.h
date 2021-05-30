#pragma once
#include <string>
#include <vector>
#include <limits>

/*
TODO KFF: add checking if section is ordered!!!!
*/

template<typename T>
void LoadBigEndian(const uint8_t* buff, T& data)
{
	data = T{};
	for (int32_t b = sizeof(data) - 1; b >= 0; --b)
	{
		uint8_t byte = *buff++;
		data += (T)byte << (8 * b);
	}
}

struct CKFFIndexPair
{
	char section_type;
	uint64_t section_pos;

	CKFFIndexPair(char section_type, uint64_t section_pos) : section_type(section_type), section_pos(section_pos) { }
};

enum class KFFDataSectionType { RAW, MINIMIZER };

struct CKFFDataSection
{
	KFFDataSectionType type;
	uint64_t nb_blocks;
	std::vector<uint8_t> minimizer;
	uint64_t data_start_pos;
};

struct CKFFVariables
{
	uint64_t kmer_size = std::numeric_limits<uint64_t>::max();
	uint64_t data_size = std::numeric_limits<uint64_t>::max(); //counter size
	uint64_t minimizer_size = std::numeric_limits<uint64_t>::max();
	uint64_t max_in_block = std::numeric_limits<uint64_t>::max();
	std::vector<CKFFDataSection> data_sections;
};

struct CKFFFileStruct
{	
	bool both_strands;	
	std::vector<CKFFVariables> scopes;
};

class CKFFInfoReader
{
	FILE* file = nullptr;

	std::string ReadVarName();

	std::vector<CKFFIndexPair> index;

	CKFFFileStruct kff_file_struct;

	void ReadVariableSection();

	void ReadRawSection();

	void ReadMinimizerSection();

public:
	explicit CKFFInfoReader(const std::string& path);

	std::vector<CKFFIndexPair> GetIndex() const
	{
		return index;
	}

	CKFFFileStruct GetKffFileStruct() const
	{
		return kff_file_struct;
	}

	~CKFFInfoReader();
};
