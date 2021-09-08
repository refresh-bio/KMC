#include "kff_writer.h"
#include <stdexcept>
#include <cstring>

constexpr uint8_t CKFFWriter::VER_MAJOR;
constexpr uint8_t CKFFWriter::VER_MINOR;
CKFFWriter::CKFFWriter(const std::string& path, uint8_t canonical, uint64_t k, uint64_t counter_size, uint64_t min_count, uint64_t max_count, uint8_t encoding) :
	k(k), counter_size(counter_size),
	min_count(min_count),
	max_count(max_count)
{
	file = fopen(path.c_str(), "wb");
	if (!file)
		throw std::runtime_error("Cannot open file " + path);
	fwrite("KFF", 1, 3, file);
	cur_file_size += 3;

	fwrite(&VER_MAJOR, 1, 1, file);
	cur_file_size += 1;

	fwrite(&VER_MINOR, 1, 1, file);
	cur_file_size += 1;

	fwrite(&encoding, 1, 1, file);
	cur_file_size += 1;

	uint8_t unique_kmers = 1;
	fwrite(&unique_kmers, 1, 1, file);
	cur_file_size += 1;

	fwrite(&canonical, 1, 1, file);
	cur_file_size += 1;

	uint32_t free_size = 0;
	std::vector<uint8_t> tmp(sizeof(uint64_t));
	StoreBigEndian(tmp.data(), free_size);

	fwrite(tmp.data(), 1, sizeof(free_size), file);
	cur_file_size += sizeof(free_size);

	//variable section
	index.push_back(cur_file_size);

	fwrite("v", 1, 1, file);
	++cur_file_size;

	std::vector<std::pair<std::string, uint64_t>> pairs;
	pairs.emplace_back("k", k);
	pairs.emplace_back("max", 1);
	pairs.emplace_back("data_size", counter_size);
	pairs.emplace_back("ordered", 1);

	uint64_t nb_vars = pairs.size();

	StoreBigEndian(tmp.data(), nb_vars);
	fwrite(tmp.data(), 1, sizeof(nb_vars), file);
	cur_file_size += sizeof(nb_vars);

	for (const auto& elem : pairs)
	{
		fwrite(elem.first.c_str(), 1, elem.first.length() + 1, file);
		cur_file_size += elem.first.length() + 1;

		StoreBigEndian(tmp.data(), elem.second);
		fwrite(tmp.data(), 1, sizeof(elem.second), file);

		cur_file_size += sizeof(elem.second);
	}

}

void CKFFWriter::StoreWholeSection(uint8_t* data, uint64_t n_recs)
{
	index.push_back(cur_file_size);

	fwrite("r", 1, 1, file);
	++cur_file_size;

	std::vector<uint8_t> tmp(sizeof(uint64_t));

	StoreBigEndian(tmp.data(), n_recs);
	fwrite(tmp.data(), 1, sizeof(n_recs), file);
	cur_file_size += sizeof(n_recs);

	uint32_t kmer_bytes = (k + 3) / 4;
	uint32_t counter_bytes = counter_size;

	fwrite(data, 1, (kmer_bytes + counter_bytes) * n_recs, file);
	cur_file_size += (kmer_bytes + counter_bytes) * n_recs;
}

void CKFFWriter::InitSection()
{
	index.push_back(cur_file_size);

	fwrite("r", 1, 1, file);
	++cur_file_size;

	section_part_state.nb_recs = 0;
	section_part_state.where_to_store_nb_recs = cur_file_size;

	std::vector<uint8_t> tmp(sizeof(uint64_t));
	
	//prepare place for nb_recs in FinishSection
	fwrite(tmp.data(), 1, sizeof(uint64_t), file);

	cur_file_size += sizeof(uint64_t);
}

void CKFFWriter::StoreSectionPart(uint8_t* data, uint64_t n_recs)
{
	uint32_t kmer_bytes = (k + 3) / 4;
	uint32_t counter_bytes = counter_size;
	
	section_part_state.nb_recs += n_recs;

	fwrite(data, 1, (kmer_bytes + counter_bytes) * n_recs, file);
	cur_file_size += (kmer_bytes + counter_bytes) * n_recs;
}
void CKFFWriter::FinishSection()
{
	fseek(file, section_part_state.where_to_store_nb_recs, SEEK_SET);
	std::vector<uint8_t> tmp(sizeof(uint64_t));
	StoreBigEndian(tmp.data(), section_part_state.nb_recs);
	fwrite(tmp.data(), 1, sizeof(uint64_t), file);
	fseek(file, 0, SEEK_END);
}

void CKFFWriter::storeIndexPair(const char* type, int64_t val, std::vector<uint8_t>& tmp)
{
	fwrite(type, 1, 1, file);
	++cur_file_size;
	StoreBigEndian(tmp.data(), val);
	fwrite(tmp.data(), 1, sizeof(int64_t), file);
	cur_file_size += sizeof(int64_t);
}

CKFFWriter::~CKFFWriter()
{
	//Index
	uint64_t nb_sections = index.size() + 1; // +1 for footer 

	uint64_t index_size = 1; // 'i'
	index_size += 8; // nb_sections
	index_size += nb_sections * (1 + 8);
	index_size += 8; //next_index

	uint64_t index_start_pos = cur_file_size;
	uint64_t index_end_pos = cur_file_size + index_size;

	fwrite("i", 1, 1, file);
	++cur_file_size;

	std::vector<uint8_t> tmp(sizeof(uint64_t));
	StoreBigEndian(tmp.data(), nb_sections);
	fwrite(tmp.data(), 1, sizeof(uint64_t), file);
	cur_file_size += sizeof(uint64_t);

	storeIndexPair("v", index[0] - index_end_pos, tmp);
	for (uint64_t i = 1; i < index.size(); ++i)
		storeIndexPair("r", index[i] - index_end_pos, tmp);

	storeIndexPair("v", 0, tmp); //footer

	int64_t next_index = 0; //only one index
	StoreBigEndian(tmp.data(), next_index);
	fwrite(tmp.data(), 1, sizeof(int64_t), file);
	cur_file_size += sizeof(uint64_t);

	//Footer
	std::vector<std::pair<std::string, uint64_t>> footer;
	footer.emplace_back("first_index", index_start_pos);
	footer.emplace_back("min_count", min_count);
	footer.emplace_back("max_count", max_count);
	footer.emplace_back("counter_size", counter_size);

	uint64_t footer_size = 1; // 'v'
	footer_size += sizeof(uint64_t); // nb_vars
	for (const auto& e : footer)
		footer_size += e.first.length() + 1 + sizeof(e.second);
	footer_size += strlen("footer_size") + 1 + sizeof(uint64_t);

	footer.emplace_back("footer_size", footer_size);

	fwrite("v", 1, 1, file);
	uint64_t nb_vars = footer.size();
	StoreBigEndian(tmp.data(), nb_vars);
	fwrite(tmp.data(), 1, sizeof(uint64_t), file);
	for (const auto& e : footer)
	{
		fwrite(e.first.c_str(), 1, e.first.length() + 1, file);
		StoreBigEndian(tmp.data(), e.second);
		fwrite(tmp.data(), 1, sizeof(e.second), file);
	}

	fwrite("KFF", 1, 3, file);
	fclose(file);
}
