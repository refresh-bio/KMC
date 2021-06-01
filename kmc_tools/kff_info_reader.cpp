#include "stdafx.h"
#include "kff_info_reader.h"
#include "defs.h"
#include <stdexcept>
#include <vector>
#include <array>
//#include <iostream>
#include <algorithm>
#include <cstring>

std::string CKFFInfoReader::ReadVarName()
{
	std::string res;
	while (true)
	{
		int c = fgetc(file);
		if (c == EOF)
			throw std::runtime_error("unexpected EOF");
		char _c = c;
		if (_c == 0)
			return res;
		res.push_back(_c);
	}
}

CKFFInfoReader::CKFFInfoReader(const std::string& path)
{
	file = my_fopen(path.c_str(), "rb");
	if (!file)
		throw std::runtime_error("Error: cannot open file " + path);

	// Check markers
	char marker[4];
	marker[3] = '\0';
	fread(marker, 1, 3, file);
	if (strncmp(marker, "KFF", 3) != 0)
		throw std::runtime_error("Error: missing KFF marker at the begining of file " + path);

	my_fseek(file, -3, SEEK_END);
	fread(marker, 1, 3, file);
	if (strncmp(marker, "KFF", 3) != 0)
		throw std::runtime_error("Error: missing KFF marker at the end of file " + path);

	my_fseek(file, -23, SEEK_END);
	char footer_size_str[12];
	fread(footer_size_str, 1, 12, file);
	bool footer_present = strcmp(footer_size_str, "footer_size") == 0;
	uint64_t first_index = std::numeric_limits<uint64_t>::max();
	std::array<uint8_t, 8> tmp;
	uint64_t nb_vars;
	char t;
	if (footer_present)
	{		
		uint64_t footer_size;
		fread(tmp.data(), 1, 8, file);
		LoadBigEndian(tmp.data(), footer_size);

		//std::cerr << "footer_size: " << footer_size << "\n";

		my_fseek(file, -((int64_t)footer_size + 3), SEEK_END);
		
		fread(tmp.data(), 1, 1, file);
		LoadBigEndian(tmp.data(), t);

		if (t != 'v')
			throw std::runtime_error("Error: footer should start as 'v' section, file " + path);

		fread(tmp.data(), 1, 8, file);
		LoadBigEndian(tmp.data(), nb_vars);

		//std::cerr << "footer nb_vars: " << nb_vars << "\n";

		for (uint64_t i = 0; i < nb_vars; ++i)
		{
			auto name = ReadVarName();
			uint64_t val;
			fread(tmp.data(), 1, 8, file);
			LoadBigEndian(tmp.data(), val);
			//std::cerr << name << ": " << val << "\n";

			kff_file_struct.footer[name] = val;

			if (name == "first_index")
				first_index = val;
		}
	}

	my_fseek(file, 3, SEEK_SET);
	uint8_t ver_minor, ver_major;
	fread(tmp.data(), 1, 1, file);
	LoadBigEndian(tmp.data(), ver_major);

	fread(tmp.data(), 1, 1, file);
	LoadBigEndian(tmp.data(), ver_minor);


	fread(tmp.data(), 1, 1, file);
	LoadBigEndian(tmp.data(), kff_file_struct.encoding);

	fread(tmp.data(), 1, 1, file);
	LoadBigEndian(tmp.data(), kff_file_struct.all_unique);
	if (kff_file_struct.all_unique == 0)
		throw std::runtime_error("Error: only unique k-mers in KFF file are supported, file " + path);

	uint8_t canonical;
	fread(tmp.data(), 1, 1, file);
	LoadBigEndian(tmp.data(), canonical);

	kff_file_struct.both_strands = canonical;

	uint32_t free_size;
	fread(tmp.data(), 1, sizeof(free_size), file);
	LoadBigEndian(tmp.data(), free_size);

	my_fseek(file, free_size, SEEK_CUR); //skip free block

	fread(tmp.data(), 1, 1, file);
	LoadBigEndian(tmp.data(), t);
	if (t == 'i')
	{
		uint64_t p = my_ftell(file);
		if (first_index != std::numeric_limits<uint64_t>::max() && first_index != p) //first index was defined in footer and is different the a real first index
			throw std::runtime_error("Error: footer defines 'first_index' but there is also an index as first section and the positions are inconsistent, file " + path);
		first_index = p;
	}

	if (first_index == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: no first_index in the footer and first section is not an index, file: " + path);


	while (first_index)
	{
		//std::cerr << "reading single index\n";

		my_fseek(file, first_index, SEEK_SET);

		fread(tmp.data(), 1, 1, file);
		LoadBigEndian(tmp.data(), t);
		if (t != 'i')
			throw std::runtime_error("Error: missing index");

		fread(tmp.data(), 1, 8, file);

		LoadBigEndian(tmp.data(), nb_vars);
		//std::cerr << "index nb_vars: " << nb_vars << "\n";
		//std::cerr << "Index: \n";

		int64_t this_index_end = my_ftell(file) + nb_vars * (sizeof(uint64_t) + 1) + sizeof(uint64_t);
		for (uint64_t i = 0; i < nb_vars; ++i)
		{
			fread(tmp.data(), 1, 1, file);
			LoadBigEndian(tmp.data(), t);
			int64_t rel_pos;

			fread(tmp.data(), 1, 8, file);
			LoadBigEndian(tmp.data(), rel_pos);

			index.emplace_back(t, this_index_end + rel_pos);
			//std::cerr << t << "\t" << rel_pos << "\n";
		}

		fread(tmp.data(), 1, 8, file);
		LoadBigEndian(tmp.data(), first_index);
	}

	std::sort(index.begin(), index.end(), [](const auto& lhs, const auto& rhs) {return lhs.section_pos < rhs.section_pos; });

	//std::cerr << "Final index:\n";
	for (const auto& e : index)
	{
		//std::cerr << e.first << "\t" << e.second << "\n";

		my_fseek(file, e.section_pos, SEEK_SET);
		fread(tmp.data(), 1, 1, file);
		LoadBigEndian(tmp.data(), t);
		if (t != e.section_type)
			throw std::runtime_error("Error: KFF index is inconsistent with file content");
	}

	for (auto e : index)
	{
		my_fseek(file, e.section_pos, SEEK_SET);
		fread(tmp.data(), 1, 1, file);
		LoadBigEndian(tmp.data(), t);
		if (t == 'i')
			continue; //skip index sections
		else if (t == 'v') //variable section
			ReadVariableSection();
		else if (t == 'r') //raw section
			ReadRawSection();
		else if (t == 'm') //minimizer section
			ReadMinimizerSection();
		else
			throw std::runtime_error(std::string("Error: unsupported section type (") + t + "), file " + path);
	}

	if (kff_file_struct.scopes.size())
	{
		if (kff_file_struct.scopes.back().data_sections.empty()) // remove prev variable section if empty
			kff_file_struct.scopes.pop_back();
	}
}

void CKFFInfoReader::ReadVariableSection()
{
	uint64_t nb_vars;	
	std::array<uint8_t, 8> tmp;
	fread(tmp.data(), 1, sizeof(uint64_t), file);
	LoadBigEndian(tmp.data(), nb_vars);
	CKFFVariables section;
	for (uint64_t i = 0; i < nb_vars; ++i)
	{
		auto var_name = ReadVarName();
		uint64_t val;
		fread(tmp.data(), 1, sizeof(uint64_t), file);
		LoadBigEndian(tmp.data(), val);
		if (var_name == "k")
			section.kmer_size = val;
		else if (var_name == "max")
			section.max_in_block = val;
		else if (var_name == "data_size")
			section.data_size = val;
		else if (var_name == "m")
			section.minimizer_size = val;
		else if (var_name == "ordered")
			section.ordered = val;
	}
	if (kff_file_struct.scopes.size())
	{
		if (kff_file_struct.scopes.back().data_sections.empty()) // remove prev variable section if empty
			kff_file_struct.scopes.pop_back();
	}
	kff_file_struct.scopes.push_back(section);
}

void CKFFInfoReader::ReadRawSection()
{
	if (kff_file_struct.scopes.empty())
		throw std::runtime_error("Error: raw section declared without variable section");
	auto& scope = kff_file_struct.scopes.back();

	if(scope.kmer_size == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: `k` variable was not defined for raw section");

	if(scope.max_in_block == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: `max` variable was not defined for raw section");

	if(scope.data_size == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: `data_size` variable was not defined for raw section");

	uint64_t nb_blocks;
	std::array<uint8_t, 8> tmp;

	fread(tmp.data(), 1, 8, file);
	LoadBigEndian(tmp.data(), nb_blocks);

	CKFFDataSection data_section;
	data_section.type = KFFDataSectionType::RAW;
	data_section.nb_blocks = nb_blocks;
	data_section.data_start_pos = my_ftell(file);
	scope.data_sections.push_back(data_section);
}

void CKFFInfoReader::ReadMinimizerSection()
{
	if (kff_file_struct.scopes.empty())
		throw std::runtime_error("Error: minimizer section declared without variable section");
	auto& scope = kff_file_struct.scopes.back();

	if (scope.kmer_size == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: `k` variable was not defined for minimizer section");

	if (scope.max_in_block == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: `max` variable was not defined for minimizer section");

	if (scope.data_size == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: `data_size` variable was not defined for minimizer section");

	if (scope.minimizer_size == std::numeric_limits<uint64_t>::max())
		throw std::runtime_error("Error: `m` variable was not defined for minimizer section");

	uint64_t nb_blocks;
	std::array<uint8_t, 8> tmp;

	uint64_t minimizer_bytes = (scope.minimizer_size + 3) / 4;

	CKFFDataSection data_section;
	data_section.minimizer.resize(minimizer_bytes);
	fread(data_section.minimizer.data(), 1, minimizer_bytes, file);

	fread(tmp.data(), 1, 8, file);
	LoadBigEndian(tmp.data(), nb_blocks);

	data_section.type = KFFDataSectionType::RAW;
	data_section.nb_blocks = nb_blocks;
	data_section.data_start_pos = my_ftell(file);
	scope.data_sections.push_back(data_section);
}

CKFFInfoReader::~CKFFInfoReader()
{
	fclose(file);
}