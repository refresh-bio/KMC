/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.1
  Date   : 2022-01-04
*/
#ifndef _SIG_TO_BIN_MAP_H
#define _SIG_TO_BIN_MAP_H
#include <fstream>
#include <string>
#include <cinttypes>
#include <vector>
#include <cassert>

class CSigToBinMap {
	//file format is:
	//KMCM - start marker (4 bytes), meanign KMC Mapping
	//signature_len (8 bytes)
	//n_bins (8 bytes)
	// array of 4^sig_len + 1 elements, each 4 B
	//KMCM - end marker (4 bytes)
	uint64_t sig_len;
	uint64_t n_bins;
	std::vector<int32_t> mapping;
	void serialize(std::ofstream& out, uint64_t num) const {
		char c[8];
		c[0] = num & 255;
		c[1] = (num >> 8) & 255;
		c[2] = (num >> 16) & 255;
		c[3] = (num >> 24) & 255;
		c[4] = (num >> 32) & 255;
		c[5] = (num >> 40) & 255;
		c[6] = (num >> 48) & 255;
		c[7] = num >> 56;
		out.write(c, 8);
	};
	void serialize32(std::ofstream& out, uint32_t num) const {
		char c[4];
		c[0] = num & 255;
		c[1] = (num >> 8) & 255;
		c[2] = (num >> 16) & 255;
		c[3] = num >> 24;
		out.write(c, 4);
	};

	void load(std::ifstream& in, uint64_t& num) {
		unsigned char c[8];
		in.read((char*)c, 8);
		num = c[0];
		num += (size_t)c[1] << 8;
		num += (size_t)c[2] << 16;
		num += (size_t)c[3] << 24;
		num += (size_t)c[4] << 32;
		num += (size_t)c[5] << 40;
		num += (size_t)c[6] << 48;
		num += (size_t)c[7] << 56;
	};

	void load32(std::ifstream& in, uint32_t& num) {
		unsigned char c[4];
		in.read((char*)c, 4);
		num = c[0];
		num += (size_t)c[1] << 8;
		num += (size_t)c[2] << 16;
		num += (size_t)c[3] << 24;
	};

	void write_marker(std::ofstream& out) const {
		out.write("KMCM", 4);
	}
	void check_marker(std::ifstream& in) const {
		char marker[4];
		in.read(marker, 4);
		if(strncmp(marker, "KMCM", 4))
			throw std::runtime_error("KMCM marker was expected");
	}
	size_t get_mapping_size() const {
		return (1ull << (sig_len * 2)) + 1;
	}
	
public:
	inline static std::string kmer_to_string(uint64_t kmer, uint8_t len) {
		std::string str_kmer;
		for (int i = len - 1; i >= 0; --i) {
			auto s = (kmer >> (2 * i)) & 3;
			str_kmer.push_back("ACGT"[s]);
		}
		return str_kmer;
	}
	uint32_t GetSigLen() const {
		return sig_len;
	}
	
	uint32_t GetNBins() const {
		return n_bins;
	}
	const int32_t* GetMapping() const {
		return mapping.data();
	}

	CSigToBinMap(const CSigToBinMap&) = delete;
	CSigToBinMap& operator=(const CSigToBinMap&) = delete;
	CSigToBinMap(const std::string& path) {
		std::ifstream in(path, std::ios::binary);
		if (!in)
			throw std::runtime_error("Cannot open file " + path);
		check_marker(in);
		load(in, sig_len);
		load(in, n_bins);
		mapping.resize(get_mapping_size());
		for (auto& x : mapping) {
			uint32_t u;
			load32(in, u);
			x = static_cast<int32_t>(u);
		}			
		check_marker(in);
	}
	CSigToBinMap(size_t sig_len, size_t n_bins, int32_t* mapping) :
		sig_len(sig_len),
		n_bins(n_bins),
		mapping(get_mapping_size())
	{
		std::copy_n(mapping, this->mapping.size(), this->mapping.data());
	}
	void Serialize(const std::string& path) const {
		std::ofstream out(path, std::ios::binary);
		if (!out)
			throw std::runtime_error("Cannot open file " + path);
		write_marker(out);
		serialize(out, sig_len);
		serialize(out, n_bins);
		assert(mapping.size() == get_mapping_size());
		for (auto x : mapping)
			serialize32(out, x);
		write_marker(out);
	}
	//dump in textual form
	void Dump(std::ostream& out) const {
		out << "sig_len: " << sig_len << "\n";
		out << "n_bins: " << n_bins << "\n";
		for (size_t i = 0; i < mapping.size() - 1; ++i)
			out << i << "\t" << kmer_to_string(i, sig_len) << "\t" << mapping[i] << "\n";
		auto last_i = mapping.size() - 1;
		out << last_i << "\t" << kmer_to_string(last_i, sig_len + 1) << "\t" << mapping[last_i] << "\n";
	}
};

#endif