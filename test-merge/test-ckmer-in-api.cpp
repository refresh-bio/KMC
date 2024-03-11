#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <memory>
#include <numeric>
#include "kmc_file.h"

template<unsigned SIZE>
void kmer_to_str(CKmer<SIZE>& kmer, uint32_t kmer_len, std::string& out) {
	auto pos = 2 * kmer_len - 2;
	for(uint32_t i = 0 ; i < kmer_len ; ++i, pos -= 2)
		out[i] = "ACGT"[kmer.get_2bits(pos)];
}

class Runner {
	CKMCFile& kmc_file_ckmer_api;
	CKMCFile& kmc_file_ckmer;
public:
	Runner(CKMCFile& kmc_file_ckmer_api, CKMCFile& kmc_file_ckmer):
		kmc_file_ckmer_api(kmc_file_ckmer_api),
		kmc_file_ckmer(kmc_file_ckmer)
	{

	}
	template<unsigned SIZE>
	void Run() {

		CKmer<SIZE> kmer;
		uint32_t kmer_len = kmc_file_ckmer_api.KmerLength();
		CKmerAPI kmer_api(kmer_len);
		uint64_t count;
		uint64_t count_api;

		std::string str_kmer(kmer_len, ' ');
		std::string str_kmer_api(kmer_len, ' ');

		auto n_bins = kmc_file_ckmer.GetNBins();

		for(uint32_t bin_id = 0 ; bin_id < n_bins ; ++bin_id)
		{
			kmc_file_ckmer.StartBin(bin_id);
			kmc_file_ckmer_api.StartBin(bin_id);

			while(true) {
				bool is_kmer = kmc_file_ckmer.ReadNextKmerFromBin(kmer, count);
				bool is_kmer_api = kmc_file_ckmer_api.ReadNextKmerFromBin(kmer_api, count_api);

				if(is_kmer != is_kmer_api) {
					std::cerr << "is_kmer != is_kmer_api!\n";
					exit(1);
				}

				if(!is_kmer)
					break;

				if(count != count_api) {
					std::cerr << "count != count_api!\n";
					exit(1);
				}

				kmer_api.to_string(str_kmer_api);
				kmer_to_str(kmer, kmer_len, str_kmer);

				if(str_kmer_api != str_kmer) {
					std::cerr << "str_kmer_api != str_kmer!\n";
					std::cerr << "str_kmer_api: " << str_kmer_api << "\n";
					std::cerr << "str_kmer: " << str_kmer << "\n";
					//exit(1);
				}

				//std::cerr << str_kmer_api << " " << str_kmer << " " << count << " " << count_api << "\n";
			}
		}
	}
};

int main(int argc, char** argv)
{
	if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <kmc_db_path>\n";
        return 0;
    }
	std::string db_path = argv[1];

	CKMCFile kmc_file_ckmer_api;
	if(!kmc_file_ckmer_api.OpenForListingWithBinOrder(db_path)) {
		std::cerr << "Error: cannot open kmc database " << db_path << "\n";
		exit(1);
	}

	CKMCFile kmc_file_ckmer;
	if(!kmc_file_ckmer.OpenForListingWithBinOrder(db_path)) {
		std::cerr << "Error: cannot open kmc database " << db_path << "\n";
		exit(1);
	}
	
	Runner runner(kmc_file_ckmer_api, kmc_file_ckmer);
	
	DispatchKmerSize<>(kmc_file_ckmer.KmerLength(), runner);

}