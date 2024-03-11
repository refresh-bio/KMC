#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include "kmc_file.h"

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <outpath> <mapping_file> <input_1> [<input_2> ...]\n";
        return 0;
    }
    std::string outpath = argv[1];
    std::string mapping = argv[2];
    std::vector<std::string> inputs;
    for (int i = 3; i < argc; ++i)
        inputs.push_back(argv[i]);

    std::cerr << "mapping: " << mapping << "\n";
    std::cerr << "inputs:\n";
    for (auto& x : inputs)
        std::cerr << "\t" << x << "\n";

    std::vector<std::unique_ptr<CKMCFile>> kmc_dbs(inputs.size());

    for (int i = 0; i < inputs.size(); ++i)
    {
        kmc_dbs[i] = std::make_unique<CKMCFile>();
        if (!kmc_dbs[i]->OpenForListingWithBinOrder(inputs[i], mapping))
        {
            std::cerr << "Error: cannot open kmc database " << inputs[i] << "\n";
            return 1;
        }
    }

    struct KmerAndCount {
        CKmerAPI kmer;
        uint64_t c;
        uint32_t db_idx;
    };

    auto k = kmc_dbs[0]->KmerLength();
    for (int i = 1; i < kmc_dbs.size(); ++i)
        if (k != kmc_dbs[i]->KmerLength())
        {
            std::cerr << "Error: inconsistent k-mer length\n";
            return 1;
        }
    std::cerr << "k: " << k << "\n";

    auto start_bin_in_all_dbs = [&]() {
        bool is_next = kmc_dbs[0]->StartBin();
        for (size_t i = 1; i < kmc_dbs.size(); ++i)
            if (is_next != kmc_dbs[i]->StartBin()) {
                std::cerr << "Error: " << __FILE__ << ": " << __LINE__ << "\n";
                exit(1);
            }
        return is_next;
    };

    std::ofstream out(outpath);
    while (start_bin_in_all_dbs()) //while there is a next bin to read from
    {
        std::vector<KmerAndCount> tops;
        CKmerAPI kmer(k);
        uint64_t c;

        //collect first (lowest) k-mer (if any) from each input sample for current bin        
        for (int i = 0; i < kmc_dbs.size(); ++i)
        {
            if (kmc_dbs[i]->ReadNextKmerFromBin(kmer, c))
            {
                tops.emplace_back();
                tops.back().kmer = kmer;
                tops.back().c = c;
                tops.back().db_idx = i;
            }
        }

        while (tops.size())
        {
            size_t min_id = 0;
            //find id of lowest k-mer in current bin
            for (size_t i = 1; i < tops.size(); ++i)
            {
                if (tops[i].kmer < tops[min_id].kmer)
                    min_id = i;
            }

            kmer = tops[min_id].kmer;
            c = 0;

            //collect all k-mers equal to the lowest k-mer and their counters
            for (size_t i = 0; i < tops.size(); ++i)
            {
                if (tops[i].kmer == kmer)
                {
                    c += tops[i].c;
                    auto db_idx = tops[i].db_idx;
                    if (kmc_dbs[db_idx]->ReadNextKmerFromBin(tops[i].kmer, tops[i].c))
                        ;
                    else
                    {
                        tops[i--] = tops.back();
                        tops.pop_back();
                    }
                }
            }            
            out << kmer.to_string() << "\t" << c << "\n";
            
        }
    }

    for (auto& d : kmc_dbs)
        d->Close();
}