#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <memory>
#include <numeric>
#include "kmc_file.h"

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <outpath> <input_1> [<input_2> ...]\n";
        return 0;
    }
    std::string outpath = argv[1];
    std::vector<std::string> inputs;
    for (int i = 2; i < argc; ++i)
        inputs.push_back(argv[i]);

    std::cerr << "inputs:\n";
    for (auto& x : inputs)
        std::cerr << "\t" << x << "\n";

    std::vector<std::unique_ptr<CKMCFile>> kmc_dbs(inputs.size());

    for (int i = 0; i < inputs.size(); ++i)
    {
        kmc_dbs[i] = std::make_unique<CKMCFile>();
        if (!kmc_dbs[i]->OpenForListingWithBinOrder(inputs[i]))
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
    auto n_bins = kmc_dbs[0]->GetNBins();
    CKMCFileInfo info;
    kmc_dbs[0]->Info(info);
    auto sig_len = info.signature_len;
    auto sss = info.signature_selection_scheme;
    for (int i = 1; i < kmc_dbs.size(); ++i) {
        if (k != kmc_dbs[i]->KmerLength())
        {
            std::cerr << "Error: inconsistent k-mer length\n";
            return 1;
        }
        if (n_bins != kmc_dbs[i]->GetNBins())
        {
            std::cerr << "Error: inconsistent number of bins\n";
            return 1;
        }
        kmc_dbs[i]->Info(info);
        if (sig_len != info.signature_len)
        {
            std::cerr << "Error: inconsistent signature length \n";
            return 1;
        }
        if (sss != info.signature_selection_scheme)
        {
            std::cerr << "Error: inconsistent signature selection scheme\n";
            return 1;
        }
    }
    std::cerr << "k: " << k << "\n";

    //get bin sizes:
    std::vector<size_t> bin_sizes(n_bins);
    for (auto& kmc_db : kmc_dbs)
        for (uint32_t bin_id = 0; bin_id < n_bins; ++bin_id)
        {
            uint64_t n_kmers{};
            if(!kmc_db->GetNKmers(bin_id, n_kmers))
            {
                std::cerr << "Error: cannot get #kmers in bin " << bin_id << "\n";
                return 1;
            }
            bin_sizes[bin_id] += n_kmers;
        }

    std::vector<std::pair<size_t /*size*/, uint32_t/*bin_id*/>> from_largest;
    from_largest.reserve(n_bins);

    std::cerr << "bin sizes:\n";
    for (uint32_t bin_id = 0; bin_id < n_bins; ++bin_id)
    {
        from_largest.emplace_back(bin_sizes[bin_id], bin_id);
        std::cerr << bin_id << "\t" << bin_sizes[bin_id] << "\n";
    }

    size_t prefix_file_buff_size_bytes = 1ull << 20;
    size_t suffix_file_buff_size_bytes = 1ull << 20;

    
    auto start_bin_in_all_dbs = [&](uint32_t bin_id) {
        for (size_t i = 0; i < kmc_dbs.size(); ++i)
            kmc_dbs[i]->StartBin(bin_id, prefix_file_buff_size_bytes, suffix_file_buff_size_bytes);
    };

    std::ofstream out(outpath);
    
    //std::vector<uint32_t> random_bin_order(n_bins);
    //std::iota(random_bin_order.begin(), random_bin_order.end(), 0ul);
    //std::shuffle(random_bin_order.begin(), random_bin_order.end(), std::mt19937{});

    std::sort(from_largest.begin(), from_largest.end(), std::greater<>{});
    
    for(uint32_t bin_id = 0 ; bin_id < n_bins ; ++bin_id)
    {
        //start_bin_in_all_dbs(bin_id); //while there is a next bin to read from
        //start_bin_in_all_dbs(random_bin_order[bin_id]); //while there is a next bin to read from
        std::cerr << "Processing bin " << from_largest[bin_id].second << " of size: " << from_largest[bin_id].first << "\n";
        start_bin_in_all_dbs(from_largest[bin_id].second);
    
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