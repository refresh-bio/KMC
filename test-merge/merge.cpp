#include <iostream>
#include <fstream>
#include <vector>
#include <mutex>
#include <string>
#include <random>
#include <atomic>
#include <thread>
#include <memory>
#include <numeric>
#include "kmc_file.h"

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <outpath> <mapping_file> <n_threads> <input_1> [<input_2> ...]\n";
        return 0;
    }
    std::string outpath = argv[1];
    std::string mapping = argv[2];
    size_t n_threads = std::atoi(argv[3]);
    std::vector<std::string> inputs;
    for (int i = 4; i < argc; ++i)
        inputs.push_back(argv[i]);

    std::cerr << "mapping: " << mapping << "\n";
    std::cerr << "inputs:\n";
    for (auto& x : inputs)
        std::cerr << "\t" << x << "\n";

    
//
    struct KmerAndCount {
        CKmerAPI kmer;
        uint64_t c;
        uint32_t db_idx;
    };

    //auto k = kmc_dbs[0]->KmerLength();

    CKMCFile to_read_config;
    to_read_config.OpenForListingWithBinOrder(inputs[0], mapping);

    //this code is not checking if all inputs have the same base config (n_bins, k, signature length), but this is still required for correct results
    auto k = to_read_config.KmerLength();
    auto n_bins = to_read_config.GetNBins();
    


    //for (int i = 1; i < kmc_dbs.size(); ++i)
    //    if (k != kmc_dbs[i]->KmerLength())
    //    {
    //        std::cerr << "Error: inconsistent k-mer length\n";
    //        return 1;
    //    }
    //std::cerr << "k: " << k << "\n";

    std::vector<std::thread> threads;
    threads.reserve(n_threads);


    std::ofstream out(outpath);
    std::mutex out_mtx;

    std::atomic<uint32_t> bin_id{};
    for (size_t i = 0 ; i < n_threads ; ++i)
        threads.emplace_back([&]{
            
            std::vector<std::unique_ptr<CKMCFile>> kmc_dbs(inputs.size());

            for (int i = 0; i < inputs.size(); ++i)
            {
                kmc_dbs[i] = std::make_unique<CKMCFile>();
                if (!kmc_dbs[i]->OpenForListingWithBinOrder(inputs[i], mapping))
                {
                    std::cerr << "Error: cannot open kmc database " << inputs[i] << "\n";
                    exit(1);
                }
            }

            auto start_bin_in_all_dbs = [&](uint32_t bin_id) {
                for (size_t i = 0; i < kmc_dbs.size(); ++i)
                    kmc_dbs[i]->StartBin(bin_id);
            };
            
            //for simplicity I keep all thread out in memory
            //should be flushed sometimes in a real implementation
            std::ostringstream my_out;

            while (true)
            {
                uint32_t my_bin_id = bin_id++;

                if (my_bin_id >= n_bins)
                    break;

                out_mtx.lock();
                std::cerr << my_bin_id << "\n";
                out_mtx.unlock();

                start_bin_in_all_dbs(my_bin_id);

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

                    //I know this is very bad to lock mutex for each out operation
                    my_out << kmer.to_string() << "\t" << c << "\n";
                    
                }
            }
            
            for (auto& d : kmc_dbs)
                d->Close();

            auto my_out_str = my_out.str();
            std::lock_guard<std::mutex> lck(out_mtx);
            out << my_out_str;
            
        });

    for (auto& t: threads)
        t.join();
}