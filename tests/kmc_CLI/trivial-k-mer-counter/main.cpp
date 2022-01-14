#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <unordered_map>

#include <zlib.h>
using namespace std;

struct Params
{
    vector<string> inputPaths;
    string outputPath;
    uint64_t kmerLen = 25;
    uint64_t cutoffMin = 2;
    uint64_t cutoffMax = 1000000000;
    uint64_t countMax = 255;
    bool canonical = true;
};


class Counter
{
    const Params& params;
    uint64_t nSeqs{};

    class LineReader
    {
        gzFile file;
        static constexpr uint64_t buff_capacity = 1ull << 24;
        uint64_t buff_pos = 0;
        uint64_t buff_size = 0;
        std::unique_ptr<uint8_t[]> buff;
        
        bool read_char(uint8_t& c)
        {
            if(buff_pos == buff_size)
            {
                buff_pos = 0;
                buff_size = gzfread(buff.get(), 1, buff_capacity, file);
            }
            if(buff_size == 0)
                return false;
            c = buff[buff_pos++];
            return true;
        }

        public:
        LineReader(gzFile file):
            file(file),
            buff(make_unique_for_overwrite<uint8_t[]>(buff_capacity))
            {
                
            }
        pair<bool, string> NextLine()
        {
            string line;
            uint8_t c;
            bool any = false;
            while (read_char(c))
            {
                if(c == '\r')
                    continue;
                any = true;
                if (c != '\n')
                    line.push_back((char)c);
                else 
                    break;
            }
            return make_pair(any, line);
        }
    };

    class SeqReader
    {
        LineReader& lineReader;
        bool is = false;
        string prevLine;
        
        string nextSeqFastq()
        {
            auto [is2, seq] = lineReader.NextLine();
            auto [is3, qual_header] = lineReader.NextLine();
            auto [is4, qual] = lineReader.NextLine();
            
            if (!is2 || !is3 || !is4)
            {
                cerr << "Error reading file " << __FILE__ << "\t" << __LINE__ << "\n";
                exit(1);
            }

            std::tie(is, prevLine) = lineReader.NextLine();

            return seq;
        }

        string nextSeqFasta()
        {
            string seq;
            while (true)
            {
                std::tie(is, prevLine) = lineReader.NextLine();
                if(!is)
                    return seq;
                
                if (prevLine.size() == 0)
                    return seq;

                if (prevLine[0] == '>')
                    return seq;

                seq += prevLine;                
            }   
        }

        public:
        SeqReader(LineReader& lineReader) :
            lineReader(lineReader)
        {
            std::tie(is, prevLine) = lineReader.NextLine();
        }
        pair<bool, string> NextSeq()
        {
            if (!is)
                return make_pair(false, "");
        
            if (prevLine.size() == 0)
            {
                cerr << "Error: something is wrong with this file " << __FILE__ << "\t" << __LINE__ << "\n";
                exit(1);
            }
            if (prevLine[0] == '@')
                return make_pair(true, nextSeqFastq());
            if (prevLine[0] == '>')
                return make_pair(true, nextSeqFasta());                
            cerr << "Error: unsupported file format. I don't know how to interpret line: " << prevLine << "\n";                
            cerr << __FILE__ << "\t" << __LINE__ << "\n";
            exit(1);        
        }
    };

    string get_rev_compl(const std::string& kmer)
    {
        string res;
        res.reserve(kmer.size());
        for (auto it = kmer.crbegin() ; it != kmer.crend() ; ++it)
        {            
            switch (*it)
            {
            case 'A': res.push_back('T'); break;
            case 'C': res.push_back('G'); break;
            case 'G': res.push_back('C'); break;
            case 'T': res.push_back('A'); break;
            default:
                cerr << "Error: wrong symbol in k-mer: " << kmer << "\n";
                exit(1);
            }
        }
        return res;
    }
    void canonicalize(string& kmer)
    {
        auto rev = get_rev_compl(kmer);
        if(rev < kmer)
            kmer = rev;
    }

    std::unordered_map<string, uint64_t> m;

    void processFile(const std::string& fname)
    {
        cerr << "Processing file " << fname << "\n";
        auto gzfile = gzopen(fname.c_str(), "r");
        if(!gzfile)
        {
            cerr << "Error: cannot open file " << fname << "\n";
            exit(1);
        }
        LineReader lineReader(gzfile);
        SeqReader seqReader(lineReader);

        for (auto [is, seq] = seqReader.NextSeq(); is; tie(is, seq) = seqReader.NextSeq())
        {
            //if (nSeqs % 100000 == 0)
            //    cout << nSeqs << "\n";
            ++nSeqs;

            if (seq.size() < params.kmerLen)
                continue;
            //cout << "seq:\n" << seq << "\n";            
            
            
            transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
            uint64_t endPos = seq.size() - params.kmerLen + 1;            
            for (uint64_t pos = 0; pos < endPos ; ++pos)
            {
                auto kmer = seq.substr(pos, params.kmerLen);
                auto er_symb = kmer.find_first_not_of("ACGT");
                if (er_symb != std::string::npos)
                {                
                    pos += er_symb;
                    continue;
                }
                if (params.canonical)
                    canonicalize(kmer);                
                ++m[kmer];
            }
        }
        gzclose(gzfile);
    }

    public:
    Counter(const Params& params):
        params(params)
        {
            for (const auto& fname : params.inputPaths)            
               processFile(fname);
            
            vector<pair<string, uint64_t>> kmers;
            cerr << "Take k-mers from unordered map to vector...";
            uint64_t nCutoffMin{}, nCutoffMax{};
            uint64_t tot_kmers{};
            for (const auto& [kmer, count] : m)
            {
                tot_kmers += count;
                if (count < params.cutoffMin)
                    ++nCutoffMin;
                else if(count > params.cutoffMax)
                    ++nCutoffMax;
                else 
                    kmers.emplace_back(kmer, count > params.countMax ? params.countMax : count);
            }
            cerr << "\n";

            ofstream stats_file(params.outputPath + ".stats");
            if(!stats_file)
            {
                cerr << "Error: cannot open file " << params.outputPath + ".stats" << "\n";
                exit(1);
            }
            stats_file << "Total unique k-mers below min cutoff:\t" << nCutoffMin << "\n";
            stats_file << "Total unique k-mers above max cutoff:\t" << nCutoffMax << "\n";            
            stats_file << "Total unique k-mers:\t" << m.size() << "\n";            
            stats_file << "Total unique counted k-mers:\t" << kmers.size() << "\n";
            stats_file << "Total k-mers:\t" << tot_kmers << "\n";            
            stats_file << "Total sequences:\t" << nSeqs << "\n";

            cerr << "Sorting k-mers...";
            sort(kmers.begin(), kmers.end());
            cerr << "\n";

            ofstream out(params.outputPath);
            if(!out)
            {
                cerr << "Error: cannot open file " << params.outputPath << "\n";
                exit(1);
            }

            for (const auto& [kmer, count] : kmers)
                out << kmer << "\t" << count << "\n";

            
        }
    
};

template<typename T>
void ParseToType(const std::string& parseFrom, T& parseTo)
{
    istringstream s(parseFrom);
    if(!(s >> parseTo))
    {
        std::cerr << "Error: cannot parse " << parseFrom << "\n";
        exit(1);
    }
}

Params parse_params(int argc, char**argv)
{
    if(argc < 3)
    {
        cerr << "Usage: " << argv[0] << " [-k<k>] [-ci<ci>] [-cx<cx>] [-cs<cs>] [-b] [@]<inputFile(s)> <outputFile>\n";
        exit(1);
    }
    Params p;
    for (int i = 1 ; i < argc - 2; ++i)
    {
        string param = argv[i];
        if (param == "-k")
            ParseToType(argv[++i], p.kmerLen);
        else if (param == "-ci")
            ParseToType(argv[++i], p.cutoffMin);
        else if (param == "-cx")
            ParseToType(argv[++i], p.cutoffMax);
        else if (param == "-cs")
            ParseToType(argv[++i], p.countMax);            
        else if (param == "-b")
            p.canonical = false;
        else
        {
            cerr << "Error: unknown parameter: " << param << "\n";
            exit(1);
        }
    }

    p.outputPath = argv[argc-1];
    string inputPath = argv[argc-2];
    if (inputPath[0] == '@')
    {
        ifstream in(inputPath);
        if(!in)
        {
            cerr << "Error: cannot open file " << inputPath <<"\n";
            exit(1);
        }
        string path;
        while (getline(in, path))
            p.inputPaths.push_back(path);
    }
    else
        p.inputPaths.push_back(inputPath);
    

    return p;
}

int main(int argc, char**argv)
{
    auto params = parse_params(argc, argv);
    Counter counter(params);
}