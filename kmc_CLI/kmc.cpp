#define _CRT_SECURE_NO_WARNINGS
#include "../kmc_core/kmc_runner.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
using namespace std;

struct CLIParams
{
	std::string jsonSummaryFileName;
	std::string estimatedHistogramFileName;
};

struct Params
{
	KMC::Stage1Params stage1Params;
	KMC::Stage2Params stage2Params;
	CLIParams cliParams;
};
//----------------------------------------------------------------------------------
// Show execution options of the software
void usage()
{
	cout << "K-Mer Counter (KMC) ver. " << KMC::CfgConsts::kmc_ver << " (" << KMC::CfgConsts::kmc_date << ")\n"
		<< "Usage:\n kmc [options] <input_file_name> <output_file_name> <working_directory>\n"
		<< " kmc [options] <@input_file_names> <output_file_name> <working_directory>\n"
		<< "Parameters:\n"
		<< "  input_file_name - single file in specified (-f switch) format (gziped or not)\n"
		<< "  @input_file_names - file name with list of input files in specified (-f switch) format (gziped or not)\n"
		<< "Options:\n"
		<< "  -v - verbose mode (shows all parameter settings); default: false\n"
		<< "  -k<len> - k-mer length (k from " << KMC::CfgConsts::min_k<< " to " << KMC::CfgConsts::max_k << "; default: 25)\n"
		<< "  -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12\n"
		<< "  -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)\n"
		<< "  -hc - count homopolymer compressed k-mers (approximate and experimental)\n"
		<< "  -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9\n"
		<< "  -f<a/q/m/bam/kmc> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC(-fkmc); default: FASTQ\n"
		<< "  -ci<value> - exclude k-mers occurring less than <value> times (default: 2)\n"
		<< "  -cs<value> - maximal value of a counter (default: 255)\n"
		<< "  -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)\n"
		<< "  -b - turn off transformation of k-mers into canonical form\n"
		<< "  -r - turn on RAM-only mode \n"
		<< "  -n<value> - number of bins \n"
		<< "  -t<value> - total number of threads (default: no. of CPU cores)\n"
		<< "  -sf<value> - number of FASTQ reading threads\n"
		<< "  -sp<value> - number of splitting threads\n"
		<< "  -sr<value> - number of threads for 2nd stage\n"
		<< "  -j<file_name> - file name with execution summary in JSON format\n"
		<< "  -w - without output\n"
		<< "  -o<kmc/kff> - output in KMC of KFF format; default: KMC\n"
		<< "  -hp - hide percentage progress (default: false)\n"
		<< "  -e<file_name> - estimage histogram of k-mers counts and work as usual\n"
		<< "  -E<file_name> - only estimage histogram of k-mers counts and stop after stage1\n"
		<< "Example:\n"
		<< "kmc -k27 -m24 NA19238.fastq NA.res /data/kmc_tmp_dir/\n"
		<< "kmc -k27 -m24 @files.lst NA.res /data/kmc_tmp_dir/\n";
}

//----------------------------------------------------------------------------------
// Check if --help or --version was used
bool help_or_version(int argc, char** argv)
{
	const string version = "--version";
	const string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

bool CanCreateFile(const string& path)
{
	FILE* f = fopen(path.c_str(), "wb");
	if (!f)
		return false;
	fclose(f);
	remove(path.c_str());
	return true;
}

bool CanCreateFileInPath(const string& path)
{
	static const string name = "kmc_test.bin"; //Some random name
	if (path.back() == '\\' || path.back() == '/')
		return CanCreateFile(path + name);
	else
		return CanCreateFile(path + '/' + name);
}
//----------------------------------------------------------------------------------
// Parse the parameters
bool parse_parameters(int argc, char* argv[], Params& params)
{
	KMC::Stage1Params& stage1Params = params.stage1Params;
	KMC::Stage2Params& stage2Params = params.stage2Params;
	CLIParams& cliParams = params.cliParams;
	int i;
	int tmp;

	bool was_sm = false;
	bool was_r = false;	
	if (argc < 4)
		return false;

	for (i = 1; i < argc; ++i)
	{
		if (argv[i][0] != '-')
			break;
		// Number of threads
		if (strncmp(argv[i], "-t", 2) == 0)
		{
			auto nThreads = atoi(&argv[i][2]);
			stage1Params.SetNThreads(nThreads); //TODO: what with stage2 in this case?
			stage2Params.SetNThreads(nThreads);
		}
		// k-mer length
		else if (strncmp(argv[i], "-k", 2) == 0)
			stage1Params.SetKmerLen(atoi(&argv[i][2]));		
		// Memory limit
		else if (strncmp(argv[i], "-m", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			stage1Params.SetMaxRamGB(tmp);
			stage2Params.SetMaxRamGB(tmp);
		}
		// Minimum counter threshold
		else if (strncmp(argv[i], "-ci", 3) == 0)
			stage2Params.SetCutoffMin(atoi(&argv[i][3]));
		// Maximum counter threshold
		else if (strncmp(argv[i], "-cx", 3) == 0)
			stage2Params.SetCutoffMax(atoll(&argv[i][3]));
		// Maximal counter value
		else if (strncmp(argv[i], "-cs", 3) == 0)
			stage2Params.SetCounterMax(atoll(&argv[i][3]));
		// Set p1
		else if (strncmp(argv[i], "-p", 2) == 0)
			stage1Params.SetSignatureLen(atoi(&argv[i][2]));
		//output type
		else if (strncmp(argv[i], "-o", 2) == 0)
		{
			if (strncmp(argv[i] + 2, "kff", 3) == 0)
				stage2Params.SetOutputFileType(KMC::OutputFileType::KFF);
			else if (strncmp(argv[i] + 2, "kmc", 3) == 0)
				stage2Params.SetOutputFileType(KMC::OutputFileType::KMC);
			else
			{
				std::cerr << "Error: unsupported output type: " << argv[i] << " (use -okff or -okmc)\n";
				exit(1);
			}
		}
		// FASTA input files
		else if (strncmp(argv[i], "-fa", 3) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::FASTA);
		// FASTQ input files
		else if (strncmp(argv[i], "-fq", 3) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::FASTQ);
		else if (strncmp(argv[i], "-fm", 3) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::MULTILINE_FASTA);
		else if (strncmp(argv[i], "-fbam", 5) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::BAM);
		else if (strncmp(argv[i], "-fkmc", 5) == 0)
			stage1Params.SetInputFileType(KMC::InputFileType::KMC);
#ifdef DEVELOP_MODE //TODO: reconsider !!!! 
		else if (strncmp(argv[i], "-vl", 3) == 0)
			Params.p_verbose_log = true;
#endif
		else if (strncmp(argv[i], "-v", 2) == 0)
		{
			static KMC::CerrVerboseLogger logger;
			stage1Params.SetVerboseLogger(&logger);			
		}
		else if (strncmp(argv[i], "-sm", 3) == 0 && strlen(argv[i]) == 3)
		{
			was_sm = true;
			stage2Params.SetStrictMemoryMode(true);
		}
		else if (strncmp(argv[i], "-hp", 3) == 0 && strlen(argv[i]) == 3)
		{
			static KMC::NullPercentProgressObserver nullPercentProgressObserver;
			stage1Params.SetPercentProgressObserver(&nullPercentProgressObserver);
		}
		else if (strncmp(argv[i], "-hc", 3) == 0 && strlen(argv[i]) == 3)
			stage1Params.SetHomopolymerCompressed(true);
		else if (strncmp(argv[i], "-r", 2) == 0)
		{
			stage1Params.SetRamOnlyMode(true);
			was_r = true;
		}
		else if (strncmp(argv[i], "-b", 2) == 0)
			stage1Params.SetCanonicalKmers(false);
		// Number of reading threads
		else if (strncmp(argv[i], "-sf", 3) == 0)
			stage1Params.SetNReaders(atoi(&argv[i][3]));
		// Number of splitting threads
		else if (strncmp(argv[i], "-sp", 3) == 0)
			stage1Params.SetNSplitters(atoi(&argv[i][3]));
		// Number of internal threads per 2nd stage
		else if (strncmp(argv[i], "-sr", 3) == 0)
			stage2Params.SetNThreads(atoi(&argv[i][3]));
		else if (strncmp(argv[i], "-n", 2) == 0)
			stage1Params.SetNBins(atoi(&argv[i][2]));
		else if (strncmp(argv[i], "-j", 2) == 0)
		{
			cliParams.jsonSummaryFileName = &argv[i][2];
			if (cliParams.jsonSummaryFileName == "")
				cerr << "Warning: file name for json summary file missed (-j switch)\n";
		}
		else if (strncmp(argv[i], "-e", 2) == 0 || strncmp(argv[i], "-E", 2) == 0)
		{
			cliParams.estimatedHistogramFileName = &argv[i][2];
			auto eE = argv[i][1]; //small or capital e
			if (cliParams.estimatedHistogramFileName == "")
				cerr << "Warning: file name for estimated histogram missed (-" << eE << " switch)\n";
			if (eE == 'e')
				stage1Params.SetEstimateHistogramCfg(KMC::EstimateHistogramCfg::ESTIMATE_AND_COUNT_KMERS);
			else // 'E'
				stage1Params.SetEstimateHistogramCfg(KMC::EstimateHistogramCfg::ONLY_ESTIMATE);
		}
		else if (strncmp(argv[i], "-w", 2) == 0)
			stage2Params.SetWithoutOutput(true);			

		if (strncmp(argv[i], "-smso", 5) == 0)
			stage2Params.SetStrictMemoryNSortingThreadsPerSorters(atoi(&argv[i][5]));		
		
		if (strncmp(argv[i], "-smun", 5) == 0)		
			stage2Params.SetStrictMemoryNUncompactors(atoi(&argv[i][5]));		
		if (strncmp(argv[i], "-smme", 5) == 0)		
			stage2Params.SetStrictMemoryNMergers(atoi(&argv[i][5]));		
	}

	if (argc - i < 3)
		return false;

	string input_file_name = string(argv[i++]);	
	stage2Params.SetOutputFileName(argv[i++]);
	stage1Params.SetTmpPath(argv[i++]);

	std::vector<std::string> input_file_names;	
	if (input_file_name[0] != '@')
		input_file_names.push_back(input_file_name);
	else
	{
		ifstream in(input_file_name.c_str() + 1);
		if (!in.good())
		{
			cerr << "Error: No " << input_file_name.c_str() + 1 << " file\n";
			return false;
		}

		string s;
		while (getline(in, s))
			if (s != "")
				input_file_names.push_back(s);

		in.close();
		std::random_shuffle(input_file_names.begin(), input_file_names.end());
	}
	stage1Params.SetInputFiles(input_file_names);

	//Validate and resolve conflicts in parameters
		
	if (was_sm && was_r)
	{
		cerr << "Error: -sm can not be used with -r\n";
		return false;
	}
	
	//Check if output files may be created and if it is possible to create file in specified tmp location
	if (!stage2Params.GetWithoutOutput())
	{		
		string pre_file_name = stage2Params.GetOutputFileName() + ".kmc_pre";
		string suff_file_name = stage2Params.GetOutputFileName() + ".kmc_suf";
		if (!CanCreateFile(pre_file_name))
		{
			cerr << "Error: Cannot create file: " << pre_file_name << "\n";
			return false;
		}
		if (!CanCreateFile(suff_file_name))
		{
			cerr << "Error: Cannot create file: " << suff_file_name << "\n";
			return false;
		}
	}
	if (!CanCreateFileInPath(stage1Params.GetTmpPath()))
	{
		cerr << "Error: Cannot create file in specified working directory: " << stage1Params.GetTmpPath() << "\n";
		return false;
	}
	return true;
}

void save_estimated_histogram(const std::string& fileName, const std::vector<uint64_t>& estimatedHistogram)
{
	if (fileName == "")
		return;

	std::ofstream out(fileName);
	if (!out)
	{
		std::cerr << "Warning: Cannot open file " << fileName << " to store estimated histogram";
		return;
	}
	for (uint32_t i = 1; i < estimatedHistogram.size(); ++i)
		out << i << "\t" << estimatedHistogram[i] << "\n";
}

void save_stats_in_json_file(const Params& params, const KMC::Stage1Results& stage1Results,	const KMC::Stage2Results& stage2Results)
{
	if (params.cliParams.jsonSummaryFileName == "")
		return;

	ofstream stats(params.cliParams.jsonSummaryFileName);

	if (!stats)
	{		
		std::cerr << "Warning: Cannot open file " << params.cliParams.jsonSummaryFileName << " to store summary of execution in JSON format";		
		return;
	}

	double time1 = stage1Results.time;
	double time2 = stage2Results.time;
	double time3 = stage2Results.timeStrictMem;

	bool display_strict_mem_stats = params.stage2Params.GetStrictMemoryMode() && stage1Results.wasSmallKOptUsed;

	stats << "{\n"
		<< "\t\"1st_stage\": \"" << time1 << "s\",\n"
		<< "\t\"2nd_stage\": \"" << time2 << "s\",\n";
	if (display_strict_mem_stats)
	{
		stats << "\t\"3rd_stage\": \"" << time3 << "s\",\n";
		stats << "\t\"Total\": \"" << (time1 + time2 + time3) << "s\",\n";
	}

	else
		stats << "\t\"Total\": \"" << (time1 + time2) << "s\",\n";


	if (display_strict_mem_stats)
	{
		stats << "\t\"Tmp_size\": \"" << stage1Results.tmpSize / 1000000 << "MB\",\n"
			<< "\t\"Tmp_size_strict_memory\": \"" << stage2Results.tmpSizeStrictMemory / 1000000 << "MB\",\n"
			<< "\t\"Tmp_total\": \"" << stage2Results.maxDiskUsage / 1000000 << "MB\",\n";
	}
	else
		stats << "\t\"Tmp_size\": \"" << stage1Results.tmpSize / 1000000 << "MB\",\n";

	stats << "\t\"Stats\": {\n";

	stats << "\t\t\"#k-mers_below_min_threshold\": " << stage2Results.nBelowCutoffMin << ",\n"
		<< "\t\t\"#k-mers_above_max_threshold\": " << stage2Results.nAboveCutoffMax << ",\n"
		<< "\t\t\"#Unique_k-mers\": " << stage2Results.nUniqueKmers << ",\n"
		<< "\t\t\"#Unique_counted_k-mers\": " << stage2Results.nUniqueKmers - stage2Results.nBelowCutoffMin - stage2Results.nAboveCutoffMax << ",\n"
		<< "\t\t\"#Total no. of k-mers\": " << stage2Results.nTotalKmers << ",\n";
	if(params.stage1Params.GetInputFileType() != KMC::InputFileType::MULTILINE_FASTA)	
		stats << "\t\t\"#Total_reads\": " << stage1Results.nSeqences << ",\n";
	else
		stats << "\t\t\"#Total_sequences\": " << stage1Results.nSeqences << ",\n";
	stats << "\t\t\"#Total_super-k-mers\": " << stage1Results.nTotalSuperKmers << "\n";

	stats << "\t}\n";
	stats << "}\n";
	stats.close();
}

void print_summary(
	const Params& params,	
	const KMC::Stage1Results& stage1Results,
	const KMC::Stage2Results& stage2Results)
{
	const KMC::Stage1Params& stage1Params = params.stage1Params;
	const KMC::Stage2Params& stage2Params = params.stage2Params;

	cout << "1st stage: " << stage1Results.time << "s\n"
		<< "2nd stage: " << stage2Results.time << "s\n";

	bool display_strict_mem_stats = stage2Params.GetStrictMemoryMode() && !stage1Results.wasSmallKOptUsed;
	if (display_strict_mem_stats)
	{
		cout << "3rd stage: " << stage2Results.timeStrictMem << "s\n";
		cout << "Total    : " << (stage1Results.time + stage2Results.time + stage2Results.timeStrictMem) << "s\n";
	}
	else
		cout << "Total    : " << (stage1Results.time + stage2Results.time) << "s\n";
	if (display_strict_mem_stats)
	{
		cout << "Tmp size : " << stage1Results.tmpSize / 1000000 << "MB\n"
			<< "Tmp size strict memory : " << stage2Results.tmpSizeStrictMemory / 1000000 << "MB\n"
			<< "Tmp total: " << stage2Results.maxDiskUsage / 1000000 << "MB\n";
	}
	else
		cout << "Tmp size : " << stage1Results.tmpSize / 1000000 << "MB\n";
	cout << "\nStats:\n"
		<< "   No. of k-mers below min. threshold : " << setw(12) << stage2Results.nBelowCutoffMin << "\n"
		<< "   No. of k-mers above max. threshold : " << setw(12) << stage2Results.nAboveCutoffMax << "\n"
		<< "   No. of unique k-mers               : " << setw(12) << stage2Results.nUniqueKmers << "\n"
		<< "   No. of unique counted k-mers       : " << setw(12) << stage2Results.nUniqueKmers - stage2Results.nBelowCutoffMin - stage2Results.nAboveCutoffMax << "\n"
		<< "   Total no. of k-mers                : " << setw(12) << stage2Results.nTotalKmers << "\n";
	if (stage1Params.GetInputFileType() != KMC::InputFileType::MULTILINE_FASTA)
		cout << "   Total no. of reads                 : " << setw(12) << stage1Results.nSeqences << "\n";
	else
		cout << "   Total no. of sequences             : " << setw(12) << stage1Results.nSeqences << "\n";
	cout << "   Total no. of super-k-mers          : " << setw(12) << stage1Results.nTotalSuperKmers << "\n";
}

//----------------------------------------------------------------------------------
// Main function
int main(int argc, char** argv)
{
	if (argc == 1 || help_or_version(argc, argv))
	{
		usage();
		return 0;
	}	
	try
	{
		Params params;		
		if (!parse_parameters(argc, argv, params))
		{
			usage();
			return 0;
		}

		KMC::Runner runner;
		auto stage1Results = runner.RunStage1(params.stage1Params);
		save_estimated_histogram(params.cliParams.estimatedHistogramFileName, stage1Results.estimatedHistogram);
		auto stage2Results = runner.RunStage2(params.stage2Params);
		print_summary(params, stage1Results, stage2Results);
		save_stats_in_json_file(params, stage1Results, stage2Results);
	}
	catch (const std::runtime_error& err)
	{
		std::cerr << err.what() << "\n";
		return 1;
	}
	
	return 0;

}
