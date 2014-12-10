/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.0
  Date   : 2014-07-04
*/

#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <functional>
#include "timer.h"
#include "kmc.h"
#include "meta_oper.h"

using namespace std;

uint64 total_reads, total_fastq_size;

void usage();
bool parse_parameters(int argc, char *argv[]);

CKMCParams Params;

//----------------------------------------------------------------------------------
// Application class
// Template parameters:
//    * KMER_TPL - k-mer class
//    * SIZE     - maximal size of the k-mer (divided by 32)
template<template<unsigned X> class KMER_TPL, unsigned SIZE, bool QUAKE_MODE> class CApplication
{
	CApplication<KMER_TPL, SIZE - 1, QUAKE_MODE> *app_1;
	CKMC<KMER_TPL<SIZE>, SIZE, QUAKE_MODE> *kmc;
	int p_k;
	bool is_selected;

public:
	CApplication(CKMCParams &Params) {
		p_k = Params.p_k;
		is_selected = p_k <= (int32) SIZE * 32 && p_k > ((int32) SIZE-1)*32;

		app_1 = new CApplication<KMER_TPL, SIZE - 1, QUAKE_MODE>(Params);
		if(is_selected)
		{			
			kmc = new CKMC<KMER_TPL<SIZE>, SIZE, QUAKE_MODE>;
			kmc->SetParams(Params);
		}
		else
		{
			kmc = NULL;
		}
	};
	~CApplication() {
		delete app_1;
		if (kmc)
			delete kmc;
	}

	void GetStats(double &time1, double &time2, double &time3, uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uint64 &_tmp_size_strict_mem, uint64 &_max_disk_usage, uint64& _n_total_super_kmers) {
		if (is_selected)
		{
			kmc->GetStats(time1, time2, time3, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, _n_reads, _tmp_size, _tmp_size_strict_mem, _max_disk_usage, _n_total_super_kmers);
		}
		else
			app_1->GetStats(time1, time2, time3, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, _n_reads, _tmp_size, _tmp_size_strict_mem, _max_disk_usage, _n_total_super_kmers);
	}

	bool Process() {
		if (is_selected)
		{
			return kmc->Process();
		}
		else
			return app_1->Process();
	}
};

//----------------------------------------------------------------------------------
// Specialization of the application class for the SIZE=1
template<template<unsigned X> class KMER_TPL, bool QUAKE_MODE> class CApplication<KMER_TPL, 1, QUAKE_MODE>
{
	CKMC<KMER_TPL<1>, 1, QUAKE_MODE> *kmc;
	int p_k;
	bool is_selected;

public:
	CApplication(CKMCParams &Params) {
		is_selected = Params.p_k <= 32;
		if(is_selected)
		{
			kmc = new CKMC<KMER_TPL<1>, 1, QUAKE_MODE>;
			kmc->SetParams(Params);
		}
		else
		{
			kmc = NULL;
		}
	};
	~CApplication() {
		if(kmc)
			delete kmc;
	};

	void GetStats(double &time1, double &time2, double &time3, uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uint64 &_tmp_size_strict_mem, uint64 &_max_disk_usage, uint64& _n_total_super_kmers) {
		if (is_selected)
		{
			if(kmc)
				kmc->GetStats(time1, time2, time3, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, _n_reads, _tmp_size, _tmp_size_strict_mem, _max_disk_usage, _n_total_super_kmers);
		}
	}

	bool Process() {
		if (is_selected)
		{
			return kmc->Process();
		}
		return false;
	}
};


//----------------------------------------------------------------------------------
// Show execution options of the software
void usage()
{
	cout << "K-Mer Counter (KMC) ver. " << KMC_VER << " (" << KMC_DATE << ")\n";
	cout << "Usage:\n kmc [options] <input_file_name> <output_file_name> <working_directory>\n";
	cout << " kmc [options] <@input_file_names> <output_file_name> <working_directory>\n";
	cout << "Parameters:\n";
	cout << "  input_file_name - single file in FASTQ format (gziped or not)\n";
	cout << "  @input_file_names - file name with list of input files in FASTQ format (gziped or not)\n";
	cout << "Options:\n";
	cout << "  -v - verbose mode (shows all parameter settings); default: false\n";
	cout << "  -k<len> - k-mer length (k from " << MIN_K << " to " << MAX_K << "; default: 25\n";
	cout << "  -m<size> - max amount of RAM in GB (from 4 to 1024); default: 12\n";
	cout << "  -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)\n";
	cout << "  -p<par> - signature length (5, 6, 7, 8); default: 7\n";
	cout << "  -f<a/q/m> - input in FASTA format (-fa), FASTQ format (-fq) or mulit FASTA (-fm); default: FASTQ\n";
	cout << "  -q[value] - use Quake's compatible counting with [value] representing lowest quality (default: 33)\n";
	cout << "  -ci<value> - exclude k-mers occurring less than <value> times (default: 2)\n";
	cout << "  -cs<value> - maximal value of a counter (default: 255)\n";
	cout << "  -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)\n";
	cout << "  -b - turn off transformation of k-mers into canonical form\n";	
	cout << "  -r - turn on RAM-only mode \n";
	cout << "  -n<value> - number of bins \n";
	cout << "  -t<value> - total number of threads (default: no. of CPU cores)\n";
	cout << "  -sf<value> - number of FASTQ reading threads\n";
	cout << "  -sp<value> - number of splitting threads\n";
	cout << "  -sr<value> - number of sorter threads\n";
	cout << "  -so<value> - number of threads per single sorter\n";	
	cout << "Example:\n";
	cout << "kmc -k27 -m24 NA19238.fastq NA.res \\data\\kmc_tmp_dir\\\n";
	cout << "kmc -k27 -q -m24 @files.lst NA.res \\data\\kmc_tmp_dir\\\n";
}

//----------------------------------------------------------------------------------
// Parse the parameters
bool parse_parameters(int argc, char *argv[])
{
	int i;
	int tmp;

	if(argc < 4)
		return false;

	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;
		// Number of threads
		if(strncmp(argv[i], "-t", 2) == 0)
			Params.p_t = atoi(&argv[i][2]);
//		else 
		// k-mer length
		if(strncmp(argv[i], "-k", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if(tmp < MIN_K || tmp > MAX_K)
			{
				cout << "Wrong parameter: k must be from range <" << MIN_K << "," << MAX_K << ">\n";
				return false;
			}
			else
				Params.p_k = tmp;
		}
		// Memory limit
		else if(strncmp(argv[i], "-m", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if(tmp < MIN_MEM)
			{
				cout << "Wrong parameret: min memory must be at least " << MIN_MEM << "GB\n";
				return false;
			}
			else
				Params.p_m = tmp;
		}
		// Minimum counter threshold
		else if(strncmp(argv[i], "-ci", 3) == 0)
			Params.p_ci = atoi(&argv[i][3]);
		// Maximum counter threshold
		else if(strncmp(argv[i], "-cx", 3) == 0)
			Params.p_cx = atoi(&argv[i][3]);
		// Maximal counter value
		else if(strncmp(argv[i], "-cs", 3) == 0)
			Params.p_cs = atoi(&argv[i][3]);
		// Quake mode
		else if(strncmp(argv[i], "-q", 2) == 0)
		{
			Params.p_quake = true;
			if(strlen(argv[i]) > 2)
				Params.p_quality = atoi(argv[i]+2);
		}
		// Set p1
		else if (strncmp(argv[i], "-p", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if (tmp < MIN_SL || tmp > MAX_SL)
			{
				cout << "Wrong parameter: p must be from range <" << MIN_SL << "," << MAX_SL << ">\n";
				return false;
			}
			else
				Params.p_p1 = tmp;
		}
		// FASTA input files
		else if(strncmp(argv[i], "-fa", 3) == 0)
			Params.p_file_type = fasta;
		// FASTQ input files
		else if(strncmp(argv[i], "-fq", 3) == 0)
			Params.p_file_type = fastq;
		else if(strncmp(argv[i], "-fm", 3) == 0)
			Params.p_file_type = multiline_fasta;
		else if (strncmp(argv[i], "-v", 2) == 0)
			Params.p_verbose = true;
		else if (strncmp(argv[i], "-sm", 3) == 0 && strlen(argv[i]) == 3)
			Params.p_strict_mem = true;
		else if (strncmp(argv[i], "-r", 2) == 0)
			Params.p_mem_mode = true;
		else if(strncmp(argv[i], "-b", 2) == 0)
			Params.p_both_strands = false;
		// Number of reading threads
		else if(strncmp(argv[i], "-sf", 3) == 0)
		{
			tmp = atoi(&argv[i][3]);
			if(tmp < MIN_SF || tmp > MAX_SF)
			{
				cout << "Wrong parameter: number of reading thread must be from range <" << MIN_SF << "," << MAX_SF << ">\n";
				return false;
			}
			else
				Params.p_sf = tmp;
		}
		// Number of splitting threads
		else if(strncmp(argv[i], "-sp", 3) == 0)
		{
			tmp = atoi(&argv[i][3]);
			if(tmp < MIN_SP || tmp > MAX_SP)
			{
				cout << "Wrong parameter: number of splitting threads must be in range <" << MIN_SP << "," << MAX_SP << "<\n";
				return false;
			}
			else
				Params.p_sp = tmp;
		}
		// Number of sorting threads
		else if(strncmp(argv[i], "-so", 3) == 0)
		{
			tmp = atoi(&argv[i][3]);
			if(tmp < MIN_SO || tmp > MAX_SO)
			{	
				cout << "Wrong parameter: number of sorter threads must be in range <" << MIN_SO << "," << MAX_SO << "\n";
				return false;
			}
			else
				Params.p_so = tmp;
		}		
		// Number of internal sorting threads (per single sorter)
		else if(strncmp(argv[i], "-sr", 3) == 0)
		{
			tmp = atoi(&argv[i][3]);
			if(tmp < MIN_SR || tmp > MAX_SR)
			{
				cout << "Wrong parameter: number of sotring threads per single sorter must be in range <" << MIN_SR << "," << MAX_SR << "\n";
				return false;
			}
			else
				Params.p_sr = tmp;
		}
		else if (strncmp(argv[i], "-n", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if (tmp < MIN_N_BINS || tmp > MAX_N_BINS)
			{
				cout << "Wrong parameter: number of bins must be in range <" << MIN_SR << "," << MAX_SR << "\n";
				return false;
			}
			else
				Params.p_n_bins = tmp;
		}
	
		if (strncmp(argv[i], "-smso", 5) == 0)
		{
			tmp = atoi(&argv[i][5]);
			if (tmp < MIN_SMSO || tmp > MAX_SMSO)
			{
				cout << "Wrong parameter: number of sorting threads per sorter in strict memory mode must be in range <" << MIN_SMSO << "," << MAX_SMSO << "\n";
				return false;
			}
			else
				Params.p_smso = tmp;
		}
		if (strncmp(argv[i], "-smun", 5) == 0)
		{
			tmp = atoi(&argv[i][5]);
			if (tmp < MIN_SMUN || tmp > MAX_SMUN)
			{
				cout << "Wrong parameter: number of uncompactor threads in strict memory mode must be in range <" << MIN_SMUN << "," << MAX_SMUN << "\n";
				return false;
			}
			else
				Params.p_smun = tmp;
		}		
		if (strncmp(argv[i], "-smme", 5) == 0)
		{
			tmp = atoi(&argv[i][5]);
			if (tmp < MIN_SMME || tmp > MAX_SMME)
			{
				cout << "Wrong parameter: number of merger threads in strict memory mode must be in range <" << MIN_SMME << "," << MAX_SMME << "\n";
				return false;
			}
			else
				Params.p_smme = tmp;
		}
	}

	if(argc - i < 3)
		return false;

	string input_file_name = string(argv[i++]);
	Params.output_file_name = string(argv[i++]);
	Params.working_directory = string(argv[i++]);

	Params.input_file_names.clear();
	if(input_file_name[0] != '@')
		Params.input_file_names.push_back(input_file_name);
	else
	{
		ifstream in(input_file_name.c_str()+1);
		if(!in.good())
		{
			cout << "Error: No " << input_file_name.c_str()+1 << " file\n";
			return false;
		}

		string s;
		while(getline(in, s))
			if(s != "")
				Params.input_file_names.push_back(s);

		in.close();
		random_shuffle(Params.input_file_names.begin(), Params.input_file_names.end());
	}


	//Validate and resolve conflicts in parameters
	if (Params.p_strict_mem && Params.p_mem_mode)
	{
		cout << "Error: -sm can not be used with -r\n";
		return false;
	}
	if (Params.p_strict_mem && Params.p_quake)
	{
		cout << "Warning: -sm is not supported in quake mode. -sm has no effect\n";
		Params.p_strict_mem = false;
	}

	return true;
}

//----------------------------------------------------------------------------------
// Main function
int _tmain(int argc, _TCHAR* argv[])
{
	CStopWatch w0, w1;
	double time1, time2, time3;
	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total, n_reads, tmp_size, tmp_size_strict_mem, max_disk_usage, n_total_super_kmers;

	omp_set_num_threads(1);

#ifdef WIN32
	_setmaxstdio(2040);
#endif

	if(!parse_parameters(argc, argv))
	{
		usage();
		return 0;
	}

	if(Params.p_quake)
	{
		CApplication<CKmerQuake, KMER_WORDS, true> *app = new CApplication<CKmerQuake, KMER_WORDS, true>(Params);

		if(!app->Process())
		{
			cout << "Not enough memory or some other error\n";
			delete app;
			return 0;
		}
		app->GetStats(time1, time2, time3, n_unique, n_cutoff_min, n_cutoff_max, n_total, n_reads, tmp_size, tmp_size_strict_mem, max_disk_usage, n_total_super_kmers);
		delete app;
	}
	else
	{
		CApplication<CKmer, KMER_WORDS, false> *app = new CApplication<CKmer, KMER_WORDS, false>(Params);

		if(!app->Process())
		{
			cout << "Not enough memory or some other error\n";
			delete app;
			return 0;
		}
		app->GetStats(time1, time2, time3, n_unique, n_cutoff_min, n_cutoff_max, n_total, n_reads, tmp_size, tmp_size_strict_mem, max_disk_usage, n_total_super_kmers);
		delete app;
	}

	cout << "1st stage: " << time1 << "s\n";
	cout << "2nd stage: " << time2  << "s\n";
	if (Params.p_strict_mem)
		cout << "3rd stage: " << time3 << "s\n";
	if (Params.p_strict_mem)
		cout << "Total    : " << (time1 + time2 + time3) << "s\n";
	else
		cout << "Total    : " << (time1+time2) << "s\n";	
	if (Params.p_strict_mem)
	{
		cout << "Tmp size : " << tmp_size / 1000000 << "MB\n";
		cout << "Tmp size strict memory : " << tmp_size_strict_mem / 1000000 << "MB\n";
		//cout << "Tmp total: " << (tmp_size + tmp_size_strict_mem) / 1000000 << "MB\n";
		cout << "Tmp total: " << max_disk_usage / 1000000 << "MB\n";
	}
	else
		cout << "Tmp size : " << tmp_size / 1000000 << "MB\n";
	cout << "\nStats:\n";
	cout << "   No. of k-mers below min. threshold : " << setw(12) << n_cutoff_min << "\n";
	cout << "   No. of k-mers above max. threshold : " << setw(12) << n_cutoff_max << "\n";
	cout << "   No. of unique k-mers               : " << setw(12) << n_unique << "\n";
	cout << "   No. of unique counted k-mers       : " << setw(12) << n_unique-n_cutoff_min-n_cutoff_max << "\n";
	cout << "   Total no. of k-mers                : " << setw(12) << n_total << "\n";
if(Params.p_file_type != multiline_fasta)
	cout << "   Total no. of reads                 : " << setw(12) << n_reads << "\n";
else
	cout << "   Total no. of sequences             : " << setw(12) << n_reads << "\n";
	cout << "   Total no. of super-k-mers          : " << setw(12) << n_total_super_kmers << "\n";
	return 0;
}

// ***** EOF
