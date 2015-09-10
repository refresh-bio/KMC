/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 2.3.0
  Date   : 2015-08-21
*/

#include "stdafx.h"
#include <iostream>
#include <vector>

#include "config.h"
#include "parser.h"
#include "timer.h"
#include "kmc1_db_reader.h"
#include "kmc2_db_reader.h"
#include "kmc1_db_writer.h"
#include "parameters_parser.h"
#include "histogram_writer.h"
#include "dump_writer.h"
#include "fastq_reader.h"
#include "fastq_filter.h"
#include "fastq_writer.h"
#ifdef ENABLE_LOGGER
#include "develop.h"
#endif
using namespace std;

template<unsigned SIZE> class CTools
{
	CParametersParser& parameters_parser;
	CConfig& config;
	bool histo()
	{
		if (!config.headers.front().IsKMC2()) //KMC1
		{
			CKMC1DbReader<SIZE> kmcdb(config.headers.front(), config.input_desc.front(), CConfig::GetInstance().percent_progress, KMCDBOpenMode::counters_only);
			CHistogramWriter<CKMC1DbReader<SIZE>> writer(kmcdb);
			return writer.Process();
		}
		else //KMC2
		{
			CKMC2DbReader<SIZE> kmcdb(config.headers.front(), config.input_desc.front(), CConfig::GetInstance().percent_progress, KMCDBOpenMode::counters_only);
			CHistogramWriter<CKMC2DbReader<SIZE>> writer(kmcdb);
			return writer.Process();
		}
	}

	bool dump()
	{
		if (!config.headers.front().IsKMC2()) //KMC1 - input is sorted
		{			
			CKMCDBForDump<CKMC1DbReader<SIZE>, SIZE, true> kmcdb_wrapper;
			CDumpWriter<decltype(kmcdb_wrapper), SIZE> writer(kmcdb_wrapper);
			return writer.Process();
		}
		else //KMC2
		{
			if (config.dump_params.sorted_output)
			{			
				CKMCDBForDump<CKMC2DbReader<SIZE>, SIZE, true> kmcdb_wrapper;
				CDumpWriter<decltype(kmcdb_wrapper), SIZE> writer(kmcdb_wrapper);
				return writer.Process();
			}
			else
			{				
				CKMCDBForDump<CKMC2DbReader<SIZE>, SIZE, false> kmcdb_wrapper;
				CDumpWriter<decltype(kmcdb_wrapper), SIZE> writer(kmcdb_wrapper);
				return writer.Process();
			}		
		}
		return true;
	}
	
	bool filter()
	{
		CFilteringParams& filtering_params = config.filtering_params;
		CFilteringQueues filtering_queues;				

		//set parameters and quques
		int32 avaiable_threads = config.avaiable_threads;
		filtering_params.n_readers = max(1, avaiable_threads / 2);

		bool gz_bz2 = false;
		vector<uint64> file_sizes;

		for (auto& p : filtering_params.input_srcs)
		{
			string ext(p.end() - 3, p.end());
			if (ext == ".gz" || ext == ".bz2")
			{
				gz_bz2 = true;				
			}
			FILE* tmp = my_fopen(p.c_str(), "rb");
			if (!tmp)
			{
				cout << "Cannot open file: " << p.c_str();
				exit(1);
			}
			my_fseek(tmp, 0, SEEK_END);
			file_sizes.push_back(my_ftell(tmp));
			fclose(tmp);
		}
		if (gz_bz2)
		{
			sort(file_sizes.begin(), file_sizes.end(), greater<uint64>());
			uint64 file_size_threshold = (uint64)(file_sizes.front() * 0.05);
			int32 n_allowed_files = 0;
			for (auto& p : file_sizes)
			if (p > file_size_threshold)
				++n_allowed_files;
			filtering_params.n_readers = MIN(n_allowed_files, MAX(1, avaiable_threads / 2));
		}
		else
			filtering_params.n_readers = 1;



		avaiable_threads -= filtering_params.n_readers;
		filtering_params.n_filters = max(1, avaiable_threads);
		
		filtering_params.fastq_buffer_size = 1 << 25; 

		filtering_params.mem_part_pmm_fastq_reader = filtering_params.fastq_buffer_size + CFastqReader::OVERHEAD_SIZE;
		filtering_params.mem_tot_pmm_fastq_reader = filtering_params.mem_part_pmm_fastq_reader * (filtering_params.n_readers + 48);

		filtering_params.mem_part_pmm_fastq_filter = filtering_params.mem_part_pmm_fastq_reader;
		filtering_params.mem_tot_pmm_fastq_filter = filtering_params.mem_part_pmm_fastq_filter * (filtering_params.n_filters + 48);

		filtering_queues.input_files_queue = new CInputFilesQueue(filtering_params.input_srcs); 
		filtering_queues.input_part_queue = new CPartQueue(filtering_params.n_readers);
		filtering_queues.filtered_part_queue = new CPartQueue(filtering_params.n_filters);

		filtering_queues.pmm_fastq_reader = new CMemoryPool(filtering_params.mem_tot_pmm_fastq_reader, filtering_params.mem_part_pmm_fastq_reader); 
		filtering_queues.pmm_fastq_filter = new CMemoryPool(filtering_params.mem_tot_pmm_fastq_filter, filtering_params.mem_part_pmm_fastq_filter);


		filtering_params.kmer_len = config.headers.front().kmer_len;

		vector<thread> readers_ths;
		vector<thread> filters_ths;
		vector<unique_ptr<CWFastqFilter>> filters;
		vector<unique_ptr<CWFastqReader>> readers;

		CKMCFile kmc_api;
		if (!kmc_api.OpenForRA(config.input_desc.front().file_src))
		{
			cout << "Error: cannot open: " << config.input_desc.front().file_src << " by KMC API\n";
			exit(1);
		}
		kmc_api.SetMinCount(config.input_desc.front().cutoff_min);
		kmc_api.SetMaxCount(config.input_desc.front().cutoff_max);

		CWFastqWriter writer(filtering_params, filtering_queues);
		thread writer_th(writer);

		for (uint32 i = 0; i < filtering_params.n_filters; ++i)
		{
			filters.push_back(make_unique<CWFastqFilter>(filtering_params, filtering_queues, kmc_api));
			filters_ths.emplace_back(ref(*filters.back()));
		}

		for (uint32 i = 0; i < filtering_params.n_readers; ++i)
		{
			readers.push_back(make_unique<CWFastqReader>(filtering_params, filtering_queues));
			readers_ths.emplace_back(ref(*readers.back()));
		}

		writer_th.join();
		for (auto& thread : filters_ths)
			thread.join();

		filters.clear();


		for (auto& thread : readers_ths)
			thread.join();

		readers.clear();
				
		delete filtering_queues.input_part_queue;
		delete filtering_queues.pmm_fastq_reader;
		delete filtering_queues.pmm_fastq_filter;
		delete filtering_queues.input_files_queue;
		delete filtering_queues.filtered_part_queue;

		return true;
	}

public:
	CTools(CParametersParser& parameters_parser) :
		parameters_parser(parameters_parser),
		config(CConfig::GetInstance())
	{
	}
	bool Process()
	{
		if (config.mode == CConfig::Mode::FILTER)
		{
			return filter();
		}
		if (config.mode == CConfig::Mode::HISTOGRAM)
		{
			return histo();
		}
		else if (config.mode == CConfig::Mode::DUMP)
		{
			return dump();
		}
		else if (config.mode == CConfig::Mode::COMPARE)
		{
			CInput<SIZE> *db1, *db2;
			if (!config.headers[0].IsKMC2())
				db1 = new CKMC1DbReader<SIZE>(config.headers[0], config.input_desc[0], CConfig::GetInstance().percent_progress, KMCDBOpenMode::sorted);
			else
				db1 = new CKMC2DbReader<SIZE>(config.headers[0], config.input_desc[0], CConfig::GetInstance().percent_progress, KMCDBOpenMode::sorted);

			if (!config.headers[1].IsKMC2())
				db2 = new CKMC1DbReader<SIZE>(config.headers[1], config.input_desc[1], CConfig::GetInstance().percent_progress, KMCDBOpenMode::sorted);
			else
				db2 = new CKMC2DbReader<SIZE>(config.headers[1], config.input_desc[1], CConfig::GetInstance().percent_progress, KMCDBOpenMode::sorted);

			CBundle<SIZE> input1(db1), input2(db2);
			CComparer<SIZE> comparer(&input1, &input2);

			bool res = comparer.Equals();

			delete db1;
			delete db2;
			std::cout << "\n";
			if (res)
			{
				cout << "DB Equals\n";
				exit(0);
			}
			else
			{
				cout << "DB Differs\n";
				exit(1);
			}
		}
		else
		{
			CExpressionNode<SIZE>* expression_root = parameters_parser.GetExpressionRoot<SIZE>();
			auto t = expression_root->GetExecutionRoot();
			delete expression_root;
			CKMC1DbWriter<SIZE> writer(t);
			writer.Process();
			delete t;
			return true;
		}
		return false;
	}

	

};


template<unsigned SIZE> class CApplication
{
	CApplication<SIZE - 1>* app_1;
	CTools<SIZE>* tools;
	bool is_selected;
	CConfig& config;
	CParametersParser& parameter_parser;
public:
	CApplication(CParametersParser& parameter_parser) :
		config(CConfig::GetInstance()), parameter_parser(parameter_parser)
	{
		is_selected = config.kmer_len <= (int32)SIZE * 32 && config.kmer_len > ((int32)SIZE - 1) * 32;

		app_1 = new CApplication<SIZE - 1>(parameter_parser);
		if (is_selected)
		{
			tools = new CTools<SIZE>(parameter_parser);
		}
		else
		{
			tools = nullptr;
		}
	}

	~CApplication()
	{
		delete app_1;
		if (is_selected)
			delete tools;
	}

	bool Process()
	{
		if (is_selected)
			return tools->Process();
		else
			return app_1->Process();
	}
};

template<> class CApplication<1>
{
	CTools<1>* tools;
	CConfig& config;
	CParametersParser& parameter_parser;
	bool is_selected;
public:
	CApplication(CParametersParser& parameter_parser) :
		config(CConfig::GetInstance()), parameter_parser(parameter_parser)
	{
		is_selected = config.kmer_len <= 32;

		if (is_selected)
			tools = new CTools<1>(parameter_parser);
		else
			tools = nullptr;
	}
	~CApplication<1>()
	{
		if (tools)
			delete tools;
	}
	bool Process() {
		if (is_selected)
		{
			return tools->Process();
		}
		return false;
	}
};


int main(int argc, char**argv)
{
#ifdef ENABLE_LOGGER
	CTimer timer;
	timer.start();
#endif
	CParametersParser params_parser(argc, argv);
	params_parser.Parse();
	if (params_parser.validate_input_dbs())
	{
		params_parser.SetThreads();
		CApplication<KMER_WORDS> app(params_parser);
		app.Process();
	}

#ifdef ENABLE_LOGGER

	cout << "RUN TIME: " << timer.get_time() <<"ms\n\n";

	CLoger::GetLogger().print_stats();

#endif	
}


// ***** EOF