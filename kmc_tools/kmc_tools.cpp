/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.0
  Date   : 2018-05-10
*/

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <utility>

#include "config.h"
#include "check_kmer.h"
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

	template<typename KMCDB>
	void ProcessTransformOper(KMCDBOpenMode open_mode)
	{
		KMCDB* db = new KMCDB(config.headers.front(), config.input_desc.front(), config.percent_progress, open_mode);

		vector<CKMC1DbWriter<SIZE>*> kmc_db_writers;
		vector<CBundle<SIZE>*> bundles;
		vector<CDumpWriterForTransform<SIZE>> dump_writers;
		vector<CHistogramWriterForTransform> histogram_writers;

		for (auto& desc : config.transform_output_desc)
		{
			switch (desc.op_type)
			{
			case CTransformOutputDesc::OpType::COMPACT:
			case CTransformOutputDesc::OpType::REDUCE:
			case CTransformOutputDesc::OpType::SORT:
			case CTransformOutputDesc::OpType::SET_COUNTS:
				kmc_db_writers.push_back(new CKMC1DbWriter<SIZE>(nullptr, desc));
				kmc_db_writers.back()->MultiOptputInit();
				bundles.push_back(new CBundle<SIZE>(nullptr));
				break;
			case CTransformOutputDesc::OpType::DUMP:
				dump_writers.emplace_back(desc);
				dump_writers.back().Init();
				break;
			case CTransformOutputDesc::OpType::HISTOGRAM:
				histogram_writers.emplace_back(desc);
				histogram_writers.back().Init();
				break;
			}
		}
		CKmer<SIZE> kmer;
		uint32 counter;
		
		if (open_mode == KMCDBOpenMode::counters_only) // only histogram oper
		{
			while (db->NextCounter(counter))
			{
				for (auto& out : histogram_writers)
					out.PutCounter(counter);
			}
		}
		else if (open_mode == KMCDBOpenMode::sequential) //historam or dump only
		{
			while (db->NextKmerSequential(kmer, counter))
			{
				for (auto& out : histogram_writers)
					out.PutCounter(counter);
				for (auto& out : dump_writers)
					out.PutKmer(kmer, counter);
			}
		}
		else 
		{
			CBundle<SIZE> tmp_bundle(db);
			db = nullptr;
			while (!tmp_bundle.Finished())
			{
				for (auto& out : histogram_writers)
					out.PutCounter(tmp_bundle.TopCounter());
				for (auto& out : dump_writers)
					out.PutKmer(tmp_bundle.TopKmer(), tmp_bundle.TopCounter());

				for (uint32 i = 0; i < bundles.size(); ++i)
				{
					bundles[i]->Insert(tmp_bundle.TopKmer(), tmp_bundle.TopCounter());
					if (bundles[i]->Full())
					{
						kmc_db_writers[i]->MultiOptputAddResultPart(*bundles[i]);
					}
				}

				tmp_bundle.Pop();				
			}			
		}

		for (uint32 i = 0; i < bundles.size(); ++i)
		{
			if (!bundles[i]->Empty())
				kmc_db_writers[i]->MultiOptputAddResultPart(*bundles[i]);
		}

		for (auto& out : histogram_writers)
			out.Finish();
		for (auto& out : dump_writers)		
			out.Finish();			
		
		for (auto& out : kmc_db_writers)
		{
			out->MultiOptputFinish();
			delete out;
		}
		for (auto& b : bundles)
			delete b;
			
		delete db;
	}

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
	
	bool info()
	{
		auto header = CConfig::GetInstance().headers.front();

		std::cout << "k                 :  " << header.kmer_len <<"\n"
				  << "total k-mers      :  " << header.total_kmers<< "\n"
				  << "cutoff max        :  " << header.max_count << "\n"
				  << "cutoff min        :  " << header.min_count << "\n"
				  << "counter size      :  " << header.counter_size << " bytes\n"
				  << "mode              :  " << (header.mode ? "quality-aware counters" : "occurrence counters") << "\n"
				  << "both strands      :  " << (header.both_strands ? "yes" : "no") << "\n"
				  << "database format   :  " << (header.IsKMC2() ? "KMC2.x" : "KMC1.x") << "\n"
				  << "signature length  :  " << header.signature_len << "\n"
				  << "number of bins    :  " << header.no_of_bins << "\n"
				  << "lut_prefix_len    :  " << header.lut_prefix_len << "\n";
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
				cerr << "Error: cannot open file: " << p.c_str();
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
			cerr << "Error: cannot open: " << config.input_desc.front().file_src << " by KMC API\n";
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

	bool simple_set()
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


		vector<COutputBundle<SIZE>*> output_bundles;
		vector<CKMC1DbWriter<SIZE>*> writers;

		for (uint32 i = 0; i < config.simple_output_desc.size(); ++i)
		{
			writers.push_back(new CKMC1DbWriter<SIZE>(nullptr, config.simple_output_desc[i]));
			writers.back()->MultiOptputInit();
			output_bundles.push_back(new COutputBundle<SIZE>(config.simple_output_desc[i].op_type, config.simple_output_desc[i].counter_op, *writers.back()));
		}


		CSimpleOperation<SIZE> op(&input1, &input2, output_bundles);
		op.Process();

		for (auto& writer : writers)
		{
			writer->MultiOptputFinish();
			delete writer;
		}

		for (auto o : output_bundles)
			delete o;

		return true;
	}

	bool check()
	{
		CKmerCheck<SIZE> checker(config.headers.front(), config.input_desc.front());
		return checker.CheckKmer();

	}

	bool transform()
	{
		bool kmers_needed = false;
		bool sort_needed = false;

		//remove not valid
		if (!config.headers.front().IsKMC2())
		{
			auto it = config.transform_output_desc.begin();
			while (it != config.transform_output_desc.end())
			{
				if (it->op_type == CTransformOutputDesc::OpType::SORT)
				{
					it = config.transform_output_desc.erase(it);
					cerr << "Warning: input database is already sorted. Each sort operation will be omitted\n";
				}
				else ++it;
			}
		}
		if (!config.transform_output_desc.size())
		{
			return false;
		}


		for (auto& desc : config.transform_output_desc)
		{
			if (desc.op_type == CTransformOutputDesc::OpType::REDUCE || desc.op_type == CTransformOutputDesc::OpType::COMPACT || desc.op_type == CTransformOutputDesc::OpType::SORT || desc.op_type == CTransformOutputDesc::OpType::SET_COUNTS)
			{
				kmers_needed = true;
				sort_needed = true;
				break;
			}
			if (desc.op_type == CTransformOutputDesc::OpType::DUMP)
			{
				kmers_needed = true;
				if (desc.sorted_output)
				{
					sort_needed = true;
					break;
				}
			}
		}

		KMCDBOpenMode open_mode;
		if (!kmers_needed)
			open_mode = KMCDBOpenMode::counters_only;
		else if (sort_needed)
			open_mode = KMCDBOpenMode::sorted;
		else
			open_mode = KMCDBOpenMode::sequential;
		if (config.headers.front().IsKMC2())
			ProcessTransformOper<CKMC2DbReader<SIZE>>(open_mode);
		else
			ProcessTransformOper<CKMC1DbReader<SIZE>>(open_mode);

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
		else if (config.mode == CConfig::Mode::INFO)
		{
			return info();
		}
		else if (config.mode == CConfig::Mode::HISTOGRAM)
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
			//std::cout << "\n";
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
		else if (config.mode == CConfig::Mode::TRANSFORM)
		{
			return transform();
		}
		else if (config.mode == CConfig::Mode::SIMPLE_SET)
		{
			return simple_set();
		}
		else if (config.mode == CConfig::Mode::CHECK)
		{
			return check();
		}
		else
		{
			CExpressionNode<SIZE>* expression_root = parameters_parser.GetExpressionRoot<SIZE>();
			auto t = expression_root->GetExecutionRoot();
			delete expression_root;
			CKMC1DbWriter<SIZE> writer(t, config.output_desc);
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

int main(int argc, char**argv)
{
#ifdef ENABLE_LOGGER 
	CTimer timer;
	timer.start();
#endif 

	if (argc == 1 || help_or_version(argc, argv))
	{
		CGeneralUsageDisplayer{}.Display();
		return 0;
	}

	CParametersParser params_parser(argc, argv);
	params_parser.Parse();
	if (params_parser.validate_input_dbs())
	{
		params_parser.SetThreads();
		CApplication<KMER_WORDS> app(params_parser);
		app.Process();
		//cout << "\n";
	}

#ifdef ENABLE_LOGGER 
	cout << "RUN TIME: " << timer.get_time() <<"ms\n\n";
	CLoger::GetLogger().print_stats(); 
#endif
}


// ***** EOF