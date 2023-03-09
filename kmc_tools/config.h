/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.2.2
  Date   : 2023-03-09
*/

#ifndef _CONFIG_H
#define _CONFIG_H

#include "defs.h"
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include "kmer_file_header.h"
#include "percent_progress.h"

enum class CounterOpType { MIN, MAX, SUM, DIFF, FROM_DB1, FROM_DB2, NONE };
enum class OutputType { KMC1, KFF1 };
enum class KmerDBOpenMode { sequential, sorted, counters_only };

class CCounterBuilder
{
public:
	static void build_counter(uint32& counter, uchar *&buffer, uint32 counter_bytes, uint32 counter_mask, bool little_endian)
	{
		if (counter_bytes == 0)
		{
			counter = 1;
			return;
		}
		memcpy(&counter, buffer, sizeof(counter));
		if (!little_endian)
			counter = _bswap_uint32(counter);
		counter &= counter_mask;
		buffer += counter_bytes;
	}
};

struct CDescBase
{
	std::string file_src;
	uint32 cutoff_min = 0; //0 means it is not set yet
	uint64 cutoff_max = 0; //0 means it is not set yet
	CDescBase(const std::string& file_src) :
		file_src(file_src)
	{
	}
	CDescBase() = default;
};

//************************************************************************************************************
// CInputDesc - description of a single input KMC database.
//************************************************************************************************************
struct CInputDesc : public CDescBase
{	
	uint32 threads = 0; //for kmc2 input
	CInputDesc(const std::string& file_src) :
		CDescBase(file_src)
	{

	}	
	CInputDesc() = default;
};

//************************************************************************************************************
// COutputDesc - description of a output KMC database.
//************************************************************************************************************
struct COutputDesc : public CDescBase
{
	uint32 counter_max = 0; //0 means it is not set yet
	uint64 counter_value = 0; //only for SET_COUNTER operation, 0 means not set yet
	uint8_t encoding = 0b00011011; //for KFF writers
	OutputType output_type = OutputType::KMC1;
	COutputDesc(const std::string& file_src) :
		CDescBase(file_src)		
	{

	}

	COutputDesc() = default;
};

struct CSimpleOutputDesc : public COutputDesc
{
	enum class OpType { INTERSECT, UNION, KMERS_SUBTRACTION, COUNTERS_SUBTRACTION, REVERSE_KMERS_SUBTRACTION, REVERSE_COUNTERS_SUBTRACTION }; 	
	OpType op_type;
	CounterOpType counter_op;
	CSimpleOutputDesc(OpType op_type) :
		op_type(op_type)		
	{		
		switch (op_type)
		{
		case OpType::INTERSECT:
			counter_op = CounterOpType::MIN; 
			break;
		case OpType::UNION:
			counter_op = CounterOpType::SUM;
			break;
		case OpType::KMERS_SUBTRACTION: //irrelevant for KMERS_SUBTRACTION and REVERSE_KMERS_SUBTRACTION
		case OpType::COUNTERS_SUBTRACTION:
		case OpType::REVERSE_COUNTERS_SUBTRACTION:			
		case OpType::REVERSE_KMERS_SUBTRACTION:
			counter_op = CounterOpType::DIFF;
			break;
		}
	}
};

struct CTransformOutputDesc : public COutputDesc
{
	enum class OpType { HISTOGRAM, DUMP, SORT, REDUCE, COMPACT, SET_COUNTS };
	OpType op_type;
	bool sorted_output = false; //only for dump operation, rest is sorted anyway (except histo which does not print k-mers at all)
	CTransformOutputDesc(OpType op_type) :op_type(op_type)
	{ }
};

struct CFilteringParams
{	
	enum class file_type { fasta, fastq };
	enum class FilterMode {normal, trim, hard_mask};
	uint32 n_readers;
	uint32 n_filters;
	int fastq_buffer_size;
	int64 mem_part_pmm_fastq_reader;
	int64 mem_tot_pmm_fastq_reader;

	int64 mem_part_pmm_fastq_filter;
	int64 mem_tot_pmm_fastq_filter;

	uint32 kmer_len;
	uint32 gzip_buffer_size = 64 << 20;

	std::vector<std::string> input_srcs;
	bool use_float_value = false;
	uint32 n_min_kmers = 2;
	uint32 n_max_kmers = 1000000000;
	float f_min_kmers = 0.0f;
	float f_max_kmers = 1.0f;
	file_type input_file_type = file_type::fastq;
	file_type output_file_type = file_type::fastq;
	std::string output_src;
	
	FilterMode filter_mode = FilterMode::normal;
};


struct CCheckParams
{
	std::string kmer;
};


//************************************************************************************************************
// CConfig - configuration of current application run. Singleton class.
//************************************************************************************************************
class CConfig
{
public:	
	enum class Mode { UNDEFINED, COMPLEX, COMPARE, FILTER, SIMPLE_SET, TRANSFORM, INFO, CHECK };
	uint32 avaiable_threads;
	uint32 kmer_len = 0;
	Mode mode = Mode::UNDEFINED;
	bool verbose = false;	
	std::vector<CInputDesc> input_desc;
	std::vector<CKmerFileHeader> headers;
		
	COutputDesc output_desc; //complex only?

	std::vector<CSimpleOutputDesc> simple_output_desc;

	CFilteringParams filtering_params; //for filter operation only	
	CCheckParams check_params; // for check operation only

	std::vector<CTransformOutputDesc> transform_output_desc;

	CPercentProgress percent_progress;

	static CConfig& GetInstance()
	{
		static CConfig config;
		return config;
	}		
	CConfig(const CConfig&) = delete;
	CConfig& operator=(const CConfig&) = delete;


	bool IsLittleEndian()
	{
		uint64 a = 0x0807060504030201;
		return ((uchar*)&a)[0] == 1 && ((uchar*)&a)[1] == 2 && ((uchar*)&a)[2] == 3 && ((uchar*)&a)[3] == 4 && ((uchar*)&a)[5] == 6 && ((uchar*)&a)[7] == 8;
	}

	bool IsSeparateThreadForMainProcessingNeeded()
	{
		return mode == Mode::SIMPLE_SET || mode == Mode::COMPLEX || mode == Mode::COMPARE;
	}

	std::string GetOperationName()
	{
		switch (mode)
		{
		case CConfig::Mode::UNDEFINED:
			return "";
		case CConfig::Mode::COMPLEX:
			return "complex";
		case CConfig::Mode::COMPARE:
			return "compare";
		case CConfig::Mode::FILTER:
			return "filter";
		case CConfig::Mode::SIMPLE_SET:
			return "simple";
		case CConfig::Mode::TRANSFORM:
			return "transform";
		case CConfig::Mode::CHECK:
			return "check";
		default:
			return "";
		}
	}

private:
	CConfig() = default;
};


class CUsageDisplayer
{
protected:
	std::string name;		
	bool is2ArgOper = false;
	bool is1ArgOper = false;
	CUsageDisplayer(const std::string& name) :name(name){}
	void Display2ArgGeneral() const
	{
		std::cout << "The '" << name << "' is two arguments' operation. General syntax:\n"
				  << " kmc_tools " << name << " <input1 [input1_params]> <input2 [input2_params]> <output [output_params]>\n"
				  << " input1, input2 - paths to databases generated by KMC \n"
				  << " output		  - path  to output database\n"
				  << " For each input there are additional parameters:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n"
				  << " For output there are additional parameters:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n"
				  << "  -cs<value> - maximal value of a counter\n";
	}

	void Display1ArgGeneral(bool output_params) const
	{
		std::cout << " The '" << name << "' is one argument operation. General syntax:\n"
				  << " kmc_tools " << name << " <input> [input_params] <output> "<< (output_params ? "[output_params]" : "") << "\n"
				  << " input - path to database generated by KMC \n"
				  << " For input there are additional parameters:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
	}
public:
	virtual void Display() const = 0;
	virtual ~CUsageDisplayer() {}
};

class CGeneralUsageDisplayer : public CUsageDisplayer
{
public:
	CGeneralUsageDisplayer() :CUsageDisplayer("")
	{}
	
	void Display() const override
	{
		std::cout << "kmc_tools ver. " << KMC_VER << " (" << KMC_DATE << ")\n"
				  << "Usage:\n kmc_tools [global parameters] <operation> [operation parameters]\n"
				  << "Available operations:\n"
				  << "  transform            - transforms single KMC's database\n"
				  << "  simple               - performs set operation on two KMC's databases\n"
				  << "  complex              - performs set operation on multiple KMC's databases\n"
				  << "  filter               - filter out reads with too small number of k-mers\n"
				  << " global parameters:\n"
				  << "  -t<value>            - total number of threads (default: no. of CPU cores)\n"
				  << "  -v                   - enable verbose mode (shows some information) (default: false)\n"
				  << "  -hp                  - hide percentage progress (default: false)\n"
				  << "Example:\n"
				  << "kmc_tools simple db1 -ci3 db2 -ci5 -cx300 union db1_union_db2 -ci10\n"
				  << "For detailed help of concrete operation type operation name without parameters:\n"
				  << "kmc_tools simple\n";
	}
};

class CSimpleOperationUsageDisplayer : public CUsageDisplayer
{
public:
	CSimpleOperationUsageDisplayer() : CUsageDisplayer("simple"){}
	void Display() const override
	{
		std::cout << "simple operation performs set operation on two input KMC's databases\n"
				  << "General syntax:\n"
				  << "kmc_tools simple <input1 [input1_params]> <input2 [input2_params]> <oper1 output1 [output_params1]> [<oper2 output2 [output_params2]> ...]\n"
				  << "input1, input2                    - paths to databases generated by KMC\n"
				  << "oper1, oper2, ..., operN          - set operations to be performed on input1 and input2\n"
				  << "output1, output2, ..., outputN    - paths to output k-mer databases\n"
				  << " Available operations:\n"
				  << "  intersect                  - output database will contains only k-mers that are present in both input sets\n"
				  << "  union                      - output database will contains each k-mer present in any of input sets\n"
				  << "  kmers_subtract             - difference of input sets based on k-mers. \n" 
				  << "                               Output database will contains only k-mers that are present in first input set but absent in the second one\n"
				  << "  counters_subtract          - difference of input sets based on k-mers and their counters (weaker version of kmers_subtract).\n"
				  << "                               Output database will contains all k-mers that are present in first input, \n"
				  << "                               beyond those for which counter operation will lead to remove (i.e. counter equal to 0 or negative number)\n"
				  << "  reverse_kmers_subtract     - same as kmers_subtract but treat input2 as first and input1 as second\n"
				  << "  reverse_counters_subtract  - same as counters_subtract but treat input2 as first and input1 as second\n"
				  << " For each input there are additional parameters:\n"
				  << "  -ci<value>  - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value>  - exclude k-mers occurring more of than <value> times\n"
				  << " For each output there are additional parameters:\n"
				  << "  -ci<value>  - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value>  - exclude k-mers occurring more of than <value> times\n"
				  << "  -cs<value>  - maximal value of a counter\n"
				  << "  -o<kmc|kff> - output in KMC or KFF format (default: kmc) \n"
				  << "  -oc<value>  - redefine counter calculation mode for equal k-mers\n"
				  << "    Available values : \n"
				  << "      min   - get lower value of a k-mer counter (default value for intersect operation)\n"
				  << "      max   - get upper value of a k-mer counter\n"
				  << "      sum   - get sum of counters from both databases\n"
				  << "      diff  - get difference between counters (default for counters_subtract operation)\n"
				  << "      left  - get counter from first database (input1)\n"
				  << "      right - get counter from second database (input2)\n"
				  << "Example:\n"
				  << "kmc_tools simple kmers1 -ci3 -cx70000 kmers2 union kmers1_kmers2_union -cs65536 -ocfirst intersect intersect_kmers1_kmers2 intersect intersect_max_kmers1_kmers2 -ocmax\n";
	}
};


class CTransformOperationUsageDisplayer : public CUsageDisplayer
{
public:
	CTransformOperationUsageDisplayer() : CUsageDisplayer("transform"){}
	void Display() const override
	{
		std::cout << "transform operation transforms single input database to output (text file or KMC database)\n"
				  << "General syntax:\n"
				  << "kmc_tools transform <input> [input_params] <oper1 [oper_params1] output1 [output_params1]> [<oper2 [oper_params2] output2 [output_params2]>...]\n"
				  << "input - path to database generated by KMC \n"
				  << "oper1, oper2, ..., operN          - transform operation name\n"
				  << "output1, output2, ..., outputN    - paths to output\n"
				  
				  << " Available operations:\n"
				  << "  sort                       - converts database produced by KMC2.x to KMC1.x database format (which contains k-mers in sorted order)\n"
				  << "  reduce                     - exclude too rare and too frequent k-mers\n"
				  << "  compact                    - remove counters of k-mers\n"
				  << "  histogram                  - produce histogram of k-mers occurrences\n"
				  << "  dump                       - produce text dump of kmc database\n"
				  << "  set_counts <value>         - set all k-mer counts to specific value\n"
				  
				  << " For input there are additional parameters:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n"
				  
				  << " For sort and reduce operations there are additional output_params:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n"
				  << "  -cs<value> - maximal value of a counter\n"

				  << " For compact, reduce, set_counts and sort operations is an additional output_params:\n"
				  << "  -o<kmc|kff> - output in KMC or KFF format (default: kmc) \n"				  
			
				  << " For histogram operation there are additional output_params:\n"
				  << "  -ci<value> - minimum value of counter to be stored in the otput file\n"
				  << "  -cx<value> - maximum value of counter to be stored in the otput file\n"

				  << " For dump operation there are additional oper_params:\n"
				  << "  -s - sorted output\n"
				  << "Example:\n"
				  << "kmc_tools transform db reduce err_kmers -cx10 reduce valid_kmers -ci11 histogram histo.txt dump dump.txt\n";
	}
};

class CComplexUsageDisplayer : public CUsageDisplayer
{
public:
	CComplexUsageDisplayer() :CUsageDisplayer("complex")
	{}
	void Display() const override
	{
		std::cout << "Complex operation allows one to define operations for more than 2 input k-mers sets. Command-line syntax:\n"
				  << "kmc_tools complex <operations_definition_file>\n"
				  << " operations_definition_file - path to file which define input sets and operations. It is text file with following syntax:\n"
				  << " ______________________________________________________________________________ \n"
				  << "|INPUT:                                                                        |\n"
				  << "|<input1>=<input1_db_path> [params]                                            |\n"
				  << "|<input2>=<input2_db_path> [params]                                            |\n"
				  << "|...                                                                           |\n"
				  << "|<inputN>=<inputN_db_path> [params]                                            |\n"
				  << "|OUTPUT:                                                                       |\n"
				  << "|<out_db>=<ref_input><oper[c_mode]><ref_input>[<oper[c_mode]><ref_input>[...]  |\n"
				  << "|[OUTPUT_PARAMS:                                                             __|\n"
				  << "|<output_params>]                                                           |  /\n"
				  << "|                                                                           | / \n"
				  << "|___________________________________________________________________________|/  \n"
				  << "input1, input2, ..., inputN - names of inputs used to define equation\n"
				  << "input1_db, input2_db_path, ..., inputN_db_path - paths to k-mers sets\n"
				  << "For each input there are additional parameters which can be set:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n"
				  << "out_db_path - path to output database\n"
				  << "ref_input - one of input1, input2, ..., inputN\n"
				  << "oper - one of {*,-,~,+}, which refers to {intersect, kmers_subtract, counters_subtract, union}\n"
				  << "operator * has the highest priority. Other operators has equal priorities. Order of operations can be changed with parentheses\n"
				  << "for {*,~,+} it is possible to redefine counter calculation mode ([c_mode]). Available values: min, max, diff, sum, left, right (detailet description available in simple help message)\n"
				  << "output_params are:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n"
				  << "  -cs<value> - maximal value of a counter\n"
				  
				  << "Example:\n"
				  << " __________________________________________________________________ \n"
				  << "|INPUT:                                                            |\n"
				  << "|set1 = kmc_o1 -ci5                                                |\n"
				  << "|set2 = kmc_o2                                                     |\n"
				  << "|set3 = kmc_o3 -ci10 -cx100                                      __|\n"
				  << "|OUTPUT:                                                        |  /\n"
				  << "|result = (set3 + min set1) * right set2                        | / \n"
				  << "|_______________________________________________________________|/  \n";
		
	}
};

class CFilterUsageDisplayer : public CUsageDisplayer
{
public:
	CFilterUsageDisplayer() : CUsageDisplayer("filter")
	{}
	void Display() const override
	{
		std::cout << " The '" << name << "' is two arguments' operation. General syntax:\n"
				  << " kmc_tools " << name << " [filter_params] <kmc_input_db> [kmc_input_db_params] <input_read_set> [input_read_set_params] <output_read_set> [output_read_set_params]\n"
				  << " filter_params:\n"
				  << " -t               - trim reads on first invalid k-mer instead of remove entirely\n"
				  << " -hm              - hard mask invalid k-mers in a read\n"
				  << " kmc_input_db     - path to database generated by KMC \n"
				  << " input_read_set   - path to input set of reads \n"
				  << " output_read_set  - path to output set of reads \n"
				  << " For k-mers' database there are additional parameters:\n"
				  << "  -ci<value> - exclude k-mers occurring less than <value> times \n"
				  << "  -cx<value> - exclude k-mers occurring more of than <value> times\n"
				  << " For input set of reads there are additional parameters:\n"
				  << " -ci<value> - remove reads containing less k-mers than value. It can be integer or floating number in range [0.0;1.0] (default: 2)\n"
				  << " -cx<value> - remove reads containing more k-mers than value. It can be integer or floating number in range [0.0;1.0] (default: 1e9)\n"
				  << " -f<a/q>    - input in FASTA format (-fa), FASTQ format (-fq); default: FASTQ\n"
				  << " For output set of reads there are additional parameters:\n"
				  << " -f<a/q>    - output in FASTA format (-fa), FASTQ format (-fq); default: same as input\n"
				  << "Example:\n"
				  << "kmc_tools filter kmc_db -ci3 input.fastq -ci0.5 -cx1.0 filtered.fastq\n"
				  << "kmc_tools filter kmc_db input.fastq -ci10 -cx100 filtered.fastq\n";
	}
};

class CUsageDisplayerFactory
{
	std::unique_ptr<CUsageDisplayer> desc;
public:
	CUsageDisplayerFactory(CConfig::Mode mode)
	{
		//TODO: add desc of compare, check and info. Check also if others are not missing
		switch (mode)
		{
		case CConfig::Mode::UNDEFINED:
			desc = std::make_unique<CGeneralUsageDisplayer>();
			break;
		case CConfig::Mode::COMPLEX:
			desc = std::make_unique<CComplexUsageDisplayer>();
			break;
		case CConfig::Mode::COMPARE:
			desc = std::make_unique<CGeneralUsageDisplayer>();
			break;
		case CConfig::Mode::FILTER:
			desc = std::make_unique<CFilterUsageDisplayer>();
			break;
		case CConfig::Mode::SIMPLE_SET:
			desc = std::make_unique<CSimpleOperationUsageDisplayer>();
			break;
		case CConfig::Mode::TRANSFORM:
			desc = std::make_unique<CTransformOperationUsageDisplayer>();
			break;
		default:
			desc = std::make_unique<CGeneralUsageDisplayer>();
			break;
		}
	}
	const CUsageDisplayer& GetUsageDisplayer()
	{
		return *desc;
	}
};

#endif


// ***** EOF
