/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#ifndef _CONFIG_H
#define _CONFIG_H

#include "defs.h"
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include "kmc_header.h"
#include "percent_progress.h"

enum class CounterOpType { MIN, MAX, SUM, DIFF, FROM_DB1, FROM_DB2, NONE };

class CCounterBuilder
{
public:
	static void build_counter(uint32& counter, uchar *&buffer, uint32 counter_bytes, uint32 counter_mask, bool little_endian)
	{
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
	uint32 cutoff_max = 0; //0 means it is not set yet
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
	enum class OpType { HISTOGRAM, DUMP, SORT, REDUCE, COMPACT };
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
	uint32 bzip2_buffer_size = 64 << 20;

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

struct CDumpParams
{
	bool sorted_output = false;
};




//************************************************************************************************************
// CConfig - configuration of current application run. Singleton class.
//************************************************************************************************************
class CConfig
{
public:	
	enum class Mode { UNDEFINED, INTERSECTION, KMERS_SUBTRACT, COUNTERS_SUBTRACT, UNION, COMPLEX, SORT, REDUCE, COMPACT, HISTOGRAM, DUMP, COMPARE, FILTER, SIMPLE_SET, TRANSFORM, INFO, CHECK }; 	
	uint32 avaiable_threads;
	uint32 kmer_len = 0;
	Mode mode = Mode::UNDEFINED;
	bool verbose = false;	
	std::vector<CInputDesc> input_desc;
	std::vector<CKMC_header> headers;
		
	COutputDesc output_desc;

	std::vector<CSimpleOutputDesc> simple_output_desc;

	CFilteringParams filtering_params; //for filter operation only	
	CDumpParams dump_params; //for dump operation only
	CCheckParams check_params; // for check operation only
	CounterOpType counter_op_type = CounterOpType::NONE; //for INTERSECTION, COUNTERS_SUBTRACT and UNION only


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

	bool Is2ArgOper()
	{
		return mode == Mode::UNION || mode == Mode::KMERS_SUBTRACT || mode == Mode::COUNTERS_SUBTRACT || mode == Mode::INTERSECTION || mode == Mode::COMPARE;
	}

	bool IsComplex()
	{
		return mode == Mode::COMPLEX;
	}
	
	bool Is1ArgOper()
	{
		return mode == Mode::SORT || mode == Mode::REDUCE || mode == Mode::COMPACT || mode == Mode::HISTOGRAM || mode == Mode::DUMP;		
	}

	std::string GetOperationName()
	{
		switch (mode)
		{
		case CConfig::Mode::UNDEFINED:
			return "";			
		case CConfig::Mode::INTERSECTION:
			return "intersect";			
		case CConfig::Mode::KMERS_SUBTRACT:
			return "kmers_subtract";
		case CConfig::Mode::COUNTERS_SUBTRACT:
			return "counters_subtract";
		case CConfig::Mode::UNION:
			return "union";			
		case CConfig::Mode::COMPLEX:
			return "complex";			
		case CConfig::Mode::SORT:
			return "sort";
		case CConfig::Mode::REDUCE:
			return "reduce";
		case CConfig::Mode::COMPACT:
			return "compact";
		case CConfig::Mode::HISTOGRAM:
			return "histogram";
		case CConfig::Mode::DUMP:
			return "dump";
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
		std::cout << "The '" << name << "' is two arguments' operation. General syntax:\n";
		std::cout << " kmc_tools " << name << " <input1 [input1_params]> <input2 [input2_params]> <output [output_params]>\n";
		std::cout << " input1, input2 - paths to databases generated by KMC \n";
		std::cout << " output		  - path  to output database\n";
		std::cout << " For each input there are additional parameters:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
		std::cout << " For output there are additional parameters:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
		std::cout << "  -cs<value> - maximal value of a counter\n";
	}

	void Display1ArgGeneral(bool output_params) const
	{
		std::cout << " The '" << name << "' is one argument operation. General syntax:\n";
		std::cout << " kmc_tools " << name << " <input> [input_params] <output> "<< (output_params ? "[output_params]" : "") << "\n";
		std::cout << " input - path to database generated by KMC \n";
		std::cout << " For input there are additional parameters:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
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
	//void Display() const override //old
	//{
	//	std::cout << "kmc_tools ver. " << KMC_VER << " (" << KMC_DATE << ")\n";
	//	std::cout << "Usage:\n kmc_tools [global parameters] <operation> [operation parameters]\n";		
	//	std::cout << "Available operations:\n";
	//	std::cout << " k-mers sets' operations for 2 KMC's databases:\n";
	//	std::cout << "  intersect           - intersection of 2 k-mers' sets\n";
	//	std::cout << "  kmers_subtract      - subtraction of 2 k-mers' sets\n";
	//	std::cout << "  counters_subtract   - counters' subtraction of 2 k-mers' sets\n";
	//	std::cout << "  union               - union of 2 k-mers' sets\n\n";
	//	std::cout << " operations for single kmc database:\n";
	//	std::cout << "  sort                - sorts k-mers from database generated by KMC2.x\n";
	//	std::cout << "  reduce              - exclude too rare and too frequent k-mers\n";
	//	std::cout << "  compact             - remove counters (store only k-mers)\n";
	//	std::cout << "  histogram           - histogram of k-mers occurrences\n";
	//	std::cout << "  dump                - dump k-mers and counters to text file\n";		
	//	std::cout << " more complex operations:\n";
	//	std::cout << "  complex             - complex operations with a number of input databases\n";
	//	std::cout << " other operatations:\n";
	//	std::cout << "  filter              - filter out reads with too small number of k-mers\n";
	//	std::cout << " global parameters:\n";
	//	std::cout << "  -t<value>           - total number of threads (default: no. of CPU cores)\n";
	//	std::cout << "  -v                  - enable verbose mode (shows some information) (default: false)\n";
	//	std::cout << "  -hp                 - hide percentage progress (default: false)\n";
	//	std::cout << "Example:\n";
	//	std::cout << "kmc_tools union db1 -ci3 db2 -ci5 -cx300 db1_union_db2 -ci10\n";
	//	std::cout << "For detailed help of concrete operation type operation name without parameters:\n";
	//	std::cout << "kmc_tools union\n";
	//}

	void Display() const override
	{
		std::cout << "kmc_tools ver. " << KMC_VER << " (" << KMC_DATE << ")\n";
		std::cout << "Usage:\n kmc_tools [global parameters] <operation> [operation parameters]\n";
		std::cout << "Available operations:\n";
		std::cout << "  transform            - transforms single KMC's database\n";
		std::cout << "  simple               - performs set operation on two KMC's databases\n";
		std::cout << "  complex              - performs set operation on multiple KMC's databases\n";
		std::cout << "  filter               - filter out reads with too small number of k-mers\n";
		std::cout << " global parameters:\n";
		std::cout << "  -t<value>            - total number of threads (default: no. of CPU cores)\n";
		std::cout << "  -v                   - enable verbose mode (shows some information) (default: false)\n";
		std::cout << "  -hp                  - hide percentage progress (default: false)\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools simple db1 -ci3 db2 -ci5 -cx300 union db1_union_db2 -ci10\n";
		std::cout << "For detailed help of concrete operation type operation name without parameters:\n";
		std::cout << "kmc_tools simple\n";
	}
};

class CSimpleOperationUsageDisplayer : public CUsageDisplayer
{
public:
	CSimpleOperationUsageDisplayer() : CUsageDisplayer("simple"){}
	void Display() const override
	{
		std::cout << "simple operation performs set operation on two input KMC's databases\n";
		std::cout << "General syntax:\n";
		std::cout << "kmc_tools simple <input1 [input1_params]> <input2 [input2_params]> <oper1 output1 [output_params1]> [<oper2 output2 [output_params2]> ...]\n";
		std::cout << "input1, input2                    - paths to databases generated by KMC\n";
		std::cout << "oper1, oper2, ..., operN          - set operations to be performed on input1 and input2\n";
		std::cout << "output1, output2, ..., outputN    - paths to output k-mer databases\n";
		std::cout << " Available operations:\n";
		std::cout << "  intersect                  - output database will contains only k-mers that are present in both input sets\n";
		std::cout << "  union                      - output database will contains each k-mer present in any of input sets\n";
		std::cout << "  kmers_subtract             - difference of input sets based on k-mers. \n" <<
					 "                               Output database will contains only k-mers that are present in first input set but absent in the second one\n";
		std::cout << "  counters_subtract          - difference of input sets based on k-mers and their counters (weaker version of kmers_subtract).\n" <<
					 "                               Output database will contains all k-mers that are present in first input, \n" <<
					 "                               beyond those for which counter operation will lead to remove (i.e. counter equal to 0 or negative number)\n";	
		std::cout << "  reverse_kmers_subtract     - same as kmers_subtract but treat input2 as first and input1 as second\n";
		std::cout << "  reverse_counters_subtract  - same as counters_subtract but treat input2 as first and input1 as second\n";
		std::cout << " For each input there are additional parameters:\n";
		std::cout << "  -ci<value>  - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value>  - exclude k-mers occurring more of than <value> times\n";
		std::cout << " For each output there are additional parameters:\n";
		std::cout << "  -ci<value>  - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value>  - exclude k-mers occurring more of than <value> times\n";
		std::cout << "  -cs<value>  - maximal value of a counter\n";
		std::cout << "  -oc<value>  - redefine counter calculation mode for equal k-mers\n";
		std::cout << "    Available values : \n";
		std::cout << "      min   - get lower value of a k-mer counter (default value for intersect operation)\n";
		std::cout << "      max   - get upper value of a k-mer counter\n";
		std::cout << "      sum   - get sum of counters from both databases\n";
		std::cout << "      diff  - get difference between counters (default for counters_subtract operation)\n";
		std::cout << "      left  - get counter from first database (input1)\n";
		std::cout << "      right - get counter from second database (input2)\n";		
		std::cout << "Example:\n";
		std::cout << "kmc_tools simple kmers1 -ci3 -cx70000 kmers2 union kmers1_kmers2_union -cs65536 -ocfirst intersect intersect_kmers1_kmers2 intersect intersect_max_kmers1_kmers2 -ocmax\n";
	}
};


class CTransformOperationUsageDisplayer : public CUsageDisplayer
{
public:
	CTransformOperationUsageDisplayer() : CUsageDisplayer("transform"){}
	void Display() const override
	{
		std::cout << "transform operation transforms single input database to output (text file or KMC database)\n";
		std::cout << "General syntax:\n";
		std::cout << "kmc_tools transform <input> [input_params] <oper1 [oper_params1] output1 [output_params1]> [<oper2 [oper_params2] output2 [output_params2]>...]\n";
		std::cout << "input - path to database generated by KMC \n";
		std::cout << "oper1, oper2, ..., operN          - transform operation name\n";
		std::cout << "output1, output2, ..., outputN    - paths to output\n";

		std::cout << " Available operations:\n";
		std::cout << "  sort                       - converts database produced by KMC2.x to KMC1.x database format (which contains k-mers in sorted order)\n";
		std::cout << "  reduce                     - exclude too rare and too frequent k-mers\n";
		std::cout << "  compact                    - remove counters of k-mers\n";
		std::cout << "  histogram                  - produce histogram of k-mers occurrences\n";
		std::cout << "  dump                       - produce text dump of kmc database\n";

		std::cout << " For input there are additional parameters:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";

		std::cout << " For sort and reduce operations there are additional output_params:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
		std::cout << "  -cs<value> - maximal value of a counter\n";

		std::cout << " For dump operation there are additional oper_params:\n";
		std::cout << "  -s - sorted output\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools transform db reduce err_kmers -cx10 reduce valid_kmers -ci11 histogram histo.txt dump dump.txt\n";
	}
};
class CUnionUsageDisplayer : public CUsageDisplayer
{
	public:
		CUnionUsageDisplayer() :CUsageDisplayer("union")
		{}
		void Display() const override
		{
			Display2ArgGeneral();
			std::cout << "The output database will contains each k-mer present in both input sets. For the same k-mers in first and second input the counter in output is equal to sum from inputs.";
			std::cout << "Example:\n";
			std::cout << "kmc - k28 file1.fastq kmers1 tmp\n";
			std::cout << "kmc - k28 file2.fastq kmers2 tmp\n";
			std::cout << "kmc_tools union kmers1 -ci3 -cx70000 kmers2 kmers1_kmers2_union -cs65536\n";
		}
};

class CIntersectUsageDisplayer : public CUsageDisplayer
{
public:
	CIntersectUsageDisplayer() :CUsageDisplayer("intersect")
	{}
	void Display() const override
	{
		Display2ArgGeneral();
		std::cout << "The output database will contains only k-mers that are present in both input sets. The counter value in output database is equal to lower counter value in input.";
		std::cout << "Example:\n";
		std::cout << "kmc - k28 file1.fastq kmers1 tmp\n";
		std::cout << "kmc - k28 file2.fastq kmers2 tmp\n";
		std::cout << "kmc_tools intersect kmers1 -ci10 -cx200 kmers2 -ci4 -cx100 kmers1_kmers2_intersect -ci20 -cx150\n";
	}
};


class CCountersSubtractUsageDisplayer : public CUsageDisplayer
{
public:
	CCountersSubtractUsageDisplayer() :CUsageDisplayer("counters_subtract")
	{}
	void Display() const override
	{
		Display2ArgGeneral();
		std::cout << "The output database will contains only k-mers that are present in first input set and have counters higher than apropriate k - mers in second set. For each k - mer the counter is equal to difference between counter in first set and counter in second set.";
		std::cout << "Example:\n";
		std::cout << "kmc -k28 file1.fastq kmers1 tmp\n";
		std::cout << "kmc -k28 file2.fastq kmers2 tmp\n";
		std::cout << "kmc_tools counters_subtract kmers1 kmers2 kmers1_kmers2_counters_subtract\n";
	}
};

class CKmersSubtractUsageDisplayer : public CUsageDisplayer
{
public:
	CKmersSubtractUsageDisplayer() :CUsageDisplayer("kmers_subtract")
	{}
	void Display() const override
	{
		Display2ArgGeneral();
		std::cout << "The output database will contains only k-mers that are present in first input set but absent in the second one. The counter value is equal to value from first input set.";
		std::cout << "Example:\n";
		std::cout << "kmc - k28 file1.fastq kmers1 tmp\n";
		std::cout << "kmc - k28 file2.fastq kmers2 tmp\n";
		std::cout << "kmc_tools kmers_subtract kmers1 kmers2 kmers1_kmers2_subtract - cs200\n";
	}
};

class CComplexUsageDisplayer : public CUsageDisplayer
{
public:
	CComplexUsageDisplayer() :CUsageDisplayer("complex")
	{}
	void Display() const override
	{
		std::cout << "Complex operation allows one to define operations for more than 2 input k-mers sets. Command-line syntax:\n";
		std::cout << "kmc_tools complex <operations_definition_file>\n";
		std::cout << " operations_definition_file - path to file which define input sets and operations. It is text file with following syntax:\n";
		std::cout << " ______________________________________________________________________________ \n";
		std::cout << "|INPUT:                                                                        |\n";
		std::cout << "|<input1>=<input1_db_path> [params]                                            |\n";
		std::cout << "|<input2>=<input2_db_path> [params]                                            |\n";
		std::cout << "|...                                                                           |\n";
		std::cout << "|<inputN>=<inputN_db_path> [params]                                            |\n";
		std::cout << "|OUTPUT:                                                                       |\n";
		std::cout << "|<out_db>=<ref_input><oper[c_mode]><ref_input>[<oper[c_mode]><ref_input>[...]  |\n";
		std::cout << "|[OUTPUT_PARAMS:                                                             __|\n";
		std::cout << "|<output_params>]                                                           |  /\n";
		std::cout << "|                                                                           | / \n";
		std::cout << "|___________________________________________________________________________|/  \n";
		std::cout << "input1, input2, ..., inputN - names of inputs used to define equation\n";
		std::cout << "input1_db, input2_db_path, ..., inputN_db_path - paths to k-mers sets\n";
		std::cout << "For each input there are additional parameters which can be set:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
		std::cout << "out_db_path - path to output database\n";
		std::cout << "ref_input - one of input1, input2, ..., inputN\n";
		std::cout << "oper - one of {*,-,~,+}, which refers to {intersect, kmers_subtract, counters_subtract, union}\n";
		std::cout << "operator * has the highest priority. Other operators has equal priorities. Order of operations can be changed with parentheses\n";
		std::cout << "for {*,~,+} it is possible to redefine counter calculation mode ([c_mode]). Available values: min, max, diff, sum, left, right (detailet description available in simple help message)\n";
		std::cout << "output_params are:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
		std::cout << "  -cs<value> - maximal value of a counter\n";

		std::cout << "Example:\n";
		std::cout << " __________________________________________________________________ \n";
		std::cout << "|INPUT:                                                            |\n";
		std::cout << "|set1 = kmc_o1 -ci5                                                |\n";
		std::cout << "|set2 = kmc_o2                                                     |\n";
		std::cout << "|set3 = kmc_o3 -ci10 -cx100                                      __|\n";		
		std::cout << "|OUTPUT:                                                        |  /\n";
		std::cout << "|result = (set3 + min set1) * right set2                        | / \n";
		std::cout << "|_______________________________________________________________|/  \n";
		
	}
};

class CSortUsageDisplayer : public CUsageDisplayer
{
public:
	CSortUsageDisplayer() :CUsageDisplayer("sort")
	{}
	void Display() const override
	{
		Display1ArgGeneral(true);
		std::cout << " For output there are additional parameters:\n";
		std::cout << "  -cs<value> - maximal value of a counter\n";
		std::cout << "Converts database produced by KMC2.x to KMC1.x database format (which contains k-mers in sorted order)\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools sort wy_kmc2 -ci3 -cx1000 wy_kmc1 -cs255\n";
	}
};

class CReduceUsageDisplayer : public CUsageDisplayer
{
public:
	CReduceUsageDisplayer() :CUsageDisplayer("reduce")
	{}
	void Display() const override
	{
		Display1ArgGeneral(true);
		std::cout << " For output there are additional parameters:\n";
		std::cout << "  -cs<value> - maximal value of a counter\n";
		std::cout << "Exclude too rare and too frequent k-mers\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools reduce wy_kmc2 -ci3 -cx1000 wy_kmc1 -cs255\n";
	}
};


class CCompactUsageDisplayer : public CUsageDisplayer
{
public:
	CCompactUsageDisplayer() :CUsageDisplayer("compact")
	{}
	void Display() const override
	{
		Display1ArgGeneral(false);		
		std::cout << "Remove counters of k-mers\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools compact wy_kmc2 -ci3 -cx1000 wy_kmc1\n";
	}
};

class CHistogramUsageDisplayer : public CUsageDisplayer
{
public:
	CHistogramUsageDisplayer() :CUsageDisplayer("histogram")
	{}
	void Display() const override
	{
		Display1ArgGeneral(false);
		std::cout << "Produce histogram of k-mers occurrences\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools histogram wy_kmc2 -ci3 -cx1000 histo.txt\n";
	}
};


class CDumpUsageDisplayer : public CUsageDisplayer
{
public:
	CDumpUsageDisplayer() :CUsageDisplayer("dump")
	{}
	void Display() const override
	{
		std::cout << " The '" << name << "' is one argument operation. General syntax:\n";
		std::cout << " kmc_tools " << name << " [dump_params] <input> [input_params] <output>\n";
		std::cout << " dump_params:\n";
		std::cout << "  -s - sorted output\n";
		std::cout << " input - path to database generated by KMC \n";
		std::cout << " For input there are additional parameters:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";


		std::cout << "Produce text dump of kmc database\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools dump wy_kmc2 -ci3 -cx1000 dump.txt\n";
	}
};

class CFilterUsageDisplayer : public CUsageDisplayer
{
public:
	CFilterUsageDisplayer() : CUsageDisplayer("filter")
	{}
	void Display() const override
	{
		std::cout << " The '" << name << "' is two arguments' operation. General syntax:\n";
		std::cout << " kmc_tools " << name << " [filter_params] <kmc_input_db> [kmc_input_db_params] <input_read_set> [input_read_set_params] <output_read_set> [output_read_set_params]\n";
		std::cout << " filter_params:\n";
		std::cout << " -t               - trim reads on first invalid k-mer instead of remove entirely\n";
		std::cout << " -hm              - hard mask invalid k-mers in a read\n";
		std::cout << " kmc_input_db     - path to database generated by KMC \n";
		std::cout << " input_read_set   - path to input set of reads \n";
		std::cout << " output_read_set  - path to output set of reads \n";
		std::cout << " For k-mers' database there are additional parameters:\n";
		std::cout << "  -ci<value> - exclude k-mers occurring less than <value> times \n";
		std::cout << "  -cx<value> - exclude k-mers occurring more of than <value> times\n";
		std::cout << " For input set of reads there are additional parameters:\n";
		std::cout << " -ci<value> - remove reads containing less k-mers than value. It can be integer or floating number in range [0.0;1.0] (default: 2)\n";
		std::cout << " -cx<value> - remove reads containing more k-mers than value. It can be integer or floating number in range [0.0;1.0] (default: 1e9)\n";
		std::cout << " -f<a/q>    - input in FASTA format (-fa), FASTQ format (-fq); default: FASTQ\n";
		std::cout << " For output set of reads there are additional parameters:\n";
		std::cout << " -f<a/q>    - output in FASTA format (-fa), FASTQ format (-fq); default: same as input\n";
		std::cout << "Example:\n";
		std::cout << "kmc_tools filter kmc_db -ci3 input.fastq -ci0.5 -cx1.0 filtered.fastq\n";
		std::cout << "kmc_tools filter kmc_db input.fastq -ci10 -cx100 filtered.fastq\n";
	}
};

class CUsageDisplayerFactory
{
	std::unique_ptr<CUsageDisplayer> desc;
public:
	CUsageDisplayerFactory(CConfig::Mode mode)
	{
		switch (mode)
		{
		case CConfig::Mode::UNDEFINED:
			desc = std::make_unique<CGeneralUsageDisplayer>();
			break;
		case CConfig::Mode::INTERSECTION:
			desc = std::make_unique<CIntersectUsageDisplayer>();
			break;
		case CConfig::Mode::KMERS_SUBTRACT:
			desc = std::make_unique<CKmersSubtractUsageDisplayer>();
			break;
		case CConfig::Mode::COUNTERS_SUBTRACT:
			desc = std::make_unique<CCountersSubtractUsageDisplayer>();
			break;
		case CConfig::Mode::UNION:
			desc = std::make_unique<CUnionUsageDisplayer>();
			break;
		case CConfig::Mode::COMPLEX:
			desc = std::make_unique<CComplexUsageDisplayer>();
			break;
		case CConfig::Mode::SORT:
			desc = std::make_unique<CSortUsageDisplayer>();
			break;
		case CConfig::Mode::REDUCE:
			desc = std::make_unique<CReduceUsageDisplayer>();
			break;
		case CConfig::Mode::COMPACT:
			desc = std::make_unique<CCompactUsageDisplayer>();
			break;
		case CConfig::Mode::HISTOGRAM:
			desc = std::make_unique<CHistogramUsageDisplayer>();
			break;
		case CConfig::Mode::DUMP:
			desc = std::make_unique<CDumpUsageDisplayer>();
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
