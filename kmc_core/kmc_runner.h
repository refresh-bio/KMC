/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.2
  Date   : 2023-03-10
*/

#ifndef _KMC_RUNNER
#define _KMC_RUNNER
#include <cinttypes>
#include <vector>
#include <string>
#include <memory>
#include <thread>


#define DEVELOP_MODE

namespace KMC
{
	class IPercentProgressObserver
	{
	public:
		virtual void SetLabel(const std::string& label) = 0;
		virtual void ProgressChanged(int newValue) = 0;
		virtual ~IPercentProgressObserver() = default;
	};
	
	//In cases where KMC does not estimate percentage, just informs that some progress is made
	class IProgressObserver
	{
	public:
		virtual void Start(const std::string& name) = 0;
		virtual void Step() = 0;
		virtual void End() = 0;
		virtual ~IProgressObserver() = default;
	};

	class ILogger
	{
	public:
		virtual void Log(const std::string& msg) = 0;
		virtual ~ILogger() = default;
	};

	class CerrPercentProgressObserver : public IPercentProgressObserver
	{
		std::string label;
	public:
		void SetLabel(const std::string& label) override;
		void ProgressChanged(int newValue) override;
	};

	class NullPercentProgressObserver : public IPercentProgressObserver
	{
	public:
		void SetLabel(const std::string& label) override;
		void ProgressChanged(int newValue) override;
	};

	class CerrProgressObserver : public IProgressObserver
	{
	public:
		void Start(const std::string& name) override;
		void Step() override;
		void End() override;
	};

	class NullProgressObserver : public IProgressObserver
	{
	public:
		void Start(const std::string& name) override;
		void Step() override;
		void End() override;
	};

	class NullLogger : public ILogger
	{
		void Log(const std::string& msg) override;
	};

	class CerrVerboseLogger : public ILogger
	{
		void Log(const std::string& msg) override;
	};

	class CerrWarningLogger: public ILogger
	{
		void Log(const std::string& msg) override;
	};

	enum class InputFileType { FASTQ, FASTA, MULTILINE_FASTA, BAM, KMC };
	enum class OutputFileType { KMC, KFF };
	
	enum class EstimateHistogramCfg { DONT_ESTIMATE, ESTIMATE_AND_COUNT_KMERS, ONLY_ESTIMATE };

	class Stage1Params
	{
		struct 
		{			
			NullLogger defaultVerboseLogger;
			CerrPercentProgressObserver defaultPercentProgressObserver;
			CerrWarningLogger defaultWarningsLogger;
			CerrProgressObserver defaultProgressObserver;
		} defaults;
		
		
		std::vector<std::string> inputFiles;
		std::string tmpPath = ".";
		uint32_t kmerLen = 25;
		uint32_t nThreads = std::thread::hardware_concurrency();
		uint32_t maxRamGB = 12;
		uint32_t signatureLen = 9;		
		bool homopolymerCompressed = false;
		InputFileType inputFileType = InputFileType::FASTQ;
		bool canonicalKmers = true;
		bool ramOnlyMode = false;
		uint32_t nBins = 512;
		uint32_t nReaders = 0;
		uint32_t nSplitters = 0;
		ILogger* verboseLogger = &defaults.defaultVerboseLogger;
		IPercentProgressObserver* percentProgressObserver = &defaults.defaultPercentProgressObserver;
		ILogger* warningsLogger = &defaults.defaultWarningsLogger;
		EstimateHistogramCfg estimateHistogramCfg = EstimateHistogramCfg::DONT_ESTIMATE;
		IProgressObserver* progressObserver = &defaults.defaultProgressObserver;
#ifdef DEVELOP_MODE
		bool developVerbose = false;
#endif
	public:		
		Stage1Params& SetInputFiles(const std::vector<std::string>& inputFiles);
		Stage1Params& SetTmpPath(const std::string& tmpPath);
		Stage1Params& SetKmerLen(uint32_t kmerLen);
		Stage1Params& SetNThreads(uint32_t nThreads);
		Stage1Params& SetMaxRamGB(uint32_t maxRamGB);
		Stage1Params& SetSignatureLen(uint32_t signatureLen);		
		Stage1Params& SetHomopolymerCompressed(bool homopolymerCompressed);
		Stage1Params& SetInputFileType(InputFileType inputFileType);
		Stage1Params& SetCanonicalKmers(bool canonicalKmers);
		Stage1Params& SetRamOnlyMode(bool ramOnlyMode);
		Stage1Params& SetNBins(uint32_t nBins);
		Stage1Params& SetNReaders(uint32_t nReaders);
		Stage1Params& SetNSplitters(uint32_t nSplitters);
		Stage1Params& SetVerboseLogger(ILogger* verboseLogger);
		Stage1Params& SetPercentProgressObserver(IPercentProgressObserver* percentProgressObserver);
		Stage1Params& SetWarningsLogger(ILogger* warningsLogger);
		Stage1Params& SetEstimateHistogramCfg(EstimateHistogramCfg estimateHistogramCfg);
		Stage1Params& SetProgressObserver(IProgressObserver* progressObserver);
#ifdef DEVELOP_MODE
		Stage1Params& SetDevelopVerbose(bool developVerbose);
#endif

		const std::vector<std::string>& GetInputFiles() const noexcept { return inputFiles; }
		const std::string& GetTmpPath() const noexcept { return tmpPath; }
		uint32_t  GetKmerLen() const noexcept { return kmerLen; }
		uint32_t GetNThreads() const noexcept { return nThreads; }
		uint32_t GetMaxRamGB() const noexcept { return maxRamGB; }
		uint32_t GetSignatureLen() const noexcept { return signatureLen; }
		bool GetHomopolymerCompressed() const noexcept { return homopolymerCompressed; }
		InputFileType GetInputFileType() const noexcept { return inputFileType; }
		bool GetCanonicalKmers() const noexcept { return canonicalKmers; }
		bool GetRamOnlyMode() const noexcept { return ramOnlyMode; }
		uint32_t GetNBins() const noexcept { return nBins; }
		uint32_t GetNReaders() const noexcept { return nReaders; }
		uint32_t GetNSplitters() const noexcept { return nSplitters; }
		ILogger* GetVerboseLogger() const noexcept { return verboseLogger; }
		IPercentProgressObserver* GetPercentProgressObserver() const noexcept { return percentProgressObserver; }
		ILogger* GetWarningsLogger() const noexcept { return warningsLogger; }
		EstimateHistogramCfg GetEstimateHistogramCfg() const noexcept { return estimateHistogramCfg; }
		IProgressObserver* GetProgressObserver() const noexcept { return progressObserver; }
#ifdef DEVELOP_MODE
		bool GetDevelopVerbose() const noexcept { return developVerbose; }
#endif
	};


	class Stage2Params
	{	
		uint32_t maxRamGB = 12;
		uint32_t nThreads = std::thread::hardware_concurrency();
		bool strictMemoryMode = false;
		uint64_t cutoffMin = 2;
		uint64_t counterMax = 255;
		uint64_t cutoffMax = 1000000000;		
		std::string outputFileName;
		OutputFileType outputFileType = OutputFileType::KMC;
		bool withoutOutput = false;
		uint32_t strictMemoryNSortingThreadsPerSorters = 0;
		uint32_t strictMemoryNUncompactors = 0;
		uint32_t strictMemoryNMergers = 0;

		
	public:	
		Stage2Params& SetMaxRamGB(uint32_t maxRamGB);
		Stage2Params& SetNThreads(uint32_t nThreads);
		Stage2Params& SetStrictMemoryMode(bool strictMemoryMode);
		Stage2Params& SetCutoffMin(uint64_t cutoffMin);
		Stage2Params& SetCounterMax(uint64_t counterMax);
		Stage2Params& SetCutoffMax(uint64_t cutoffMax);		
		Stage2Params& SetOutputFileName(const std::string& outputFileName);
		Stage2Params& SetOutputFileType(OutputFileType outputFileType);
		Stage2Params& SetWithoutOutput(bool withoutOutput);		
		Stage2Params& SetStrictMemoryNSortingThreadsPerSorters(uint32_t strictMemoryNSortingThreadsPerSorters);
		Stage2Params& SetStrictMemoryNUncompactors(uint32_t strictMemoryNUncompactors);
		Stage2Params& SetStrictMemoryNMergers(uint32_t strictMemoryNMergers);

		uint32_t GetMaxRamGB() const noexcept { return maxRamGB; }
		uint32_t GetNThreads() const noexcept { return nThreads; }
		bool GetStrictMemoryMode() const noexcept { return strictMemoryMode; }
		uint64_t GetCutoffMin() const noexcept { return cutoffMin; }
		uint64_t GetCounterMax() const noexcept { return counterMax; }
		uint64_t GetCutoffMax() const noexcept { return cutoffMax; }
		const std::string& GetOutputFileName() const noexcept { return outputFileName; }
		OutputFileType GetOutputFileType() const noexcept { return outputFileType; }
		bool GetWithoutOutput() const noexcept { return withoutOutput; }
		uint32_t GetStrictMemoryNSortingThreadsPerSorters() const noexcept { return strictMemoryNSortingThreadsPerSorters; }
		uint32_t GetStrictMemoryNUncompactors() const noexcept { return strictMemoryNUncompactors; }
		uint32_t GetStrictMemoryNMergers() const noexcept { return strictMemoryNMergers; }
	};

	struct Stage1Results
	{
		double time{};
		uint64_t nSeqences{};
		bool wasSmallKOptUsed = false;
		uint64_t nTotalSuperKmers{};
		uint64_t tmpSize{};
		std::vector<uint64_t> estimatedHistogram;
	};

	struct Stage2Results
	{
		double time{};
		double timeStrictMem{};		
		uint64_t tmpSizeStrictMemory{};
		uint64_t maxDiskUsage{};
		uint64_t nBelowCutoffMin{};
		uint64_t nAboveCutoffMax{};
		uint64_t nTotalKmers{}; //TODO: this can be get after first stage, maybe changed
		uint64_t nUniqueKmers{};
	};

	
	class Runner
	{
		class RunnerImpl;
		const std::unique_ptr<RunnerImpl> pImpl;
	public:
		Runner();
		~Runner();
		Stage1Results RunStage1(const Stage1Params& params);
		Stage2Results RunStage2(const Stage2Params& params);
	};

	struct CfgConsts
	{
		static const std::string kmc_ver;
		static const std::string kmc_date;
		static const uint32_t min_k;
		static const uint32_t max_k;
		//static const uint32_t min_mem;
		//static const uint32_t min_sl;
		//static const uint32_t max_sl;
	};
	


}
#endif // !_KMC_RUNNER
