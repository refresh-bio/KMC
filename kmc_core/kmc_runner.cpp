/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.4
  Date   : 2024-02-09
*/

#include "kmc_runner.h"
#include "kmc.h"
#include "defs.h"

#include <string>
#include <exception>
#include <memory>

namespace KMC
{	
	//----------------------------------------------------------------------------------
// Application class
// Template parameters:
//    * SIZE     - maximal size of the k-mer (divided by 32)
	template<unsigned SIZE> class CApplication
	{
		std::unique_ptr<CApplication<SIZE - 1>> app_1;
		std::unique_ptr<CKMC<SIZE>> kmc;
		int p_k;
		bool is_selected;

	public:

		CApplication(uint32_t p_k)
		{
			this->p_k = p_k;
			is_selected = p_k <= (int32)SIZE * 32 && p_k > ((int32)SIZE - 1) * 32;

			app_1 = std::make_unique<CApplication<SIZE - 1>>(p_k);
			if (is_selected)
			{
				kmc = std::make_unique<CKMC<SIZE>>();
			}
		}

		KMC::Stage1Results ProcessStage1(const KMC::Stage1Params& stage1Params)
		{
			if (is_selected)
			{
				kmc->SetParamsStage1(stage1Params);
				return kmc->ProcessStage1();
			}
			else
				return app_1->ProcessStage1(stage1Params);
		}

		KMC::Stage2Results ProcessStage2(const KMC::Stage2Params& stage2Params)
		{
			if (is_selected)
			{
				kmc->SetParamsStage2(stage2Params);
				return kmc->ProcessStage2();
			}
			else
				return app_1->ProcessStage2(stage2Params);
		}
	};

	//----------------------------------------------------------------------------------
	// Specialization of the application class for the SIZE=1
	template<> class CApplication<1>
	{
		std::unique_ptr<CKMC<1>> kmc;
		int p_k;
		bool is_selected;

	public:
		CApplication(uint32_t p_k) {
			this->p_k = p_k;
			is_selected = p_k <= 32;
			if (is_selected)
			{
				kmc = std::make_unique<CKMC<1>>();
			}
		};

		KMC::Stage1Results ProcessStage1(const KMC::Stage1Params& stage1Params)
		{
			if (is_selected)
			{
				kmc->SetParamsStage1(stage1Params);
				return kmc->ProcessStage1();
			}
			else
				throw std::runtime_error("Running stage 1 failed");
		}

		KMC::Stage2Results ProcessStage2(const KMC::Stage2Params& stage2Params)
		{
			if (is_selected)
			{
				kmc->SetParamsStage2(stage2Params);
				return kmc->ProcessStage2();
			}
			else
				throw std::runtime_error("Running stage 2 failed");
		}
	};

	void CerrPercentProgressObserver::SetLabel(const std::string& label)
	{
		this->label = label;		
	}
	void CerrPercentProgressObserver::ProgressChanged(int newValue)
	{
		std::cerr << "\r" << label << newValue << "%";
		if (newValue == 100)
			std::cerr << "\n";
		std::cerr.flush();
	}

	void NullPercentProgressObserver::SetLabel(const std::string& label)
	{

	}

	void NullPercentProgressObserver::ProgressChanged(int newValue)
	{
		
	}

	void CerrProgressObserver::Start(const std::string& /*name*/)
	{

	}

	void CerrProgressObserver::Step()
	{
		std::cerr << '*';
	}

	void CerrProgressObserver::End()
	{
		std::cerr << '\n';
	}

	void NullProgressObserver::Start(const std::string& /*name*/)
	{

	}

	void NullProgressObserver::Step()
	{

	}

	void NullProgressObserver::End()
	{

	}



	void NullLogger::Log(const std::string& msg)
	{

	}

	void CerrVerboseLogger::Log(const std::string& msg)
	{
		std::cerr << msg;
	}

	void CerrWarningLogger::Log(const std::string& msg)
	{
		std::cerr << "Warning: " << msg << "\n";
	}

	Stage1Params& Stage1Params::SetInputFiles(const std::vector<std::string>& inputFiles)
	{
		this->inputFiles = inputFiles;
		return *this;
	}
	Stage1Params& Stage1Params::SetTmpPath(const std::string& tmpPath)
	{
		this->tmpPath = tmpPath;
		return *this;
	}
	Stage1Params& Stage1Params::SetKmerLen(uint32_t kmerLen)
	{
		if (kmerLen < MIN_K || kmerLen > MAX_K)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: k must be from range <" << MIN_K << "," << MAX_K << ">";
			throw std::runtime_error(err_msg.str());			
		}
		this->kmerLen = kmerLen;
		return *this;
	}
	Stage1Params& Stage1Params::SetNThreads(uint32_t nThreads)
	{
		this->nThreads = nThreads;
		return *this;
	}
	Stage1Params& Stage1Params::SetMaxRamGB(uint32_t maxRamGB)
	{
		if (maxRamGB < MIN_MEM)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameret: min memory must be at least " << MIN_MEM << "GB";
			throw std::runtime_error(err_msg.str());
		}
		this->maxRamGB = maxRamGB;
		return *this;
	}
	Stage1Params& Stage1Params::SetSignatureLen(uint32_t signatureLen)
	{
		if (signatureLen < MIN_SL || signatureLen > MAX_SL)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: p must be from range <" << MIN_SL << "," << MAX_SL << ">";
			throw std::runtime_error(err_msg.str());
		}
		this->signatureLen = signatureLen;
		return *this;
	}
	
	Stage1Params& Stage1Params::SetHomopolymerCompressed(bool homopolymerCompressed)
	{
		this->homopolymerCompressed = homopolymerCompressed;
		return *this;
	}
	Stage1Params& Stage1Params::SetInputFileType(InputFileType inputFileType)
	{
		this->inputFileType = inputFileType;
		return *this;
	}
	Stage1Params& Stage1Params::SetCanonicalKmers(bool canonicalKmers)
	{
		this->canonicalKmers = canonicalKmers;
		return *this;
	}
	Stage1Params& Stage1Params::SetRamOnlyMode(bool ramOnlyMode)
	{
		this->ramOnlyMode = ramOnlyMode;
		return *this;
	}
	Stage1Params& Stage1Params::SetNBins(uint32_t nBins)
	{
		if (nBins < MIN_N_BINS || nBins > MAX_N_BINS)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: number of bins must be in range <" << MIN_N_BINS << "," << MAX_N_BINS << ">";
			throw std::runtime_error(err_msg.str());			
		}
		this->nBins = nBins;
		return *this;
	}
	Stage1Params& Stage1Params::SetNReaders(uint32_t nReaders)
	{
		if (nReaders < MIN_SF || nReaders > MAX_SF)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: number of reading thread must be from range <" << MIN_SF << "," << MAX_SF << ">";
			throw std::runtime_error(err_msg.str());			
		}
		this->nReaders = nReaders;
		return *this;
	}
	Stage1Params& Stage1Params::SetNSplitters(uint32_t nSplitters)
	{		
		if (nSplitters < MIN_SP || nSplitters > MAX_SP)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: number of splitting threads must be in range <" << MIN_SP << "," << MAX_SP << ">";
			throw std::runtime_error(err_msg.str());			
		}
		this->nSplitters = nSplitters;
		return *this;
	}
	Stage1Params& Stage1Params::SetVerboseLogger(ILogger* verboseLogger)
	{
		this->verboseLogger = verboseLogger;
		return *this;
	}
	Stage1Params& Stage1Params::SetPercentProgressObserver(IPercentProgressObserver* percentProgressObserver)
	{
		this->percentProgressObserver = percentProgressObserver;
		return *this;
	}
	Stage1Params& Stage1Params::SetWarningsLogger(ILogger* warningsLogger)
	{
		this->warningsLogger = warningsLogger;
		return *this;
	}
	Stage1Params& Stage1Params::SetEstimateHistogramCfg(EstimateHistogramCfg estimateHistogramCfg)
	{
		this->estimateHistogramCfg = estimateHistogramCfg;
		return *this;
	}
	Stage1Params& Stage1Params::SetProgressObserver(IProgressObserver* progressObserver)
	{
		this->progressObserver = progressObserver;
		return *this;
	}

	Stage1Params& Stage1Params::SetDevelopVerbose(bool developVerbose)
	{
		this->developVerbose = developVerbose;
		return *this;
	}

	Stage1Params& Stage1Params::SetSigToBinMappingPath(const std::string& sigToBinMappingPath)
	{
		this->sigToBinMappingPath = sigToBinMappingPath;
		return *this;
	}

	Stage2Params& Stage2Params::SetMaxRamGB(uint32_t maxRamGB)
	{
		if (maxRamGB < MIN_MEM)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameret: min memory must be at least " << MIN_MEM << "GB\n";
			throw std::runtime_error(err_msg.str());
		}
		this->maxRamGB = maxRamGB;
		return *this;
	}
	Stage2Params& Stage2Params::SetNThreads(uint32_t nThreads)
	{
		if (nThreads < MIN_SR || nThreads > MAX_SR)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: number of threads for 2nd stage must be in range <" << MIN_SR << "," << MAX_SR << ">";
			throw std::runtime_error(err_msg.str());					
		}
		this->nThreads = nThreads;
		return *this;
	}
	Stage2Params& Stage2Params::SetStrictMemoryMode(bool strictMemoryMode)
	{
		this->strictMemoryMode = strictMemoryMode;
		return *this;
	}
	Stage2Params& Stage2Params::SetCutoffMin(uint64_t cutoffMin)
	{
		this->cutoffMin = cutoffMin;
		return *this;
	}
	Stage2Params& Stage2Params::SetCounterMax(uint64_t counterMax)
	{
		this->counterMax = counterMax;
		return *this;
	}
	Stage2Params& Stage2Params::SetCutoffMax(uint64_t cutoffMax)
	{
		this->cutoffMax = cutoffMax;
		return *this;
	}	
	Stage2Params& Stage2Params::SetOutputFileName(const std::string& outputFileName)
	{
		this->outputFileName = outputFileName;
		return *this;
	}
	Stage2Params& Stage2Params::SetOutputFileType(OutputFileType outputFileType)
	{
		this->outputFileType = outputFileType;
		return *this;
	}
	Stage2Params& Stage2Params::SetWithoutOutput(bool withoutOutput)
	{
		this->withoutOutput = withoutOutput;
		return *this;
	}	
	Stage2Params& Stage2Params::SetStrictMemoryNSortingThreadsPerSorters(uint32_t strictMemoryNSortingThreadsPerSorters)
	{
		if (strictMemoryNSortingThreadsPerSorters < MIN_SMSO || strictMemoryNSortingThreadsPerSorters > MAX_SMSO)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: number of sorting threads per sorter in strict memory mode must be in range <" << MIN_SMSO << "," << MAX_SMSO << ">";
			throw std::runtime_error(err_msg.str());
		}

		this->strictMemoryNSortingThreadsPerSorters = strictMemoryNSortingThreadsPerSorters;
		return *this;
	}
	Stage2Params& Stage2Params::SetStrictMemoryNUncompactors(uint32_t strictMemoryNUncompactors)
	{
		if (strictMemoryNUncompactors < MIN_SMUN || strictMemoryNUncompactors > MAX_SMUN)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: number of uncompactor threads in strict memory mode must be in range <" << MIN_SMUN << "," << MAX_SMUN << ">";
			throw std::runtime_error(err_msg.str());
		}
		this->strictMemoryNUncompactors = strictMemoryNUncompactors;
		return *this;
	}
	Stage2Params& Stage2Params::SetStrictMemoryNMergers(uint32_t strictMemoryNMergers)
	{
		if (strictMemoryNMergers < MIN_SMME || strictMemoryNMergers > MAX_SMME)
		{
			std::ostringstream err_msg;
			err_msg << "Wrong parameter: number of merger threads in strict memory mode must be in range <" << MIN_SMME << "," << MAX_SMME << ">";
			throw std::runtime_error(err_msg.str());
		}
		this->strictMemoryNMergers = strictMemoryNMergers;
		return *this;
	}

	class Runner::RunnerImpl
	{
		std::unique_ptr<CApplication<KMER_WORDS>> app;
		bool stage1WasCalled = false;
	public:
		Stage1Results RunStage1(const Stage1Params& stage1Params)
		{
			stage1WasCalled = true;
#ifdef _WIN32
			_setmaxstdio(2040);
#endif			
			app = std::make_unique<CApplication<KMER_WORDS>>(stage1Params.GetKmerLen());
			return app->ProcessStage1(stage1Params);
		}

		Stage2Results RunStage2(const Stage2Params& stage2Params)
		{
			if (!stage1WasCalled)
				throw std::runtime_error("Cannot run stage 2 when stage 1 was not run");
			return app->ProcessStage2(stage2Params);
		}
	};

	Runner::Runner() : pImpl(std::make_unique<RunnerImpl>()){ };
	Runner::~Runner() = default; //must be defined here to make pImpl idiom work
	

	Stage1Results Runner::RunStage1(const Stage1Params& params)
	{		
		return pImpl->RunStage1(params);
	}

	Stage2Results Runner::RunStage2(const Stage2Params& params)	
	{		
		return pImpl->RunStage2(params);
	}

	const std::string CfgConsts::kmc_ver = KMC_VER;
	const std::string CfgConsts::kmc_date = KMC_DATE;
	const uint32_t CfgConsts::min_k = MIN_K;
	const uint32_t CfgConsts::max_k = MAX_K;
}