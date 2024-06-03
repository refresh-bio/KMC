#ifndef _DB_READER_FACTORY_H
#define _DB_READER_FACTORY_H

#include "bundle.h"
#include "kmer_file_header.h"
#include "kmc1_db_reader.h"
#include "kmc2_db_reader.h"
#include "kff_db_reader.h"

template<unsigned SIZE>
CInput<SIZE>* db_reader_factory(const CKmerFileHeader& header, const CInputDesc& input_desc, KmerDBOpenMode kmer_db_open_mode)
{
	switch (header.kmer_file_type)
	{
		case KmerFileType::KMC1:
			return new CKMC1DbReader<SIZE>(header, input_desc, CConfig::GetInstance().percent_progress, kmer_db_open_mode);
		case KmerFileType::KMC2:
			return new CKMC2DbReader<SIZE>(header, input_desc, CConfig::GetInstance().percent_progress, kmer_db_open_mode);
		case KmerFileType::KFF1:
			return new CKFFDbReader<SIZE>(header, input_desc, CConfig::GetInstance().percent_progress, kmer_db_open_mode);
		default:
		{
			std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << "\t" << __LINE__ << "\n";
			exit(1);
		}
	}
}

#endif // !_DB_READER_FACTORY_H



