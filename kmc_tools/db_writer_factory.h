#ifndef _DB_WRITER_FACTORY_H
#define _DB_WRITER_FACTORY_H

#include "kmc1_db_writer.h"
#include "kff_db_writer.h"

template<unsigned SIZE>
CDbWriter<SIZE>* db_writer_factory(COutputDesc& output_desc, OutputType output_type)
{
	switch (output_type)
	{
		case OutputType::KMC1:
			return new CKMC1DbWriter<SIZE>(nullptr, output_desc);
		case OutputType::KFF1:
			return new CKFFDbWriter<SIZE>(nullptr, output_desc);
		default:
		{
			std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << "\t" << __LINE__ << "\n";
			exit(1);
		}
	}
}


#endif