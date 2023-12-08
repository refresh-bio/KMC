/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Marek Kokot

Version: 3.2.3
Date   : 2023-12-08
*/

//TODO KFF: consider separate thread for writing linke in kmc1 db writer

#ifndef _KFF_DB_WRITER_H
#define _KFF_DB_WRITER_H

#include "../kmc_core/kff_writer.h"
#include "defs.h"


template<unsigned SIZE> class CKFFDbWriter : public CDbWriter<SIZE>
{
	uint64 counter_size;
	uint64 kmer_bytes;
	uint64 rec_bytes;
	CBundle<SIZE>* bundle;
	CKFFWriter kff_writer;
	COutputDesc& output_desc;
	uchar* buff;	
	uint64 buff_size;
	uint64 buff_pos;
	
	void store_buff()
	{
		kff_writer.StoreSectionPart(buff, buff_pos);
		buff_pos = 0;
	}

	void ProcessBundleData(CBundleData<SIZE>& bundle_data)
	{
		while (!bundle_data.Empty())
		{
			add_kmer(bundle_data.TopKmer(), bundle_data.TopCounter());
			bundle_data.Pop();
		}
		bundle_data.Clear();
	}

	void start()
	{
		kff_writer.InitSection();
	}

	void finish()
	{
		if (buff_pos != 0)
			store_buff();
		kff_writer.FinishSection();		
	}

	bool canonical()
	{
		auto& headers = CConfig::GetInstance().headers;
		bool both_stands = true;

		for (auto& input : headers)
			both_stands = both_stands && input.both_strands; //if any input database is in both strands, output is also in both strands

		return both_stands;
	}

public:
	
	CKFFDbWriter(CBundle<SIZE>* bundle, COutputDesc& output_desc) :
		counter_size(MIN(BYTE_LOG(output_desc.counter_max), BYTE_LOG(output_desc.cutoff_max))),
		kmer_bytes((CConfig::GetInstance().kmer_len + 3) / 4),
		rec_bytes(kmer_bytes + counter_size),
		bundle(bundle),
		kff_writer(
			output_desc.file_src + ".kff",
			canonical(),
			CConfig::GetInstance().kmer_len,
			counter_size,
			output_desc.cutoff_min,
			output_desc.cutoff_max,
			output_desc.encoding),
		output_desc(output_desc)
	{
		buff = new uchar[KFF_DB_WRITER_BUFF_BYTES];
		rec_bytes = kmer_bytes + counter_size;
		buff_size = KFF_DB_WRITER_BUFF_BYTES / rec_bytes;
		buff_pos = 0;
	}
	~CKFFDbWriter() 
	{
		delete[] buff;
	}

	void add_kmer(CKmer<SIZE>& kmer, uint32 counter)
	{
		//if specific counter value is set use it as counter value (set_counts operation), do not check if counter is valid in term of cutoffs and counter max
		if (output_desc.counter_value)
			counter = static_cast<uint32>(output_desc.counter_value);
		else
		{
			if (counter < output_desc.cutoff_min || counter > output_desc.cutoff_max)
				return;
			if (counter > output_desc.counter_max)
				counter = output_desc.counter_max;
		}
		uchar* rec = buff + buff_pos * rec_bytes;
		kmer.store(rec, kmer_bytes);

		for (int32 j = (int32)counter_size - 1; j >= 0; --j)
			*rec++ = (counter >> (j * 8)) & 0xFF;
		++buff_pos;
		if (buff_size == buff_pos)
			store_buff();		
	}

	bool Process() override 
	{
		start();
		while (!bundle->Finished())
			ProcessBundleData(bundle->Data());			
		
		finish();
		return true;
	}

	void MultiOptputInit() override 
	{
		start();
	}

	void MultiOptputAddResultPart(COutputBundle<SIZE>& bundle) override 
	{
		ProcessBundleData(bundle.Data());
	}
	
	void MultiOptputAddResultPart(CBundle<SIZE>& bundle) override 
	{
		ProcessBundleData(bundle.Data());
	}
	
	void MultiOptputFinish() override 
	{
		finish();
	}
};

#endif