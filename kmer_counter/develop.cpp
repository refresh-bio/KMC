#include "stdafx.h"

#include "develop.h"
#include "params.h"
#include <iostream>
#include <functional>
using namespace std;

void map_log(uint32 signature_len, uint32 map_size, int32* signature_map)
{
	#ifdef MAP_LOG_SRC
		FILE* mapLogFile = fopen(MAP_LOG_SRC, "w");
		char ACGT[10];
		ACGT[signature_len] = '\0';
		char symbols[] = { 'A', 'C', 'G', 'T' };
		if (!mapLogFile)
		{
			cerr << "Error: cannot save map log to file";
			exit(1);
		}
		fprintf(mapLogFile, "SIGNMATURE | ACGT | BIN NO\n");
		for (uint32 i = 0; i < map_size; ++i)
		{
	
			for (int j = signature_len - 1; j >= 0; --j)
				ACGT[signature_len - j - 1] = symbols[(i >> 2 * j) & 3];
	
			if (signature_map[i] >= 0)
				fprintf(mapLogFile, "%i\t\t%s\t%i\n", i, ACGT, signature_map[i]);
			else
				fprintf(mapLogFile, "%i\t\t%s\tDISABLED_SIGNATURE\n", i, ACGT);
		}
	
		fclose(mapLogFile);
	#endif
}



void save_bins_stats(CKMCQueues& Queues, CKMCParams& Params, uint32 kmer_size, uint64 n_reads, uint32 /*signature_len*/, uint32 map_size, int32* signature_map)
{
#ifdef KMERS_PER_BIN_LOG_FILE
	int32 bin_id;
	CMemDiskFile *file;
	string name;
	uint64 n_rec;
	uint64 n_plus_x_recs;
	uint64 n_super_kmers;
	uint64 size;

	Queues.bd->reset_reading();
	FILE* stats_file = fopen(KMERS_PER_BIN_LOG_FILE, "w");
	uint64 sum_size, sum_n_rec, sum_n_plus_x_recs, sum_n_super_kmers;
	sum_size = sum_n_rec = sum_n_plus_x_recs = sum_n_super_kmers = 0;
	if (!stats_file)
	{
		cerr << "Error: cannot open file to store kmers per bin: " << KMERS_PER_BIN_LOG_FILE << "\n";
		exit(1);
	}
	fprintf(stats_file, "%s;%s;%s;%s;%s;%s\n", "bin_id", "n_rec", "n_super_kmers", "size", "2nd stage MEM", "n_singatures");
	while ((bin_id = Queues.bd->get_next_sort_bin()) >= 0)
	{
		Queues.bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, n_super_kmers);

		// Reserve memory necessary to process the current bin at all next stages
		uint64 input_kmer_size;
		int64 kxmer_counter_size;
		uint32 kxmer_symbols;
		if (Params.max_x)
		{
			input_kmer_size = n_plus_x_recs * kmer_size;
			kxmer_counter_size = n_plus_x_recs * sizeof(uint32);
			kxmer_symbols = Params.kmer_len + Params.max_x + 1;
		}
		else
		{
			input_kmer_size = n_rec * kmer_size;
			kxmer_counter_size = 0;
			kxmer_symbols = Params.kmer_len;
		}
		uint64 max_out_recs = (n_rec + 1) / max(Params.cutoff_min, 1);

		std::function<int64(int64)> round_up_to_alignment = [](int64 x){ return (x + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT; };


		uint64 counter_size = min(BYTE_LOG(Params.cutoff_max), BYTE_LOG(Params.counter_max));		

		uint32 kmer_symbols = Params.kmer_len - Params.lut_prefix_len;
		uint64 kmer_bytes = kmer_symbols / 4;
		uint64 out_buffer_size = max_out_recs * (kmer_bytes + counter_size);

		uint32 rec_len = (kxmer_symbols + 3) / 4;

		uint64 lut_recs = 1 << (2 * Params.lut_prefix_len);
		uint64 lut_size = lut_recs * sizeof(uint64);


		size = round_up_to_alignment(size);
		input_kmer_size = round_up_to_alignment(input_kmer_size);
		out_buffer_size = round_up_to_alignment(out_buffer_size);
		kxmer_counter_size = round_up_to_alignment(kxmer_counter_size);
		lut_size = round_up_to_alignment(lut_size);


		int64 part1_size;
		int64 part2_size;

		if (rec_len % 2 == 0)
		{
			part1_size = input_kmer_size + kxmer_counter_size;
			part2_size = max(max(size, input_kmer_size), out_buffer_size + lut_size);
		}
		else
		{
			part1_size = max(input_kmer_size + kxmer_counter_size, size);
			part2_size = max(input_kmer_size, out_buffer_size + lut_size);
		}
		int64 req_size = part1_size + part2_size;

		uint64 n_signatures = 0;
		for (uint32 i = 0; i < map_size; ++i)
		{
			if (signature_map[i] == bin_id)
				++n_signatures;
		}

		fprintf(stats_file, "%i;%llu;%llu;%llu;%llu;%llu\n", bin_id, n_rec, n_super_kmers, size, (uint64)req_size, n_signatures);
		sum_size += size;
		sum_n_rec += n_rec;
		sum_n_plus_x_recs += n_plus_x_recs;
		sum_n_super_kmers += n_super_kmers;
	}

	fprintf(stats_file, "%s;%llu;%llu;%llu\n", "SUMMARY", sum_n_rec, sum_n_super_kmers, sum_size);
	fprintf(stats_file, "n_reads: %llu\n", n_reads);

	fclose(stats_file);

	Queues.bd->reset_reading();
	exit(1);
#endif
}

// ***** EOF