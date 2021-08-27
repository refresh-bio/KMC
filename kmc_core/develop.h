#ifndef _DEVELOP_H
#define _DEVELOP_H

#include "defs.h"

#define MAP_LOG_SRC "map.log"
#define KMERS_PER_BIN_LOG_FILE "kmers_per_bin.log"

void map_log(uint32 signature_len, uint32 map_size, int32* signature_map);

struct CKMCQueues;
struct CKMCParams;
void save_bins_stats(CKMCQueues& Queues, CKMCParams& Params, uint32 kmer_size, uint64 n_reads, uint32 signature_len, uint32 map_size, int32* signature_map);


#endif

// ***** EOF