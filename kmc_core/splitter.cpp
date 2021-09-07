/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.1
Date   : 2019-05-19
*/

#include "stdafx.h"
#include "splitter.h"

//************************************************************************************************************
// CSplitter class - splits kmers into bins according to their signatures
//************************************************************************************************************

//uint32 CSplitter::MAX_LINE_SIZE = 1 << 14;
uint32 CSplitter::MAX_LINE_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
// Assigns queues
CSplitter::CSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	//mm = Queues.mm;
	file_type = Params.file_type;
	both_strands = Params.both_strands;

	bin_part_queue = Queues.bpq.get();
	pmm_reads = Queues.pmm_reads.get();
	kmer_len = Params.kmer_len;
	signature_len = Params.signature_len;

	mem_part_pmm_bins = Params.mem_part_pmm_bins;

	mem_part_pmm_reads = Params.mem_part_pmm_reads; 

	s_mapper = Queues.s_mapper.get();

	part = nullptr;

	// Prepare encoding of symbols
	for (int i = 0; i < 256; ++i)
		codes[i] = -1;
	codes['A'] = codes['a'] = 0;
	codes['C'] = codes['c'] = 1;
	codes['G'] = codes['g'] = 2;
	codes['T'] = codes['t'] = 3;

	n_reads = 0;

	homopolymer_compressed = Params.homopolymer_compressed;

	ntHashEstimator = Queues.ntHashEstimator.get();
}

void CSplitter::InitBins(CKMCParams &Params, CKMCQueues &Queues)
{
	n_bins = Params.n_bins;
	uint32 buffer_size = Params.bin_part_size;
	// Create objects for all bin
	bins.resize(n_bins);
	for (uint32 i = 0; i < n_bins; ++i)
	{
		bins[i] = std::make_unique<CKmerBinCollector>(Queues, Params, buffer_size, i);
	}
}

//----------------------------------------------------------------------------------
// Parse long read, header_merker is '@' or '>'
bool CSplitter::GetSeqLongRead(char *seq, uint32 &seq_size, uchar header_marker)
{	
	uint32 pos = 0;
	//long read may or may not contain header
	if (part_pos == 0 && part[0] == header_marker)
	{
		++n_reads;
		for (; part[part_pos] != '\n' && part[part_pos] != '\r'; ++part_pos)
			;
	}
	while (pos < mem_part_pmm_reads && part_pos < part_size)
		seq[pos++] = codes[part[part_pos++]];
	seq_size = pos;
	if (part_pos < part_size)
		part_pos -= kmer_len - 1;
	return true;
}



//----------------------------------------------------------------------------------
// Return a single record from FASTA/FASTQ data
bool CSplitter::GetSeq(char *seq, uint32 &seq_size, ReadType read_type)
{
	if (part_pos >= part_size)
		return false;

	uchar c = 0;
	uint32 pos = 0;

	if (file_type == InputType::FASTA || file_type == InputType::KMC)
	{		
		if (read_type == ReadType::long_read)
			return GetSeqLongRead(seq, seq_size, '>');
		if (curr_read_len == 0)
		{
			// Title
			c = part[part_pos++];
			if (c != '>')
				return false;
			++n_reads;

			for (; part_pos < part_size;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
			}
			if (part_pos >= part_size)
				return false;

			c = part[part_pos++];
			if (c >= 32 || c == part[part_pos - 2]) //read may be empty
				part_pos--;
			else if (part_pos >= part_size)
				return false;

			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}

			seq_size = pos;

			if (part_pos >= part_size)
				return true;

			curr_read_len = pos;

			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}
		else // we are inside read
		{
			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}

			seq_size = pos;

			if (part_pos >= part_size)
				return true;

			curr_read_len += pos - kmer_len + 1;

			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}
		
		curr_read_len = 0;

		//end of last record 
		if (part_pos >= part_size)
			return true;

		if (part[part_pos++] >= 32)
			part_pos--;
		else if (part_pos >= part_size)
			return true;
	}
	else if (file_type == InputType::FASTQ)
	{
		if (read_type == ReadType::long_read)		
			return GetSeqLongRead(seq, seq_size, '@');	
	
		if (curr_read_len == 0)
		{
			// Title
			c = part[part_pos++];
			if (c != '@')
				return false;
			++n_reads;

			for (; part_pos < part_size;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
			}
			if (part_pos >= part_size)
				return false;

			c = part[part_pos++];
			if (c >= 32 || c == part[part_pos - 2]) //read may be empty
				part_pos--;
			else if (part_pos >= part_size)
				return false;

			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}
			if (part_pos >= part_size)
				return false;

			seq_size = pos;
			curr_read_len = pos;

			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}
		else // we are inside read
		{
			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}
			if (part_pos >= part_size)
				return false;

			seq_size = pos;
			curr_read_len += pos - kmer_len + 1;
			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}

		c = part[part_pos++];
		if (c >= 32)
			part_pos--;
		else if (part_pos >= part_size)
			return false;

		// Plus
		c = part[part_pos++];
		if (part_pos >= part_size)
			return false;
		if (c != '+')
			return false;
		for (; part_pos < part_size;)
		{
			c = part[part_pos++];
			if (c < 32)					// newliners
				break;
		}
		if (part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if (c >= 32 || c == part[part_pos - 2]) //qual may be empty
			part_pos--;
		else if (part_pos >= part_size)
			return false;

		// Quality
		part_pos += curr_read_len; //skip quality
		
		curr_read_len = 0;

		if (part_pos >= part_size)
			return false;
		c = part[part_pos++];

		//end of last record 
		if (part_pos >= part_size)
			return true;

		//may be additional EOL character 
		if (part[part_pos++] >= 32)
			part_pos--;
		else if (part_pos >= part_size)
			return true;
		
	}
	else if (file_type == InputType::MULTILINE_FASTA)
	{
		if (part[part_pos] == '>')//need to ommit header
		{
			++n_reads;
			for (; part_pos < part_size && part[part_pos] != '\n' && part[part_pos] != '\r'; ++part_pos);//find EOF
			++part_pos;
			if (part[part_pos] == '\n' || part[part_pos] == '\r')
				++part_pos;
		}
		for (; part_pos < part_size && pos < mem_part_pmm_reads && part[part_pos] != '>';)
		{
			seq[pos++] = codes[part[part_pos++]];
		}
		seq_size = pos;
		if (part_pos < part_size && part[part_pos] != '>')//need to copy last k-1 kmers 
		{
			part_pos -= kmer_len - 1;
		}
		return true;

	}
	else if (file_type == InputType::BAM)
	{
		while (true)
		{
			if (part_pos >= part_size)
				return false;

			int32_t block_size;
			read_int32_t(block_size, part, part_pos);

			uint64_t start_pos = part_pos;

			part_pos += 8;

			uint32_t bin_mq_nl;
			read_uint32_t(bin_mq_nl, part, part_pos);

			uint32_t l_read_name = (bin_mq_nl & ((1 << 8) - 1));
			uint32_t flag_nc;
			read_uint32_t(flag_nc, part, part_pos);
			uint32_t n_cigar_op = flag_nc & ((1ul << 16) - 1);
			int32_t l_seq;
			read_int32_t(l_seq, part, part_pos);

			part_pos += 12;

			uint32_t flags = flag_nc >> 16;

			bool exclude_read = ((flags >> 8) & 1) || ((flags >> 11) & 1); //TODO: I think that is the way samtools filter out some reads (secondary and supplementary)

			part_pos += l_read_name; // skip read name

			part_pos += 4 * n_cigar_op;
			if (!exclude_read)
			{
				bool is_rev_comp = (flags >> 4) & 1;

				if (!both_strands && is_rev_comp) //if read is reversed and kmc was run to count all (not only canonical) kmers read must be transformed back
				{
					//static const char rev_maping[] = "=TGMCRSVAWYHKDBN";
					static const char rev_maping[] = { -1, 3, 2, -1, 1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1 };// "=TGMCRSVAWYHKDBN";
					uint32 n_bytes = l_seq / 2;
					uint64_t pos_after = pos + l_seq;
					pos = pos_after;
					for (uint32_t ii = 0; ii < n_bytes; ++ii)
					{
						unsigned char byte = part[part_pos++];
						seq[--pos_after] = rev_maping[byte >> 4];
						seq[--pos_after] = rev_maping[byte & 15];
					}

					if (l_seq & 1) //odd
					{
						unsigned char byte = part[part_pos++];
						seq[--pos_after] = rev_maping[byte >> 4];
					}
				}
				else
				{
					static const char maping[] = { -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1 };//"=ACMGRSVTWYHKDBN";
					uint32 n_bytes = l_seq / 2;
					for (uint32_t ii = 0; ii < n_bytes; ++ii)
					{
						unsigned char byte = part[part_pos++];
						seq[pos++] = maping[byte >> 4];
						seq[pos++] = maping[byte & 15];
					}

					if (l_seq & 1) //odd
					{
						unsigned char byte = part[part_pos++];
						seq[pos++] = maping[byte >> 4];
					}
				}
				seq_size = pos;
			}
			else
			{
				part_pos += (l_seq + 1) / 2;
			}

			//move to next record
			uint64_t readed = part_pos - start_pos;
			uint64_t remaining = block_size - readed;
			part_pos += remaining;

			if (!exclude_read) //if readed successfuly return		
			{
				++n_reads;
				return true;
			}
		}

	}
	return (c == '\n' || c == '\r');
}

//----------------------------------------------------------------------------------
void CSplitter::HomopolymerCompressSeq(char* seq, uint32 &seq_size)
{
	if (seq_size <= 1)
		return;

	uint32 read_pos = 1;
	uint32 write_pos = 0;
	for (; read_pos < seq_size; ++read_pos)
		if (seq[read_pos] != seq[write_pos])
			seq[++write_pos] = seq[read_pos];
	seq_size = write_pos + 1;
}

//----------------------------------------------------------------------------------
// Calculate statistics of m-mers
void CSplitter::CalcStats(uchar* _part, uint64 _part_size, ReadType read_type, uint32* _stats)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char *seq;
	uint32 seq_size;
	pmm_reads->reserve(seq);

	uint32 signature_start_pos;
	CMmer current_signature(signature_len), end_mmer(signature_len);

	uint32 i;
	uint32 len;//length of extended kmer

	while (GetSeq(seq, seq_size, read_type))
	{
		if (homopolymer_compressed)
			HomopolymerCompressSeq(seq, seq_size);
		i = 0;
		len = 0;
		while (i + kmer_len - 1 < seq_size)
		{
			bool contains_N = false;
			//building first signature after 'N' or at the read begining
			for (uint32 j = 0; j < signature_len; ++j, ++i)
				if (seq[i] < 0)//'N'
				{
					contains_N = true;
					break;
				}
			//signature must be shorter than k-mer so if signature contains 'N', k-mer will contains it also
			if (contains_N)
			{
				++i;
				continue;
			}
			len = signature_len;
			signature_start_pos = i - signature_len;
			current_signature.insert(seq + signature_start_pos);
			end_mmer.set(current_signature);
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)//'N'
				{
					if (len >= kmer_len)
						_stats[current_signature.get()] += 1 + len - kmer_len;
					len = 0;
					++i;
					break;
				}
				end_mmer.insert(seq[i]);
				if (end_mmer < current_signature)//signature at the end of current k-mer is lower than current
				{
					if (len >= kmer_len)
					{
						_stats[current_signature.get()] += 1 + len - kmer_len;
						len = kmer_len - 1;
					}
					current_signature.set(end_mmer);
					signature_start_pos = i - signature_len + 1;
				}
				else if (end_mmer == current_signature)
				{
					current_signature.set(end_mmer);
					signature_start_pos = i - signature_len + 1;
				}
				else if (signature_start_pos + kmer_len - 1 < i)//need to find new signature
				{
					_stats[current_signature.get()] += 1 + len - kmer_len;
					len = kmer_len - 1;
					//looking for new signature
					++signature_start_pos;
					//building first signature in current k-mer
					end_mmer.insert(seq + signature_start_pos);
					current_signature.set(end_mmer);
					for (uint32 j = signature_start_pos + signature_len; j <= i; ++j)
					{
						end_mmer.insert(seq[j]);
						if (end_mmer <= current_signature)
						{
							current_signature.set(end_mmer);
							signature_start_pos = j - signature_len + 1;
						}
					}
				}
				++len;
			}
		}
		if (len >= kmer_len)//last one in read
			_stats[current_signature.get()] += 1 + len - kmer_len;
	}
	pmm_reads->free(seq);
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part, but only for k-mer occurence estimation
bool CSplitter::ProcessReadsOnlyEstimate(uchar* _part, uint64 _part_size, ReadType read_type)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char* seq;
	uint32 seq_size;
	pmm_reads->reserve(seq);

	while (GetSeq(seq, seq_size, read_type))
		ntHashEstimator->Process(seq, seq_size);

	pmm_reads->free(seq);

	return true;
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part
bool CSplitter::ProcessReads(uchar *_part, uint64 _part_size, ReadType read_type)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char *seq;
	uint32 seq_size;
	pmm_reads->reserve(seq);

	uint32 signature_start_pos;
	CMmer current_signature(signature_len), end_mmer(signature_len);
	uint32 bin_no;

	uint32 i;
	uint32 len;//length of extended kmer

	while (GetSeq(seq, seq_size, read_type))
	{		
		if (ntHashEstimator)
			ntHashEstimator->Process(seq, seq_size);

		if (homopolymer_compressed)
			HomopolymerCompressSeq(seq, seq_size);
		//if (file_type != multiline_fasta && file_type != fastq) //read conting moved to GetSeq
		//	n_reads++;
		i = 0;
		len = 0;
		while (i + kmer_len - 1 < seq_size)
		{
			bool contains_N = false;
			//building first signature after 'N' or at the read begining
			for (uint32 j = 0; j < signature_len; ++j, ++i)
				if (seq[i] < 0)//'N'
				{
					contains_N = true;
					break;
				}
			//signature must be shorter than k-mer so if signature contains 'N', k-mer will contains it also
			if (contains_N)
			{
				++i;
				continue;
			}
			len = signature_len;
			signature_start_pos = i - signature_len;
			current_signature.insert(seq + signature_start_pos);
			end_mmer.set(current_signature);
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)//'N'
				{
					if (len >= kmer_len)
					{
						bin_no = s_mapper->get_bin_id(current_signature.get());
						bins[bin_no]->PutExtendedKmer(seq + i - len, len);
					}
					len = 0;
					++i;
					break;
				}
				end_mmer.insert(seq[i]);
				if (end_mmer < current_signature)//signature at the end of current k-mer is lower than current
				{
					if (len >= kmer_len)
					{
						bin_no = s_mapper->get_bin_id(current_signature.get());
						bins[bin_no]->PutExtendedKmer(seq + i - len, len);
						len = kmer_len - 1;
					}
					current_signature.set(end_mmer);
					signature_start_pos = i - signature_len + 1;
				}
				else if (end_mmer == current_signature)
				{
					current_signature.set(end_mmer);
					signature_start_pos = i - signature_len + 1;
				}
				else if (signature_start_pos + kmer_len - 1 < i)//need to find new signature
				{
					bin_no = s_mapper->get_bin_id(current_signature.get());
					bins[bin_no]->PutExtendedKmer(seq + i - len, len);
					len = kmer_len - 1;
					//looking for new signature
					++signature_start_pos;
					//building first signature in current k-mer
					end_mmer.insert(seq + signature_start_pos);
					current_signature.set(end_mmer);
					for (uint32 j = signature_start_pos + signature_len; j <= i; ++j)
					{
						end_mmer.insert(seq[j]);
						if (end_mmer <= current_signature)
						{
							current_signature.set(end_mmer);
							signature_start_pos = j - signature_len + 1;
						}
					}
				}
				++len;
				if (len == kmer_len + 255) //one byte is used to store counter of additional symbols in extended k-mer
				{
					bin_no = s_mapper->get_bin_id(current_signature.get());
					bins[bin_no]->PutExtendedKmer(seq + i + 1 - len, len);
					i -= kmer_len - 2;
					len = 0;
					break;
				}

			}
		}
		if (len >= kmer_len)//last one in read
		{
			bin_no = s_mapper->get_bin_id(current_signature.get());
			bins[bin_no]->PutExtendedKmer(seq + i - len, len);
		}
	}

	pmm_reads->free(seq);

	return true;
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part in small k optimization mode
template<typename COUNTER_TYPE>
bool CSplitter::ProcessReadsSmallK(uchar *_part, uint64 _part_size, ReadType read_type, CSmallKBuf<COUNTER_TYPE>& small_k_buf)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char *seq;
	uint32 seq_size;
	int omit_next_n_kmers;
	CKmer<1> kmer_str, kmer_rev, kmer_can;
	uint32 i;
	CKmer<1> kmer_mask;
	pmm_reads->reserve(seq);
	kmer_mask.set_n_1(2 * kmer_len);

	uint32 kmer_len_shift = (kmer_len - 1) * 2;

	if (both_strands)
		while (GetSeq(seq, seq_size, read_type))
		{
			if (homopolymer_compressed)
				HomopolymerCompressSeq(seq, seq_size);
			//if (file_type != multiline_fasta)
			//	n_reads++;

			// Init k-mer
			kmer_str.clear();
			kmer_rev.clear();

			// Process first k-1 symbols of a read
			uint32 str_pos = kmer_len_shift - 2;
			uint32 rev_pos = 2;

			omit_next_n_kmers = 0;

			for (i = 0; i < kmer_len - 1; ++i, str_pos -= 2, rev_pos += 2)
			{
				if (seq[i] < 0)
				{
					seq[i] = 0;
					omit_next_n_kmers = i + 1;
				}
				kmer_str.set_2bits(seq[i], str_pos);
				kmer_rev.set_2bits(3 - seq[i], rev_pos);
			}

			// Process next part of a read
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)		// N in a read
				{
					seq[i] = 0;
					omit_next_n_kmers = kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
				}
				kmer_str.SHL_insert_2bits(seq[i]);
				kmer_str.mask(kmer_mask);
				kmer_rev.SHR_insert_2bits(3 - seq[i], kmer_len_shift);

				// If necessary ommit next symbols
				if (omit_next_n_kmers > 0)
				{
					omit_next_n_kmers--;
					continue;
				}

				// Find canonical kmer representation
				kmer_can = (kmer_str < kmer_rev) ? kmer_str : kmer_rev;

				++small_k_buf.buf[kmer_can.data];
				++total_kmers;
			}
		}
	else
		while (GetSeq(seq, seq_size, read_type))
		{
			if (homopolymer_compressed)
				HomopolymerCompressSeq(seq, seq_size);
			//if (file_type != multiline_fasta)
			//	n_reads++;

			// Init k-mer
			kmer_str.clear();

			// Process first k-1 symbols of a read
			uint32 str_pos = kmer_len_shift - 2;

			omit_next_n_kmers = 0;

			for (i = 0; i < kmer_len - 1; ++i, str_pos -= 2)
			{
				if (seq[i] < 0)
				{
					seq[i] = 0;
					omit_next_n_kmers = i + 1;
				}
				kmer_str.set_2bits(seq[i], str_pos);
			}

			// Process next part of a read
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)		// N in a read
				{
					seq[i] = 0;
					omit_next_n_kmers = kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
				}
				kmer_str.SHL_insert_2bits(seq[i]);
				kmer_str.mask(kmer_mask);

				// If necessary ommit next symbols
				if (omit_next_n_kmers > 0)
				{
					omit_next_n_kmers--;
					continue;
				}

				++small_k_buf.buf[kmer_str.data];
				++total_kmers;
			}
		}

	pmm_reads->free(seq);
	return true;
}

//----------------------------------------------------------------------------------
// Finish the processing of input file
void CSplitter::Complete()
{
	for (auto& bin : bins)
		if (bin)
			bin->Flush();
}

//************************************************************************************************************
// CWSplitter class - wrapper for multithreading purposes
//************************************************************************************************************
//----------------------------------------------------------------------------------
// Constructor
CWSplitter::CWSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	pq = Queues.part_queue.get();
	bpq = Queues.bpq.get();
	pmm_fastq = Queues.pmm_fastq.get();
	spl = std::make_unique<CSplitter>(Params, Queues);
	spl->InitBins(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
void CWSplitter::operator()()
{
	// Splitting parts
	while (!pq->completed())
	{
		uchar *part;
		uint64 size;
		ReadType read_type;
		if (pq->pop(part, size, read_type))
		{			
			spl->ProcessReads(part, size, read_type);
			pmm_fastq->free(part);
		}
	}
	spl->Complete();
	bpq->mark_completed();

	spl->GetTotal(n_reads);

	spl.reset();
}

//----------------------------------------------------------------------------------
// Destructor
CWSplitter::~CWSplitter()
{
}

//----------------------------------------------------------------------------------
// Return statistics
void CWSplitter::GetTotal(uint64 &_n_reads)
{
	if (spl)
		spl->GetTotal(n_reads);

	_n_reads = n_reads;
}


//************************************************************************************************************
// CWStatsSplitter class - wrapper for multithreading purposes
//************************************************************************************************************


//----------------------------------------------------------------------------------
// Constructor
CWStatsSplitter::CWStatsSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	spq = Queues.stats_part_queue.get();
	pmm_fastq = Queues.pmm_fastq.get();
	pmm_stats = Queues.pmm_stats.get();
	progressObserver = Params.progressObserver;
	spl = std::make_unique<CSplitter>(Params, Queues);

	signature_len = Params.signature_len;
	pmm_stats->reserve(stats);
	fill_n(stats, (1 << signature_len * 2) + 1, 0);
}

//----------------------------------------------------------------------------------
// Destructor
CWStatsSplitter::~CWStatsSplitter()
{
	pmm_stats->free(stats);
}

//----------------------------------------------------------------------------------
// Execution
void CWStatsSplitter::operator()()
{
	// Splitting parts
	while (!spq->completed())
	{
		uchar *part;
		uint64 size;
		ReadType read_type;
		if (spq->pop(part, size, read_type))
		{
			spl->CalcStats(part, size, read_type, stats);
			progressObserver->Step();
			pmm_fastq->free(part);
		}
	}

	spl.reset();
}

//----------------------------------------------------------------------------------
void CWStatsSplitter::GetStats(uint32* _stats)
{
	uint32 size = (1 << signature_len * 2) + 1;
	for (uint32 i = 0; i < size; ++i)
		_stats[i] += stats[i];
}


//************************************************************************************************************
// CWSmallKSplitter class - wrapper for multithreading purposes
//************************************************************************************************************


//----------------------------------------------------------------------------------
// Constructor
template <typename COUNTER_TYPE> CWSmallKSplitter<COUNTER_TYPE>::CWSmallKSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	pq = Queues.part_queue.get();
	pmm_fastq = Queues.pmm_fastq.get();
	pmm_small_k = Queues.pmm_small_k_buf.get();
	kmer_len = Params.kmer_len;
	spl = std::make_unique<CSplitter>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
template <typename COUNTER_TYPE> CWSmallKSplitter<COUNTER_TYPE>::~CWSmallKSplitter()
{
}

//----------------------------------------------------------------------------------
// Execution
template <typename COUNTER_TYPE> void CWSmallKSplitter<COUNTER_TYPE>::operator()()
{
	pmm_small_k->reserve(small_k_buf.buf);
	memset(small_k_buf.buf, 0, (1ull << 2 * kmer_len) * sizeof(*small_k_buf.buf));

	// Splitting parts
	while (!pq->completed())
	{
		uchar *part;
		uint64 size;
		ReadType read_type;
		if (pq->pop(part, size, read_type))
		{
			spl->ProcessReadsSmallK(part, size, read_type, small_k_buf);
			pmm_fastq->free(part);
		}
	}
	spl->Complete();

	spl->GetTotal(n_reads);
	total_kmers = spl->GetTotalKmers();
	spl.reset();
}

//----------------------------------------------------------------------------------
// Return statistics
template <typename COUNTER_TYPE> void CWSmallKSplitter<COUNTER_TYPE>::GetTotal(uint64 &_n_reads)
{
	if (spl)
		spl->GetTotal(n_reads);
	_n_reads = n_reads;
}




//************************************************************************************************************
// CWSplitter class - wrapper for multithreading purposes
//************************************************************************************************************
//----------------------------------------------------------------------------------
// Constructor
CWEstimateOnlySplitter::CWEstimateOnlySplitter(CKMCParams& Params, CKMCQueues& Queues)
{
	pq = Queues.part_queue.get();
	pmm_fastq = Queues.pmm_fastq.get();
	spl = std::make_unique<CSplitter>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
void CWEstimateOnlySplitter::operator()()
{
	// Splitting parts
	while (!pq->completed())
	{
		uchar* part;
		uint64 size;
		ReadType read_type;
		if (pq->pop(part, size, read_type))
		{
			spl->ProcessReadsOnlyEstimate(part, size, read_type);
			pmm_fastq->free(part);
		}
	}

	spl->GetTotal(n_reads);

	spl.reset();
}

//----------------------------------------------------------------------------------
// Destructor
CWEstimateOnlySplitter::~CWEstimateOnlySplitter()
{
}

//----------------------------------------------------------------------------------
// Return statistics
void CWEstimateOnlySplitter::GetTotal(uint64& _n_reads)
{
	if (spl)
		spl->GetTotal(n_reads);

	_n_reads = n_reads;
}





//instantiate some templates
template bool CSplitter::ProcessReadsSmallK(uchar *_part, uint64 _part_size, ReadType read_type, CSmallKBuf<uint32>& small_k_buf);
template bool CSplitter::ProcessReadsSmallK(uchar *_part, uint64 _part_size, ReadType read_type, CSmallKBuf<uint64>& small_k_buf);
template class CWSmallKSplitter<uint32>;
template class CWSmallKSplitter<uint64>;

// ***** EOF