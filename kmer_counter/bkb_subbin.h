/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _BKB_SUBBIN_H
#define _BKB_SUBBIN_H

//************************************************************************************************************
// CSubBin - sorted k-mers (part of some bin), used in strict memory mode 
//************************************************************************************************************
template<unsigned SIZE>
class CSubBin
{
	CDiskLogger* disk_logger;
	uchar* raw_lut;
	uint64* lut;
	uint32 current_prefix;
	uchar* suff_buff;
	uint64 suff_buff_size, max_in_suff_buff, lut_start_pos_in_file;	
	uint64 left_to_read, n_kmers, in_current_prefix;
	uint32 kmer_len, lut_size, lut_buff_recs, lut_offset, cur_in_suff_buff, suff_buff_pos;
	string name;
	FILE* file;		
	uint32 suff_rec_len, lut_prefix_len, counter_size, suffix_bytes;	
	uint64 size;
	void read_next_lut_part();
public:
	bool get_min(CKmer<SIZE>& kmer, uint32& count);
	CSubBin(CDiskLogger* _disk_logger)
	{
		lut_size = 0;
		disk_logger = _disk_logger;
	}
	void init(FILE* _file, uint64 _size, uint32 _lut_prefix_len, uint64 _n_kmers, string _name, uint32 _kmer_len, uchar* _lut_buff, uint32 _lut_buff_size, uchar* _suff_buff, uint64 _suff_buff_size);	
};

//--------------------------------------------------------------------------
template<unsigned SIZE>
void CSubBin<SIZE>::read_next_lut_part()
{
	uint32 to_read = MIN(lut_size - lut_offset, lut_buff_recs);
	lut_offset += lut_buff_recs;
	if (to_read)
	{
		uint64 prev_pos = my_ftell(file);
		my_fseek(file, lut_start_pos_in_file + (lut_offset - lut_buff_recs) * sizeof(uint64), SEEK_SET);
		if (fread(lut, sizeof(uint64), to_read, file) != to_read)
		{
			cerr << "Error while reading file : " << name << "\n";
			exit(1);
		}
		my_fseek(file, prev_pos, SEEK_SET);
	}
}

//--------------------------------------------------------------------------
template<unsigned SIZE>
bool CSubBin<SIZE>::get_min(CKmer<SIZE>& kmer, uint32& count)
{
	while (true)
	{
		if (current_prefix >= lut_offset)
		{
			read_next_lut_part();
		}
		if (in_current_prefix >= lut[current_prefix + lut_buff_recs - lut_offset])
		{
			++current_prefix;
			in_current_prefix = 0;
		}
		else
		{
			++in_current_prefix;
			break;
		}
		if (current_prefix >= lut_size)
		{
			fclose(file);
			remove(name.c_str());
			disk_logger->log_remove(size);
			return false;
		}
	}

	uchar *suf_rec = suff_buff + suff_buff_pos * suff_rec_len;
	uint32 tmp = current_prefix;
	uint32 pos = suffix_bytes;
	kmer.load(suf_rec, suffix_bytes);
	while (tmp)
	{
		kmer.set_byte(pos++, (uchar)tmp & 0xFF);
		tmp >>= 8;
	}

	count = 0;
	for (uint32 i = 0; i < counter_size; ++i)
		count += (*suf_rec++) << (8 * i);

	suff_buff_pos++;

	if (suff_buff_pos >= cur_in_suff_buff)
	{
		cur_in_suff_buff = (uint32)fread(suff_buff, 1, MIN(suff_rec_len * max_in_suff_buff, left_to_read), file) / suff_rec_len;
		suff_buff_pos = 0;
		left_to_read -= cur_in_suff_buff * suff_rec_len;
	}
	return true;
}

//--------------------------------------------------------------------------
template<unsigned SIZE>
void CSubBin<SIZE>::init(FILE* _file, uint64 _size, uint32 _lut_prefix_len, uint64 _n_kmers, string _name, uint32 _kmer_len, uchar* _lut_buff, uint32 _lut_buff_size, uchar* _suff_buff, uint64 _suff_buff_size)
{
	size = _size;
	lut = (uint64*)_lut_buff;
	lut_buff_recs = _lut_buff_size / sizeof(uint64);
	suff_buff = _suff_buff;
	suff_buff_size = _suff_buff_size;
	lut_offset = 0;

	lut_prefix_len = _lut_prefix_len;
	kmer_len = _kmer_len;
	suffix_bytes = (kmer_len - lut_prefix_len) / 4;
	file = _file;
	n_kmers = _n_kmers;
	name = _name;
	counter_size = sizeof(uint32);

	lut_size = (1 << lut_prefix_len * 2);

	suff_rec_len = (kmer_len - lut_prefix_len) / 4 + counter_size;
	left_to_read = suff_rec_len * n_kmers;
	max_in_suff_buff = suff_buff_size / suff_rec_len;
	lut_start_pos_in_file = n_kmers * suff_rec_len;
	rewind(file);
	read_next_lut_part();
	cur_in_suff_buff = (uint32)fread(suff_buff, 1, MIN(max_in_suff_buff * suff_rec_len, left_to_read), file) / suff_rec_len;
	left_to_read -= cur_in_suff_buff * suff_rec_len;
	current_prefix = 0;
	in_current_prefix = 0;
	suff_buff_pos = 0;
}

#endif

// ***** EOF