/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.4
  Date   : 2024-02-09
 */

#ifndef _KMC_FILE_H
#define _KMC_FILE_H

#include "kmer_defs.h"
#include "kmer_api.h"
#include "../kmc_api/sig_to_bin_map.h"
#include "../kmc_core/kmc_runner.h"
#include "../kmc_core/kmer.h"
#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <sstream>
#include <string>
#include <limits>

struct CKMCFileInfo
{
	uint32 kmer_length;
	uint32 mode;
	uint32 counter_size;
	uint32 lut_prefix_length;
	uint32 signature_len;
	KMC::SignatureSelectionScheme signature_selection_scheme;
	uint32_t n_bins;
	uint32 min_count;
	uint64 max_count;
	bool both_strands;
	uint64 total_kmers;
};

class CKMCFile
{
	//a FILE* opened for read binary wrapper that may be configured to reopen for
	//each operation and then close
	class FILEWrapper
	{
		bool reopen_each_time;
		std::string path;
		using read_pos_t = decltype(my_ftell(nullptr));
		size_t file_size{};
		read_pos_t read_pos{}; // if we reopen file at each operation, for reading we need to know where last read ended
		FILE* file{};
		bool is_opened = false; //not necesarly opened, because if we are reopening, its not opened per se, but it means the Open was called and the Close was not

		void close_internal()
		{
			if(file)
			{
				fclose(file);
				file = nullptr;
			}
		}
		bool Open()
		{
			assert(!file);
			file = fopen(path.c_str(), "rb");
			if (!file)
				return false;
			setvbuf(file, nullptr, _IONBF, 0);
//			std::cerr << "File " << path << " opened and move to pos " << read_pos << "\n";
			my_fseek(file, read_pos, SEEK_SET);
			return true;
		}
	public:
		bool Open(const std::string& path)
		{
			this->path = path;
			read_pos = 0;

			if (!Open())
				return false;

			my_fseek(file, 0, SEEK_END);
			file_size = my_ftell(file);
			my_fseek(file, 0, SEEK_SET);

			if (reopen_each_time)
				close_internal();

			is_opened = true;

			return true;
		}

		//offset from start
		void Seek(long long offset)
		{
			if (reopen_each_time)
			{
				assert(!file);
				read_pos = offset;
			}
			else
			{
				assert(file);
				my_fseek(file, offset, SEEK_SET);
			}
		}

		size_t Size() const
		{
			return file_size;
		}

		int Getc()
		{
			if (reopen_each_time)
			{
				Open();
				auto res = fgetc(file);
				close_internal();
				return res;
			}
			else
				return fgetc(file);
		}

		void Rewind()
		{
			if (reopen_each_time)
			{
				read_pos = 0;
			}
			else
			{
				assert(file);
				rewind(file);
			}
		}

		bool IsOpened() const
		{
			return is_opened;
		}

		size_t Read(void* ptr, size_t size, size_t count)
		{
			if (reopen_each_time)
			{
				//it should not fail, but who know, maybe someone deleted the file? ro changed permissions
				if (!Open())
				{
					std::cerr << "Error: Cannot open file " << path << "\n";
					exit(1);
				}
				auto res = fread(ptr, size, count, file);
				read_pos = my_ftell(file);
				close_internal();
				return res;
			}
			else
				return fread(ptr, size, count, file);
		}

		void Close()
		{
			close_internal();
			path = "";
			read_pos = 0;
			is_opened = false;
		}

		FILEWrapper(bool reopen_each_time):reopen_each_time(reopen_each_time)
		{

		}

		FILEWrapper(const FILEWrapper& rhs) = delete;
		FILEWrapper& operator=(const FILEWrapper& rhs) = delete;

		FILEWrapper(FILEWrapper&& rhs) noexcept:
			reopen_each_time(rhs.reopen_each_time),
			path(std::move(rhs.path)),
			read_pos(rhs.read_pos),
			file(rhs.file)
		{
			rhs.file = nullptr;
		}

		~FILEWrapper()
		{
			close_internal();
		}
	};

	class CPrefixFileBufferForListingMode
	{
		const uint64_t buffCapacity = 1 << 22;
		uint64_t* buff{};
		uint64_t buffPosInFile{};
		uint64_t buffSize{};
		uint64_t posInBuf{};
		uint64_t leftToRead{};
		uint64 prefixMask; //for kmc2 db
		FILEWrapper& file;
		bool isKMC1 = false;
		uint64_t totalKmers; //for

		void reload()
		{
			assert(leftToRead);
			buffPosInFile += buffSize;
			buffSize = (std::min)(buffCapacity, leftToRead);
			auto readed = file.Read(buff, 1, 8 * buffSize);
			assert(readed == 8 * buffSize);

			if (isKMC1 && buffSize == leftToRead) //last read, in case of KMC1 guard must be added, fread will read `k` from db instead of guard, fixes #180
				buff[buffSize - 1] = totalKmers;

			leftToRead -= buffSize;
			posInBuf = 0;
		}
	public:
		CPrefixFileBufferForListingMode(FILEWrapper& file, uint64_t wholeLutSize, uint64_t lutPrefixLen, bool isKMC1, uint64_t totalKmers)
			:
			buff(new uint64_t[buffCapacity]),
			leftToRead(wholeLutSize),
			prefixMask((1ull << (2 * lutPrefixLen)) - 1),
			file(file),
			isKMC1(isKMC1),
			totalKmers(totalKmers)
		{
			file.Seek(4 + 8); //	skip KMCP and LUT[0] (always = 0)
		}

		//no control if next prefix exists here, responsibility to the caller
		uint64_t GetPrefix(uint64_t suffix_number)
		{
			while(true)
			{
				if (posInBuf >= buffSize)
					reload();

				if (suffix_number != buff[posInBuf])
					break;
				else
					++posInBuf;
			}
			return (buffPosInFile + posInBuf) & prefixMask;
		}

		~CPrefixFileBufferForListingMode()
		{
			delete[] buff;
		}
	};

	class OrderedBinReading {

		template<typename T>
		class buffered_scanning
		{
			FILEWrapper* file;
			std::vector<T> buffer;
			size_t buf_pos;
			size_t file_byte_pos;
			size_t file_byte_end_pos;

			bool reload()
			{
				assert(file_byte_end_pos >= file_byte_pos);

				auto left_in_file = file_byte_end_pos - file_byte_pos;

				assert(left_in_file % sizeof(T) == 0);

				if (!left_in_file)
					return false;

				auto bytes_to_read = buffer.size() * sizeof(T);

				if (left_in_file < bytes_to_read) {
					bytes_to_read = left_in_file;
					buffer.resize(bytes_to_read / sizeof(T));
				}

				file->Read(buffer.data(), sizeof(T), buffer.size());
				file_byte_pos += bytes_to_read;
				buf_pos = 0;
				return true;
			}
		public:
			void reset(FILEWrapper* file, size_t buff_size_bytes, size_t file_byte_pos, size_t file_byte_end_pos)
			{
				if (file_byte_end_pos - file_byte_pos < buff_size_bytes)
					buff_size_bytes = file_byte_end_pos - file_byte_pos;

				this->file = file;
				file->Seek(file_byte_pos);
				this->buffer.resize(buff_size_bytes / sizeof(T));
				this->buf_pos = buffer.size();
				this->file_byte_pos = file_byte_pos;
				this->file_byte_end_pos = file_byte_end_pos;
			}

			bool next_elem(T& out) {
				if (this->buf_pos == this->buffer.size())
					if (!this->reload())
						return false;
				out = this->buffer[this->buf_pos++];
				return true;
			}

			//caller responsible to call only when there is data
			//and to make how_many such that it will at some point land in buffer.size()
			T* read_many(size_t how_many) {
				assert(this->buf_pos <= this->buffer.size());
				if (this->buf_pos == this->buffer.size())
					if (!this->reload())
						assert(false);
				T* res = &this->buffer[this->buf_pos];
				this->buf_pos += how_many;
				return res;
			}
		};

		FILEWrapper& pre_file;
		size_t pre_file_data_start_pos;
		FILEWrapper& suf_file;
		size_t suf_file_data_start_pos;

		uint32_t guard;
		uint32_t cur_bin_id;
		uint32_t n_bins;
		uint32_t single_LUT_size;
		uint32_t suf_rec_size_bytes;
		std::vector<uint32_t> bin_map; //maps bin ids from kmc run to bin ids in kmc database

		std::vector<uint64_t> bin_sizes; //indexed with bin ids in kmc database, i.e. having bin_id from kmc run one should use bin_sizes[bin[map[id]], how many k-mers are in each bin
		uint64_t kmers_left_in_cur_bin{};

		//starting positions in file byte pos, guard at the end
		std::vector<size_t> bins_starts_in_pre;
		std::vector<size_t> bins_starts_in_suf;

		void calc_bin_ranges()
		{
			bins_starts_in_pre[0] = pre_file_data_start_pos;

			uint64_t x;
			pre_file.Seek(bins_starts_in_pre[0]);
			pre_file.Read(&x, sizeof(uint64), 1);

			bins_starts_in_suf[0] = suf_file_data_start_pos + suf_rec_size_bytes * x;

			for (size_t bin_id = 1; bin_id <= n_bins; ++bin_id)
			{
				bins_starts_in_pre[bin_id] = bins_starts_in_pre[bin_id - 1] + (single_LUT_size * sizeof(uint64_t));

				pre_file.Seek(bins_starts_in_pre[bin_id]);
				pre_file.Read(&x, sizeof(uint64), 1);

				bins_starts_in_suf[bin_id] = suf_file_data_start_pos + suf_rec_size_bytes * x;

				bin_sizes[bin_id - 1] = (bins_starts_in_suf[bin_id] - bins_starts_in_suf[bin_id - 1]) / suf_rec_size_bytes;
			}
		}

		uint64_t current_prefix;
		uint64_t last_val_from_prefix_file;
		uint64_t left_in_current_prefix;
		buffered_scanning<uint64_t> scan_prefix;
		buffered_scanning<unsigned char> scan_suffix;

		bool next_in_current_bin(uint32_t bin_id)
		{
			if (bin_id == guard)
				return false;

			bin_id = bin_map[bin_id];

			//if (!bin_sizes[bin_id])
			if (!kmers_left_in_cur_bin)
				return false;

			while (!left_in_current_prefix) {
				++current_prefix;

				uint64_t new_pre_val;
				bool have_elem = scan_prefix.next_elem(new_pre_val);
				assert(have_elem);
				left_in_current_prefix = new_pre_val - last_val_from_prefix_file;
				last_val_from_prefix_file = new_pre_val;
			}
			return true;
		}
	public:
		//for kmc signature selection scheme
		OrderedBinReading(FILEWrapper& pre_file, size_t pre_file_data_start_pos, FILEWrapper& suf_file, size_t suf_file_data_start_pos, const std::string& file_name, uint32_t n_bins, uint32 signature_len, uint32* signature_map, uint32 signature_map_size, uint32_t single_LUT_size, uint32_t suf_rec_size_bytes) :
			pre_file(pre_file),
			pre_file_data_start_pos(pre_file_data_start_pos),
			suf_file(suf_file),
			suf_file_data_start_pos(suf_file_data_start_pos),
			guard((uint32_t)-1/* std::numeric_limits<uint32_t>::max()*/), //VS have problems with numeric_limits::max(), probably defines max somewhere...
			cur_bin_id(guard),
			n_bins(n_bins),
			single_LUT_size(single_LUT_size),
			suf_rec_size_bytes(suf_rec_size_bytes),
			bin_map(n_bins, guard),
			bin_sizes(n_bins),
			bins_starts_in_pre(n_bins + 1),
			bins_starts_in_suf(n_bins + 1)
		{
			CSigToBinMap stbm(file_name);
			//stbm.Dump(std::cerr);
			if (stbm.GetNBins() != n_bins) {
				std::ostringstream oss;
				oss << "Number of bins (" << stbm.GetNBins() << ") in mapping file (" << file_name << ") is different than in KMC database (" << n_bins << ")\n";
				throw std::runtime_error(oss.str());
			}
			if (stbm.GetSigLen() != signature_len) {
				std::ostringstream oss;
				oss << "Signature length (" << stbm.GetSigLen() << ") in mapping file (" << file_name << ") is different than in KMC database (" << signature_len << ")\n";
				throw std::runtime_error(oss.str());
			}

			assert((uint32_t)stbm.GetMapping()[signature_map_size - 1] == n_bins - 1);

			for (uint32_t x = 0; x < signature_map_size; ++x)
			{
				if (stbm.GetMapping()[x] == -1) //disabled signature
					continue;
				if (bin_map[stbm.GetMapping()[x]] == guard)
					bin_map[stbm.GetMapping()[x]] = signature_map[x];
				else if (bin_map[stbm.GetMapping()[x]] != signature_map[x])
				{
					std::ostringstream oss;
					oss << "Signature to bin mapping from file " << file_name << " does not match signature to bin mapping in KMC database";
					throw std::runtime_error(oss.str());
				}
			}
			calc_bin_ranges();
		}

		//min_hash signature selection scheme
		OrderedBinReading(FILEWrapper& pre_file, size_t pre_file_data_start_pos, FILEWrapper& suf_file, size_t suf_file_data_start_pos, uint32_t n_bins, uint32_t single_LUT_size, uint32_t suf_rec_size_bytes, const std::vector<uint32_t>& bin_map) :
			pre_file(pre_file),
			pre_file_data_start_pos(pre_file_data_start_pos),
			suf_file(suf_file),
			suf_file_data_start_pos(suf_file_data_start_pos),
			guard((uint32_t)-1/* std::numeric_limits<uint32_t>::max()*/), //VS have problems with numeric_limits::max(), probably defines max somewhere...
			cur_bin_id(guard),
			n_bins(n_bins),
			single_LUT_size(single_LUT_size),
			suf_rec_size_bytes(suf_rec_size_bytes),
			bin_map(bin_map),
			bin_sizes(n_bins),
			bins_starts_in_pre(n_bins + 1),
			bins_starts_in_suf(n_bins + 1)
		{
			calc_bin_ranges();
		}


		uint32_t get_n_bins() const
		{
			return n_bins;
		}
		void start_bin(size_t bin_id, size_t prefix_file_buff_size_bytes, size_t suffix_file_buff_size_bytes)
		{
			cur_bin_id = bin_id;
			bin_id = bin_map[bin_id];

			current_prefix = (uint64_t)-1;// std::numeric_limits<uint64_t>::max(); //will be incremented
			left_in_current_prefix = 0;
			kmers_left_in_cur_bin = bin_sizes[bin_id];
			scan_prefix.reset(&pre_file, prefix_file_buff_size_bytes, bins_starts_in_pre[bin_id], bins_starts_in_pre[bin_id + 1] + sizeof(uint64_t));

			//we need to adjust suffix_file_buff_size_bytes such that it is divisible by suf_rec_size_bytes
			suffix_file_buff_size_bytes = suffix_file_buff_size_bytes / suf_rec_size_bytes * suf_rec_size_bytes;
			if (suffix_file_buff_size_bytes == 0)
				suffix_file_buff_size_bytes = suf_rec_size_bytes;

			scan_suffix.reset(&suf_file, suffix_file_buff_size_bytes, bins_starts_in_suf[bin_id], bins_starts_in_suf[bin_id + 1]);

			bool have_elem = scan_prefix.next_elem(last_val_from_prefix_file);
			assert(have_elem);
		}

		//go to the next bin if exists (even if empty)
		bool next_bin(size_t prefix_file_buff_size_bytes, size_t suffix_file_buff_size_bytes)
		{
			if (cur_bin_id == n_bins)
				return false;
			++cur_bin_id;
			if (cur_bin_id == n_bins)
				return false;

			start_bin(cur_bin_id, prefix_file_buff_size_bytes, suffix_file_buff_size_bytes);

			return true;
		}

		template<unsigned SIZE>
		bool read_from_cur_bin(CKmer<SIZE>& kmer, uint64& count, uint32_t lut_prefix_length, uint32_t suffix_size, uint32_t counter_size, uint32_t min_count, uint32_t max_count)
		{
			while (true)
			{
				if (!next_in_current_bin(cur_bin_id))
					return false;

				//--bin_sizes[bin_map[cur_bin_id]];
				--kmers_left_in_cur_bin;
				--left_in_current_prefix;
				auto suf_rec = scan_suffix.read_many(suf_rec_size_bytes);

				kmer.clear();

				uint32 a = 0;
				uint32_t kmer_byte_pos = suffix_size - 1;
				for (; a < suffix_size; ++a, --kmer_byte_pos)
					kmer.set_byte(kmer_byte_pos, suf_rec[a]);

				auto lut_prefix_len_bits = 2 * lut_prefix_length;
				auto prefix_bit_pos = suffix_size * 8;
				kmer.set_bits(prefix_bit_pos, lut_prefix_len_bits, current_prefix);

				//read counter:
				if (counter_size == 0) {
					count = 1;
					return true;
				}
				else
				{
					count = suf_rec[a++];
					for (uint32 b = 1; b < counter_size; b++)
					{
						uint64 aux = 0x000000ff & suf_rec[a++];
						aux = aux << 8 * (b);
						count = aux | count;
					}
				}

				if (counter_size == 0)
					return true;
				if (count >= min_count && count <= max_count)
					return true;
			}
			return false;
		}

		bool read_from_cur_bin(CKmerAPI& kmer, uint64& count, uint32_t off, uint32_t suffix_size, uint32_t counter_size, uint32_t min_count, uint32_t max_count)
		{
			while (true)
			{
				if (!next_in_current_bin(cur_bin_id))
					return false;
				
				//--bin_sizes[bin_map[cur_bin_id]];
				--kmers_left_in_cur_bin;
				--left_in_current_prefix;
				auto suf_rec = scan_suffix.read_many(suf_rec_size_bytes);

				uint64 temp_prefix = current_prefix << off;	// shift prefix towards MSD. "& prefix_mask" necessary for kmc2 db format

				kmer.kmer_data[0] = temp_prefix;			// store prefix in an object CKmerAPI

				for (uint32 i = 1; i < kmer.no_of_rows; i++)
					kmer.kmer_data[i] = 0;

				//read sufix:
				uint32 row_index = 0;
				uint64 suf = 0;

				off = off - 8;
				uint32 a = 0;
				for (; a < suffix_size; a++)
				{
					suf = suf_rec[a];

					suf = suf << off;
					kmer.kmer_data[row_index] = kmer.kmer_data[row_index] | suf;

					if (off == 0)				//the end of a word in kmer_data
					{
						off = 56;
						row_index++;
					}
					else
						off -= 8;
				}
				//read counter:
				if (counter_size == 0) {
					count = 1;
					return true;
				}
				else
				{
					count = suf_rec[a++];
					for (uint32 b = 1; b < counter_size; b++)
					{
						uint64 aux = 0x000000ff & suf_rec[a++];
						aux = aux << 8 * (b);
						count = aux | count;
					}
				}

				if (counter_size == 0)
					return true;
				if (count >= min_count && count <= max_count)
					return true;
			}
			return false;
		}

		uint64_t get_n_kmers(size_t bin_id) const
		{
			return bin_sizes[bin_map[bin_id]];
		}
	};
	std::unique_ptr<OrderedBinReading> ordered_bin_reading;
protected:
	enum open_mode {closed, opened_for_RA, opened_for_listing, opened_for_listing_with_bin_order};
	open_mode is_opened;
	uint64 suf_file_left_to_read = 0; // number of bytes that are yet to read in a listing mode
	uint64 suffix_file_total_to_read = 0; // number of bytes that constitutes records in kmc_suf file
	bool end_of_file;

	FILEWrapper file_pre;
	FILEWrapper file_suf;

	uint64* prefix_file_buf; //only for random access mode
	uint64 prefix_file_buf_size; //only for random access mode
	std::unique_ptr<CPrefixFileBufferForListingMode> prefixFileBufferForListingMode;

	uint64 prefix_index;			// The current prefix's index in an array "prefix_file_buf", readed from *.kmc_pre
	uint32 single_LUT_size{};			// The size of a single LUT (in no. of elements)

	uint32* signature_map;
	uint32_t n_bins{};				//only for KMC2
	std::vector<uint32_t> bins_order; //order of bins id in this file
	std::vector<uint32_t> bin_id_to_pos; //reverse relation to bin_order, bin_id_to_pos[i] is a position of bin $i$ in the kmc database
	uint32 signature_map_size{};
	
	uchar* sufix_file_buf;
	uint64 sufix_number;			// The sufix's number to be listed
	uint64 index_in_partial_buf;	// The current byte's number in an array "sufix_file_buf", for listing mode

	uint32 kmer_length;
	uint32 mode;
	uint32 counter_size;
	uint32 lut_prefix_length;
	uint32 signature_len;
	uint32 min_count;
	uint64 max_count;
	uint64 total_kmers;
	bool both_strands;
	KMC::SignatureSelectionScheme signature_selection_scheme;

	uint32 kmc_version;
	uint32 sufix_size;		// sufix's size in bytes 
	uint32 sufix_rec_size;  // sufix_size + counter_size

	uint32 original_min_count;
	uint64 original_max_count;

	static uint64 part_size; // the size of a block readed to sufix_file_buf, in listing mode 
	
	bool BinarySearch(int64 index_start, int64 index_stop, const CKmerAPI& kmer, uint64& counter, uint32 pattern_offset);

	// Open a file, recognize its size and check its marker. Auxiliary function.
	bool OpenASingleFile(const std::string &file_name, FILEWrapper &file_handler, uint64 &size, char marker[]);

	// Recognize current parameters. Auxiliary function.
	bool ReadParamsFrom_prefix_file_buf(uint64 &size, open_mode _open_mode, const std::string& bin_order_file_name = "");

	// Reload a contents of an array "sufix_file_buf" for listing mode. Auxiliary function. 
	void Reload_sufix_file_buf();

	// Implementation of GetCountersForRead for kmc1 database format for both strands
	bool GetCountersForRead_kmc1_both_strands(const std::string& read, std::vector<uint32>& counters);

	// Implementation of GetCountersForRead for kmc1 database format without choosing canonical k-mer
	bool GetCountersForRead_kmc1(const std::string& read, std::vector<uint32>& counters);		

	using super_kmers_t = std::vector<std::tuple<uint32, uint32, uint32>>;//start_pos, len, bin_no

	template<typename mmer_t>
	void GetSuperKmers(const std::string& transformed_read, super_kmers_t& super_kmers);

	// Implementation of GetCountersForRead for kmc2 database format for both strands
	bool GetCountersForRead_kmc2_both_strands(const std::string& read, std::vector<uint32>& counters);

	// Implementation of GetCountersForRead for kmc2 database format
	bool GetCountersForRead_kmc2(const std::string& read, std::vector<uint32>& counters);
public:
		
	CKMCFile(bool reopen_each_time = false);
	~CKMCFile();

	// Open files *.kmc_pre & *.kmc_suf, read them to RAM, close files. *.kmc_suf is opened for random access
	bool OpenForRA(const std::string &file_name);

	// Open files *kmc_pre & *.kmc_suf, read *.kmc_pre to RAM, *.kmc_suf is buffered
	bool OpenForListing(const std::string& file_name);

	//bin_order_file_name may be empty for min_hash signature selection scheme
	bool OpenForListingWithBinOrder(const std::string& file_name, const std::string& bin_order_file_name = "");

	//mkokot_TODO: I will use 0x201 for new format that I support different kind of minimizers
	// Return true if kmc is in KMC2 compatiblie format
	bool IsKMC2() const noexcept { return kmc_version == 0x201; }

	// Return next kmer in CKmerAPI &kmer. Return its counter in uint64 &count. Return true if not EOF
	bool ReadNextKmer(CKmerAPI &kmer, uint64 &count); //for small k-values when counter may be longer than 4bytes
	
	bool ReadNextKmer(CKmerAPI &kmer, uint32 &count);

	//only when oppened in listing with bin order
	uint32_t GetNBins() const;

	//only when oppened in listing with bin order
	bool StartBin(size_t prefix_file_buff_size_bytes = 1ull << 25 , size_t suffix_file_buff_size_bytes = 1ull << 25);

	//only when oppened in listing with bin order
	bool StartBin(uint32_t bin_id, size_t prefix_file_buff_size_bytes = 1ull << 25, size_t suffix_file_buff_size_bytes = 1ull << 25);

	//only when oppened in listing with bin order
	bool GetNKmers(uint32_t bin_id, uint64_t& n_kmers);

	//only when oppened in listing with bin order
	bool ReadNextKmerFromBin(CKmerAPI& kmer, uint64& count);

	//only when oppened in listing with bin order
	template<unsigned SIZE>
	bool ReadNextKmerFromBin(CKmer<SIZE>& kmer, uint64& count);

	// Release memory and close files in case they were opened 
	bool Close();

	// Set the minimal value for a counter. Kmers with counters below this theshold are ignored
	bool SetMinCount(uint32 x);

	// Return a value of min_count. Kmers with counters below this theshold are ignored 
	uint32 GetMinCount(void);

	// Set the maximal value for a counter. Kmers with counters above this theshold are ignored
	bool SetMaxCount(uint32 x);

	// Return a value of max_count. Kmers with counters above this theshold are ignored 
	uint64 GetMaxCount(void);
	
	//Return true if kmc was run without -b switch.
	bool GetBothStrands(void);

	// Return the total number of kmers between min_count and max_count
	uint64 KmerCount(void);

	// Return the length of kmers
	uint32 KmerLength(void);

	// Set initial values to enable listing kmers from the begining. Only in listing mode
	bool RestartListing(void);

	// Return true if all kmers are listed
	bool Eof(void);

	// Return true if kmer exists. In this case return kmer's counter in count
	bool CheckKmer(CKmerAPI &kmer, uint32 &count);

	bool CheckKmer(CKmerAPI &kmer, uint64 &count);

	// Return true if kmer exists
	bool IsKmer(CKmerAPI &kmer);

	// Set original (readed from *.kmer_pre) values for min_count and max_count
	void ResetMinMaxCounts(void);

	// Get current parameters from kmer_database
	bool Info(uint32 &_kmer_length, uint32 &_mode, uint32 &_counter_size, uint32 &_lut_prefix_length, uint32 &_signature_len, uint32 &_min_count, uint64 &_max_count, uint64 &_total_kmers);
	
	// Get current parameters from kmer_database
	bool Info(CKMCFileInfo& info);

	// Get counters for all k-mers in read
	bool GetCountersForRead(const std::string& read, std::vector<uint32>& counters);
	template<typename Runner>
	void DispatchKmerSize(Runner& runner);
private:
	uint32 count_for_kmer_kmc1(CKmerAPI& kmer);
	uint32 count_for_kmer_kmc2(CKmerAPI& kmer, uint32 bin_start_pos);
};

//-----------------------------------------------------------------------------------------------
template<unsigned SIZE>
bool CKMCFile::ReadNextKmerFromBin(CKmer<SIZE>& kmer, uint64& count)
{
	if (is_opened != opened_for_listing_with_bin_order)
		return false;

	return ordered_bin_reading->read_from_cur_bin(kmer, count, lut_prefix_length, sufix_size, counter_size, min_count, max_count);
}


template<unsigned SIZE, typename Runner>
struct KmerSizeDispatcher
{
	static void Dispatch(uint32_t kmer_len, Runner& runner)
	{
		auto min_k = 32 * (SIZE - 1);
		auto max_k = min_k + 32;
		if (kmer_len > min_k && kmer_len <= max_k)
			runner.template Run<SIZE>();
		else
			KmerSizeDispatcher<SIZE - 1, Runner>::Dispatch(kmer_len, runner);
	}
};

template<typename Runner>
struct KmerSizeDispatcher<0, Runner>
{
	static void Dispatch(uint32_t kmer_len, Runner& runner)
	{
		std::cerr << "Error: k-mer size dispatcher failed!\n";
		exit(1);
	}
};

template<typename Runner>
void DispatchKmerSize(uint32_t kmer_length, Runner& runner)
{
	constexpr auto max_no_of_uint64_t_for_kmer = (MAX_K + 31) / 32;

	if (kmer_length > MAX_K)
	{
		std::cerr << "Error: k too large, extend compilation constant MAX_K\n";
		exit(1);
	}
	KmerSizeDispatcher<max_no_of_uint64_t_for_kmer, Runner>::Dispatch(kmer_length, runner);
}

template<typename Runner>
void CKMCFile::DispatchKmerSize(Runner& runner)
{
	DispatchKmerSize<>(kmer_length, runner);
}


#endif

// ***** EOF
