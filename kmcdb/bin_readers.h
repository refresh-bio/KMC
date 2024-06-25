#ifndef BIN_READERS_H_
#define BIN_READERS_H_
#include "utils.h"
#include "stream_names.h"
#include "bin_metadata.h"

#include <vector>
#include <string>
#include <numeric>

namespace kmcdb
{
	namespace detail
	{
		template<unsigned SIZE, typename VALUE_T>
		bool check_kmer_binary_search(const CKmer<SIZE>& suffix, const uint8_t* begin, const uint8_t* end, VALUE_T* values, uint64_t num_values, const std::array<uint64_t, values_size<VALUE_T>()>& num_bytes_single_value, int32_t suffix_len, uint64_t single_suf_elem_bytes) noexcept
		{
			if (begin == end)
			{
				detail::SetZeros(values, num_values);
				return false;
			}
			end -= single_suf_elem_bytes;
			int32_t pos_a = (suffix_len + 3) / 4 - 1;

			std::ptrdiff_t begin_idx = 0;
			std::ptrdiff_t end_idx = (end - begin) / single_suf_elem_bytes;

			while (begin_idx <= end_idx)
			{
				std::ptrdiff_t mid_idx = (begin_idx + end_idx) / 2;
				const uint8_t* ptr = begin + mid_idx * single_suf_elem_bytes;
				int32_t pos = pos_a;

				uint8_t byte{};
				uint8_t x{};
				while (pos >= 0)
				{
					byte = suffix.get_byte(pos--);
					x = *ptr++;
					if (byte != x)
						break;
				}
				if (byte == x)
				{
					assert(pos + 1 == 0);
					detail::LoadValues(values, num_values, ptr, num_bytes_single_value);
					return true;
				}
				if (byte > x)
					begin_idx = mid_idx + 1;
				else
					end_idx = mid_idx - 1;
			}
			detail::SetZeros(values, num_values);
			return false;
		}

		template<typename VALUE_T>
		class BinReaderBase
		{
		protected:
			uint64_t bin_id;
			archive_input_t* archive;

			uint64_t kmer_len;
			uint64_t num_values;

			std::array<uint64_t, detail::values_size<VALUE_T>()> num_bytes_single_value;

			BinMetadata bin_metadata{};
		public:
			BinReaderBase(uint64_t bin_id,
				archive_input_t* archive,
				uint64_t kmer_len,
				uint64_t num_values,
				const std::vector<uint64_t>& num_bytes_single_value) :
				bin_id(bin_id),
				archive(archive),
				kmer_len(kmer_len),
				num_values(num_values),
				num_bytes_single_value(detail::vec_to_array<uint64_t, detail::values_size<VALUE_T>()>(num_bytes_single_value))
			{
				const int metadata_stream_id = archive->get_stream_id(stream_names::BinMetadata(bin_id));
				if (metadata_stream_id == -1)
					throw std::runtime_error("Cannot find " + stream_names::BinMetadata(bin_id) + " stream");

				std::vector<uint8_t> serialized_bin_metadata;
				size_t tmp;
				if (!archive->get_part(metadata_stream_id, serialized_bin_metadata, tmp))
					throw std::runtime_error("Stream" + stream_names::BinMetadata(bin_id) + " is empty");

				bin_metadata.load(serialized_bin_metadata);

				if (archive->get_part(metadata_stream_id, serialized_bin_metadata, tmp))
					throw std::runtime_error("Unexpected data in stream" + stream_names::BinMetadata(bin_id));
			}

			const BinMetadata& GetBinMetadata() const
			{
				return this->bin_metadata;
			}
		};

		template<typename VALUE_T>
		class BinReaderSortedPlainBase : public BinReaderBase<VALUE_T>
		{
		protected:

			int stream_id;
			uint64_t bytes_for_kmer;
			uint64_t single_elem_bytes;

		public:
			BinReaderSortedPlainBase(const BinReaderSortedPlainBase&) = delete;
			BinReaderSortedPlainBase(BinReaderSortedPlainBase&&) = delete;
			const BinReaderSortedPlainBase& operator=(const BinReaderSortedPlainBase&) = delete;
			const BinReaderSortedPlainBase& operator=(BinReaderSortedPlainBase&&) = delete;
			~BinReaderSortedPlainBase() = default;

			BinReaderSortedPlainBase(uint64_t bin_id,
				archive_input_t* archive,
				uint64_t kmer_len,
				uint64_t num_values,
				const std::vector<uint64_t>& num_bytes_single_value) :
				detail::BinReaderBase<VALUE_T>(bin_id, archive, kmer_len, num_values, num_bytes_single_value),
				stream_id(archive->get_stream_id(stream_names::Bin(bin_id))),
				bytes_for_kmer((kmer_len + 3) / 4),
				single_elem_bytes(bytes_for_kmer + std::accumulate(num_bytes_single_value.begin(), num_bytes_single_value.end(), 0ull) * num_values)
			{
				if (stream_id == -1)
					throw std::runtime_error("Cannot find " + stream_names::Bin(bin_id) + " stream");
			}
		};

		template<typename VALUE_T>
		class BinReaderSortedWithLUTBase : public BinReaderBase<VALUE_T>
		{
		protected:
			int stream_id_suf;
			int stream_id_LUT;

			uint64_t lut_prefix_len;
			uint64_t bytes_for_kmer_suffix;
			uint64_t single_suf_elem_bytes;

		public:
			BinReaderSortedWithLUTBase(uint64_t bin_id,
				archive_input_t* archive,
				uint64_t kmer_len,
				uint64_t num_values,
				uint64_t lut_prefix_len,
				const std::vector<uint64_t>& num_bytes_single_value) :
				detail::BinReaderBase<VALUE_T>(bin_id, archive, kmer_len, num_values, num_bytes_single_value),
				stream_id_suf(archive->get_stream_id(stream_names::BinSufData(bin_id))),
				stream_id_LUT(archive->get_stream_id(stream_names::BinLut(bin_id))),
				lut_prefix_len(lut_prefix_len),
				bytes_for_kmer_suffix((kmer_len - lut_prefix_len + 3) / 4),
				single_suf_elem_bytes(bytes_for_kmer_suffix + std::accumulate(num_bytes_single_value.begin(), num_bytes_single_value.end(), 0ull) * num_values)
			{
				//mkokot_TODO: there is some code duplication here -> refactor
				if (stream_id_suf == -1)
					throw std::runtime_error("Cannot find " + stream_names::BinSufData(bin_id) + " stream");

				if (stream_id_LUT == -1)
					throw std::runtime_error("Cannot find " + stream_names::BinLut(bin_id) + " stream");
			}
		};
	}

	template<typename VALUE_T>
	class BinReaderSortedPlainForListing : public detail::BinReaderSortedPlainBase<VALUE_T>
	{
		struct current_pack_t
		{
			std::vector<uint8_t> data{};
			size_t read_pos{};
		};
		current_pack_t current_pack{};

		size_t max_single_read_size;

		//such that each part contains full records
		//mkokot_TODO: this is code repetition -> fix it!
		static size_t adjust_max_part_size(size_t val, size_t single_elem_bytes)
		{
			if (val < single_elem_bytes)
				return single_elem_bytes;

			if (single_elem_bytes == 0)
				return val;

			return val / single_elem_bytes * single_elem_bytes;
		}

		using detail::BinReaderSortedPlainBase<VALUE_T>::archive;
		using detail::BinReaderSortedPlainBase<VALUE_T>::stream_id;
		using detail::BinReaderSortedPlainBase<VALUE_T>::single_elem_bytes;
		using detail::BinReaderSortedPlainBase<VALUE_T>::kmer_len;
		using detail::BinReaderSortedPlainBase<VALUE_T>::num_values;
		using detail::BinReaderSortedPlainBase<VALUE_T>::num_bytes_single_value;
	public:
		BinReaderSortedPlainForListing(uint64_t bin_id,
			archive_input_t* archive,
			uint64_t kmer_len,
			uint64_t num_values,
			const std::vector<uint64_t>& num_bytes_single_value,
			size_t max_single_read_size) :
			detail::BinReaderSortedPlainBase<VALUE_T>(bin_id, archive, kmer_len, num_values, num_bytes_single_value),
			max_single_read_size(adjust_max_part_size(max_single_read_size, single_elem_bytes))
		{

		}

		template<unsigned SIZE>
		bool NextKmer(CKmer<SIZE>& kmer, VALUE_T* values)
		{
			if (current_pack.data.size() == current_pack.read_pos)
			{
				int part_id;
				int sub_part_id;
				uint64_t meta;
				if (!archive->get_sub_part(stream_id, max_single_read_size, part_id, sub_part_id, current_pack.data, meta))
					return false;
				current_pack.read_pos = 0;
			}
			const uint8_t* ptr = current_pack.data.data() + current_pack.read_pos;
			kmer.load_from_left_aligned(ptr, kmer_len);
			current_pack.read_pos += ptr - (current_pack.data.data() + current_pack.read_pos);

			detail::LoadValues(values, num_values, current_pack.data, current_pack.read_pos, this->num_bytes_single_value);

			return true;
		}
	};
	//mkokot_TODO: tak swoja droga to moze w ogole nie ma sensu w zapisie binarnym robic podzialu na lut i suffix tylko zrobic zawsze ta kompresje z l. wspolnych symboli z nastepnym k-merem
	//nie wiem czy to bedzie lepsze czy nie, ale mialoby ta zalete ze czasem i tak nie jestesmy w stanie okreslic jaki lut jest optymalny przy zapisie bo nie wiemy z gory ile bedzie zapisanych k-merow
	//natomiast ten podzial na lut dla random access to spoko mozna by zrobic bo wtedy wiemy ile jest k-merow, a i tak musimy dane odczytac i zdeserializowac
	//dobra na razie robie tak, ze po prostu bedzie binary searchem w surowej pamieci to robione

	//This is for sorted representation
	template<typename VALUE_T>
	class BinReaderSortedPlainForRandomAccess : public detail::BinReaderSortedPlainBase<VALUE_T>
	{
		std::unique_ptr<uint8_t[]> flat_data;
	public:
		BinReaderSortedPlainForRandomAccess(uint64_t bin_id,
			archive_input_t* archive,
			uint64_t kmer_len,
			uint64_t num_values,
			const std::vector<uint64_t>& num_bytes_single_value) :
			detail::BinReaderSortedPlainBase<VALUE_T>(bin_id, archive, kmer_len, num_values, num_bytes_single_value),
			flat_data(std::make_unique_for_overwrite<uint8_t[]>(this->bin_metadata.total_kmers* this->single_elem_bytes))
		{
			if (!flat_data)
				throw std::runtime_error("Cannot allocate memory for bin random access");
			uint8_t* ptr = flat_data.get();
			size_t size{};
			uint64_t metadata;
			while (archive->get_part(this->stream_id, ptr, size, metadata))
				ptr += size;

			assert(ptr - flat_data.get() == static_cast<std::ptrdiff_t>(this->bin_metadata.total_kmers * this->single_elem_bytes));
		}

		template<unsigned SIZE>
		bool CheckKmer(const CKmer<SIZE>& kmer, VALUE_T* values) const noexcept
		{
			auto start_ptr = flat_data.get();
			auto end_ptr = flat_data.get() + this->bin_metadata.total_kmers * this->single_elem_bytes;

			//mkokot_TODO: since as for now my binary representation of CKmer is right aligned i need to SHL sometimes
			auto copy = kmer;
			if (this->kmer_len % 4)
				copy.SHL(4 - this->kmer_len % 4);

			//but when this will change the code below should just work and the `copy` variable and SHL is not needed
			//return detail::check_kmer_binary_search(kmer, start_ptr, end_ptr, values, this->num_values, this->num_bytes_single_value, static_cast<int32_t>(this->kmer_len), this->single_elem_bytes);

			return detail::check_kmer_binary_search(copy, start_ptr, end_ptr, values, this->num_values, this->num_bytes_single_value, static_cast<int32_t>(this->kmer_len), this->single_elem_bytes);
		}
	};

	template<typename VALUE_T>
	class BinReaderSortedWithLUTForListing : public detail::BinReaderSortedWithLUTBase<VALUE_T>
	{
		uint64_t current_prefix{};
		uint64_t already_readed_kmers{};
		struct current_pack_lut_t
		{
			std::vector<uint64_t> data{};
			size_t read_pos{};
			void resize(size_t new_size)
			{
				data.resize(new_size);
			}
		};
		current_pack_lut_t current_pack_lut{};
		struct current_pack_suf_t
		{
			std::vector<uint8_t> data{};
			size_t read_pos{};
		};
		current_pack_suf_t current_pack_suf{};

		size_t max_single_read_size_suf;
		size_t max_single_read_size_lut;

		//such that each part contains full records
		//mkokot_TODO: this is code repetition -> fix it!
		static size_t adjust_max_part_size(size_t val, size_t single_elem_bytes)
		{
			if (val < single_elem_bytes)
				return single_elem_bytes;

			if (single_elem_bytes == 0)
				return val;

			return val / single_elem_bytes * single_elem_bytes;
		}

		using detail::BinReaderSortedWithLUTBase<VALUE_T>::archive;

		using detail::BinReaderSortedWithLUTBase<VALUE_T>::stream_id_suf;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::stream_id_LUT;

		using detail::BinReaderSortedWithLUTBase<VALUE_T>::kmer_len;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::num_values;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::lut_prefix_len;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::bytes_for_kmer_suffix;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::num_bytes_single_value;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::single_suf_elem_bytes;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::bin_metadata;

		bool read_lut_pack()
		{
			int part_id;
			int sub_part_id;
			uint64_t meta;
			assert(max_single_read_size_lut % sizeof(uint64_t) == 0);
			current_pack_lut.resize(max_single_read_size_lut / sizeof(uint64_t));

			uint8_t* ptr = reinterpret_cast<uint8_t*>(current_pack_lut.data.data());
			size_t readed{};
			if (!archive->get_sub_part(stream_id_LUT, max_single_read_size_lut, part_id, sub_part_id, ptr, readed, meta))
				return false;

			assert(readed % sizeof(uint64_t) == 0);
			assert(readed <= max_single_read_size_lut);

			current_pack_lut.resize(readed / sizeof(uint64_t));

			if constexpr (std::endian::native == std::endian::little)
				;
			else if constexpr (std::endian::native == std::endian::big)
				for (uint64_t& lut_elem : current_pack_lut.data)
					lut_elem = refresh::serialization::bswap_uint64(lut_elem);
			else
				static_assert(!sizeof(VALUE_T), "Unknown endian");

			return true;
		}

	public:
		BinReaderSortedWithLUTForListing(uint64_t bin_id,
			archive_input_t* archive,
			uint64_t kmer_len,
			uint64_t num_values,
			uint64_t lut_prefix_len,
			const std::vector<uint64_t>& num_bytes_single_value,
			size_t max_single_read_size) :
			detail::BinReaderSortedWithLUTBase<VALUE_T>(
				bin_id,
				archive,
				kmer_len,
				num_values,
				lut_prefix_len,
				num_bytes_single_value),
			max_single_read_size_suf(adjust_max_part_size(max_single_read_size, single_suf_elem_bytes)),
			max_single_read_size_lut(adjust_max_part_size(max_single_read_size, 2 * sizeof(uint64_t))) //we need at least two uint64_t, if subpart returns less it will not work (although could be possibly implemented calling subpart many times)
		{
			//read first lut part
			if (!read_lut_pack())
				throw std::runtime_error("Cannot read LUT from stream " +
					stream_names::BinMetadata(bin_id));

			//LUT begins with zero
			++current_pack_lut.read_pos;
		}

		template<unsigned SIZE>
		bool NextKmer(CKmer<SIZE>& kmer, VALUE_T* values)
		{
			if (already_readed_kmers == bin_metadata.total_kmers)
				return false;
			while (true)
			{
				if (already_readed_kmers < current_pack_lut.data[current_pack_lut.read_pos])
					break;
				++current_pack_lut.read_pos;
				++current_prefix;
				if (current_pack_lut.data.size() == current_pack_lut.read_pos)
				{
					//size_t tmp;
					if (!read_lut_pack())
					{
						assert(already_readed_kmers == bin_metadata.total_kmers);
						return false;
					}
					current_pack_lut.read_pos = 0;
				}
			}

			if (current_pack_suf.data.size() == current_pack_suf.read_pos && single_suf_elem_bytes != 0)
			{
				int part_id;
				int sub_part_id;
				uint64_t meta;
				if (!archive->get_sub_part(stream_id_suf, max_single_read_size_suf, part_id, sub_part_id, current_pack_suf.data, meta))
					assert(false); //this should never happen because return false will be due to lut reading

				current_pack_suf.read_pos = 0;
			}

			const uint8_t* ptr = current_pack_suf.data.data() + current_pack_suf.read_pos;
			kmer.load_from_left_aligned(ptr, kmer_len - lut_prefix_len);
			current_pack_suf.read_pos += ptr - (current_pack_suf.data.data() + current_pack_suf.read_pos);

			kmer.set_prefix(current_prefix, static_cast<uint32_t>((kmer_len - lut_prefix_len) * 2));
			detail::LoadValues(values, num_values, current_pack_suf.data, current_pack_suf.read_pos, this->num_bytes_single_value);
			
			++already_readed_kmers;
			return true;
		}

		const BinMetadata& GetBinMetadata() const
		{
			return bin_metadata;
		}
	};

	template<typename VALUE_T>
	class BinReaderSortedWithLUTForRandomAccess : public detail::BinReaderSortedWithLUTBase<VALUE_T>
	{
		std::unique_ptr<uint8_t[]> flat_data;
		std::unique_ptr<uint64_t[]> LUT;

		using detail::BinReaderSortedWithLUTBase<VALUE_T>::archive;

		using detail::BinReaderSortedWithLUTBase<VALUE_T>::stream_id_suf;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::stream_id_LUT;

		using detail::BinReaderSortedWithLUTBase<VALUE_T>::kmer_len;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::num_values;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::lut_prefix_len;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::bytes_for_kmer_suffix;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::num_bytes_single_value;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::single_suf_elem_bytes;
		using detail::BinReaderSortedWithLUTBase<VALUE_T>::bin_metadata;
	public:
		BinReaderSortedWithLUTForRandomAccess(uint64_t bin_id,
			archive_input_t* archive,
			uint64_t kmer_len,
			uint64_t num_values,
			uint64_t lut_prefix_len,
			const std::vector<uint64_t>& num_bytes_single_value) :
			detail::BinReaderSortedWithLUTBase<VALUE_T>(
				bin_id,
				archive,
				kmer_len,
				num_values,
				lut_prefix_len,
				num_bytes_single_value),
			flat_data(std::make_unique_for_overwrite<uint8_t[]>(
				this->bin_metadata.total_kmers* this->single_suf_elem_bytes)),
			LUT(std::make_unique_for_overwrite<uint64_t[]>(
				(1ull << (2 * lut_prefix_len)) + 1))
		{
			//mkokot_TODO: te sprawdzenia bez sensu bo make_unique rzuci bad alloc!!!!
			if (!flat_data)
				throw std::runtime_error("Cannot allocate memory for bin random access (suffix)");

			if (!LUT)
				throw std::runtime_error("Cannot allocate memory for bin random access (LUT)");

			uint8_t* ptr = flat_data.get();
			size_t size{};
			uint64_t metadata;
			while (archive->get_part(stream_id_suf, ptr, size, metadata))
				ptr += size;

			assert(ptr - flat_data.get() == static_cast<std::ptrdiff_t>(
				this->bin_metadata.total_kmers * this->single_suf_elem_bytes));

			//load LUT
			size_t lut_elems = (1ull << (2 * lut_prefix_len)) + 1; //+1 guard

			size_t size_loaded{};
			while (true)
			{
				uint64_t tmp;
				size_t size_loaded_part;
				if (!archive->get_part(stream_id_LUT, reinterpret_cast<uint8_t*>(
						LUT.get()) + size_loaded, size_loaded_part, tmp))
					break;

				size_loaded += size_loaded_part;
			}

			assert(size_loaded == lut_elems * sizeof(uint64_t));

			if constexpr (std::endian::native == std::endian::little)
				;
			else if constexpr (std::endian::native == std::endian::big)
				for (size_t i = 0; i < lut_elems; ++i)
					refresh::serialization::bswap_uint64(LUT[i]);
			else
				static_assert(!sizeof(VALUE_T), "Unknown endian");
		}

		template<unsigned SIZE>
		bool CheckKmer(const CKmer<SIZE>& kmer, VALUE_T* values) const noexcept
		{
			uint64_t prefix = 0;
			if (lut_prefix_len)
				prefix = kmer.remove_suffix(static_cast<uint32_t>
					(2 * (kmer_len - lut_prefix_len)));

			assert(prefix <= ((1ull << (2 * lut_prefix_len)) + 1));

			//whole k-mer is in lut
			if (this->kmer_len == this->lut_prefix_len)
			{
				if (LUT[prefix] == LUT[prefix + 1])
				{
					detail::SetZeros(values, this->num_values);
					return false;
				}
				//found
				const uint8_t* ptr = flat_data.get() + single_suf_elem_bytes * LUT[prefix];
				detail::LoadValues(values, this->num_values, ptr, this->num_bytes_single_value);
				return true;
			}

			auto start_ptr = flat_data.get() + single_suf_elem_bytes * LUT[prefix];
			auto end_ptr = flat_data.get() + single_suf_elem_bytes * LUT[prefix + 1];

			CKmer<SIZE> mask;
			mask.set_n_1(static_cast<uint32_t>(2 * (kmer_len - lut_prefix_len)));
			CKmer<SIZE> suffix = kmer;
			suffix.mask(mask);

			return detail::check_kmer_binary_search(suffix, start_ptr, end_ptr, values, this->num_values, this->num_bytes_single_value, static_cast<int32_t>(kmer_len - lut_prefix_len), single_suf_elem_bytes);
		}
	};
}
#endif // ! BIN_READERS_H_
