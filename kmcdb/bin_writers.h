#ifndef BIN_WRITERS_H_
#define BIN_WRITERS_H_

#include "utils.h"
#include "bin_metadata.h"
#include <numeric>

namespace kmcdb
{
	//simple bin writer, no lut, no unitigs, just simple k-mers
	template<typename VALUE_T>
	class BinWriterSortedPlain
	{
		uint64_t bin_id;
		archive_t* archive;
		int stream_id;

		uint64_t kmer_len;
		uint64_t num_values;
		uint64_t bytes_for_kmer;
		std::array<uint64_t, detail::values_size<VALUE_T>()> num_bytes_single_value;
		uint64_t single_elem_bytes;

		detail::serialization_buffer serialized_data;

		BinMetadata bin_metadata{};

		bool is_closed = false;
	public:

		//mkokot_TODO: this should be never called from outside, only created in Writer and other classes
		//consider some techiniques to hide this (private and friend is complex because of the fact make_unique creates object...)
		BinWriterSortedPlain(
			uint64_t bin_id,
			archive_t* archive,
			uint64_t kmer_len,
			uint64_t num_values,
			const std::vector<uint64_t>& num_bytes_single_value,
			size_t max_part_size) :
			bin_id(bin_id),
			archive(archive),
			stream_id(archive->register_stream(stream_names::Bin(bin_id))),
			kmer_len(kmer_len),
			num_values(num_values),
			bytes_for_kmer((kmer_len + 3) / 4),
			num_bytes_single_value(detail::vec_to_array<uint64_t, detail::values_size<VALUE_T>()>(num_bytes_single_value)),
			single_elem_bytes(bytes_for_kmer + std::accumulate(num_bytes_single_value.begin(), num_bytes_single_value.end(), 0ull) * num_values),
			serialized_data((std::max)(max_part_size, single_elem_bytes))
		{

		}

		//must add unique k-mers and in sorted order
		template<unsigned SIZE>
		void AddKmer(const CKmer<SIZE>& kmer, const VALUE_T* values)
		{
			assert(!is_closed);
			if (single_elem_bytes + serialized_data.size() > serialized_data.capacity())
			{
				assert(serialized_data.size());
				auto res = archive->add_part(stream_id, serialized_data.get(), serialized_data.size());
				if (!res)
					throw std::runtime_error("Cannot store k-mers");

				serialized_data.clear();

				assert(single_elem_bytes + serialized_data.size() <= serialized_data.capacity());
			}

			kmer.store_left_aligned(serialized_data.get_end(), kmer_len);

			detail::SerializeValues(values, num_values, serialized_data.get_end(), num_bytes_single_value);

			++bin_metadata.total_kmers;
		}

		void Close()
		{
			if (is_closed)
				return;
			is_closed = true;

			if (serialized_data.size())
			{
				auto res = archive->add_part(stream_id, serialized_data.get(), serialized_data.size());
				if (!res)
					throw std::runtime_error("Cannot store k-mers");

				serialized_data.clear();
			}
			std::vector<uint8_t> serialized_bin_metadata;
			bin_metadata.serialize(serialized_bin_metadata);
			auto res = archive->add_part(archive->register_stream(stream_names::BinMetadata(bin_id)), serialized_bin_metadata);
			if (!res)
				throw std::runtime_error("Cannot store bin metadata");
		}

		//Scott Meyers, Effective C++, Item 11: Prevent exceptions from leaving destructors.
		~BinWriterSortedPlain() noexcept
		{
			try
			{
				Close();
			}
			catch (...) {}
		}
	};

	template<typename VALUE_T>
	class BinWriterSortedWithLUT;


	template<typename VALUE_T>
	class BinWriterSortedWithLUTRaw
	{
		uint64_t bin_id;
		archive_t* archive;
		int stream_id_suf;
		int stream_id_LUT;

		uint64_t kmer_len;
		uint64_t num_values;
		uint64_t lut_prefix_len;
		uint64_t bytes_for_kmer_suffix;
		std::array<uint64_t, detail::values_size<VALUE_T>()> num_bytes_single_value;
		uint64_t single_suf_elem_bytes;
		BinMetadata bin_metadata{};

		friend class BinWriterSortedWithLUT<VALUE_T>;

		bool is_closed = false;
		bool any_write_started = false; //mkokot_TODO: this is for error detection purposes only
	public:
		//mkokot_TODO: this should be never called from outside, only created in Writer and other classes
		//consider some techniques to hide this (private and friend is complex because of the fact make_unique creates object...)
		BinWriterSortedWithLUTRaw(
			uint64_t bin_id,
			archive_t* archive,
			uint64_t kmer_len,
			uint64_t num_values,
			uint64_t lut_prefix_len,
			const std::vector<uint64_t>& num_bytes_single_value) :
			bin_id(bin_id),
			archive(archive),
			stream_id_suf(archive->register_stream(stream_names::BinSufData(bin_id))),
			stream_id_LUT(archive->register_stream(stream_names::BinLut(bin_id))),
			kmer_len(kmer_len),
			num_values(num_values),
			num_bytes_single_value(detail::vec_to_array<uint64_t, detail::values_size<VALUE_T>()>(num_bytes_single_value))
		{
			ChangeLutPrefixLen(lut_prefix_len);
		}

		void AddSufAndData(const uint8_t* begin, const uint8_t* end)
		{
			any_write_started = true;
			assert(!is_closed);
			assert((end - begin) % single_suf_elem_bytes == 0);
			bin_metadata.total_kmers += (end - begin) / single_suf_elem_bytes;

			auto res = archive->add_part(stream_id_suf, begin, end - begin);
			if (!res)
				throw std::runtime_error("Cannot store k-mer suffix + data");
		}

		//this should be called when all suffix data was added
		void AddLUT(uint64_t* LUT)
		{
			any_write_started = true;
			assert(!is_closed);
			uint64_t lut_recs = 1ull << (2 * lut_prefix_len);

			if (single_suf_elem_bytes == 0)
				bin_metadata.total_kmers = LUT[lut_recs];

			assert(LUT[lut_recs] == bin_metadata.total_kmers);

			if constexpr (std::endian::native == std::endian::big)
				for (size_t pos = 0; pos < lut_recs + 1; ++pos)
					LUT[pos] = refresh::serialization::bswap_uint64(LUT[pos]);
			else if constexpr (std::endian::native == std::endian::little)
				;//if the machine is little endian we may just cast it to uint_8*
			else
				static_assert(!sizeof(VALUE_T), "Unknown endian");

			//now I can just get the raw pointer, and the endianness is fine
			uint8_t* ptr = reinterpret_cast<uint8_t*>(LUT);

			auto res = archive->add_part(stream_id_LUT, ptr, (lut_recs + 1) * sizeof(uint64_t));
			if (!res)
				throw std::runtime_error("Cannot store k-mer LUT");
		}

		void Close()
		{
			if (is_closed)
				return;
			is_closed = true;

			std::vector<uint8_t> serialized_bin_metadata;
			bin_metadata.serialize(serialized_bin_metadata);
			bool res = archive->add_part(archive->register_stream(stream_names::BinMetadata(bin_id)), serialized_bin_metadata);
			if (!res)
				throw std::runtime_error("Cannot store bin metadata");
		}

		void ChangeLutPrefixLen(uint64_t new_lut_prefix_len)
		{
			if (any_write_started)
				throw std::runtime_error("Cannot change lut prefix len after writes started");

			lut_prefix_len = new_lut_prefix_len;
			bytes_for_kmer_suffix = (kmer_len - lut_prefix_len + 3) / 4;
			single_suf_elem_bytes = bytes_for_kmer_suffix + std::accumulate(num_bytes_single_value.begin(), num_bytes_single_value.end(), 0ull) * num_values;

		}

		~BinWriterSortedWithLUTRaw() noexcept
		{
			try
			{
				Close();
			}
			catch (...) {}
		}
	};

	//mkokot_TODO: add common base class for this and Raw ?
	//this variant just provides AddKmer
	template<typename VALUE_T>
	class BinWriterSortedWithLUT
	{
		BinWriterSortedWithLUTRaw<VALUE_T> impl;

		std::vector<uint64_t> LUT;

		detail::serialization_buffer serialized_data;

		bool is_closed = false;
	public:

		//mkokot_TODO: this should be never called from outside, only created in Writer and other classes
		//consider some techniques to hide this (private and friend is complex because of the fact make_unique creates object...)
		//mkokot_TODO:I just started this class, a lot to be done here
		BinWriterSortedWithLUT(
			uint64_t bin_id,
			archive_t* archive,
			uint64_t kmer_len,
			uint64_t num_values,
			uint64_t lut_prefix_len,
			const std::vector<uint64_t>& num_bytes_single_value,
			size_t max_part_size) :
			impl(bin_id,
				archive,
				kmer_len,
				num_values,
				lut_prefix_len,
				num_bytes_single_value),
			LUT((1ull << (2 * lut_prefix_len)) + 1), // + 1 for guard
			serialized_data((std::max)(max_part_size, impl.single_suf_elem_bytes))
		{

		}

		//must add unique k-mers and in sorted order
		template<unsigned SIZE>
		void AddKmer(const CKmer<SIZE>& kmer, const VALUE_T* values)
		{
			assert(!is_closed);
			//mkokot_TODO: probably some changes will be needed here, when we change the k-mer representation to left aligned
			uint64_t prefix = 0;

			//if its 0, there is an issue with shift for example for k=32
			if (impl.lut_prefix_len)
				prefix = kmer.remove_suffix(static_cast<uint32_t>(2 * (impl.kmer_len - impl.lut_prefix_len)));

			++LUT[prefix];

			if (impl.single_suf_elem_bytes + serialized_data.size() > serialized_data.capacity())
			{
				assert(serialized_data.size());

				impl.AddSufAndData(serialized_data.get(), serialized_data.get() + serialized_data.size());
				serialized_data.clear();

				assert(impl.single_suf_elem_bytes + serialized_data.size() <= serialized_data.capacity());
			}

			kmer.store_left_aligned(serialized_data.get_end(), impl.kmer_len - impl.lut_prefix_len);

			detail::SerializeValues(values, impl.num_values, serialized_data.get_end(), impl.num_bytes_single_value);
		}

		void Close()
		{
			if (is_closed)
				return;
			is_closed = true;

			if (serialized_data.size())
				impl.AddSufAndData(serialized_data.get(), serialized_data.get() + serialized_data.size());

			uint64_t lut_recs = 1ull << (2 * impl.lut_prefix_len);
			uint64_t prev = 0;
			for (uint64_t i = 0; i < lut_recs; ++i)
			{
				auto x = LUT[i];
				LUT[i] = prev;
				prev += x;
			}
			LUT[lut_recs] = prev; //guard
			impl.AddLUT(LUT.data());
			impl.Close();
		}

		//Scott Meyers, Effective C++, Item 11: Prevent exceptions from leaving destructors.
		//but maybe I should consider this https://stackoverflow.com/questions/130117/if-you-shouldnt-throw-exceptions-in-a-destructor-how-do-you-handle-errors-in-i
		//and especially https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4152.pdf
		//https://en.cppreference.com/w/cpp/error/uncaught_exception
		//but then need to remember to noexcept(false)
		~BinWriterSortedWithLUT() noexcept
		{
			try
			{
				Close();
			}
			catch (...) {}
		}
	};

}
#endif // ! BIN_WRITERS_H_