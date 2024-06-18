#ifndef READERS_H_
#define READERS_H_

#include "utils.h"
#include "stream_names.h"
#include "metadata.h"
#include "bin_readers.h"
#include "mmer.h"
#include "kmer_to_bin_id.h"
#include "metadata_reader.h"
#include <array>

namespace kmcdb
{
	namespace detail
	{
		inline archive_t* get_archive_from_metadata_reader(MetadataReader& metadata_reader)
		{
			return &metadata_reader.archive;
		}
		inline Metadata& get_metadata_from_metadata_reader(MetadataReader& metadata_reader)
		{
			return metadata_reader.metadata;
		}

		class HistoryReader
		{
			archive_t* archive;
			int stream_id;
		public:
			HistoryReader(MetadataReader& metadata_reader) :
			archive(detail::get_archive_from_metadata_reader(metadata_reader)),
			stream_id(archive->get_stream_id(stream_names::HISTORY))
			{
				if (stream_id == -1)
					throw std::runtime_error("Cannot find " + stream_names::HISTORY + " stream");
			}

			size_t GetNumHistoryItems() const
			{
				return archive->get_no_parts(stream_id);
			}

			bool Ith(HistoryItem& history_item, int i)
			{
				std::vector<uint8_t> serialized;
				uint64_t meta;
				if (!archive->get_part(stream_id, i, serialized, meta))
					return false;
				history_item.load(serialized);
				return false;
			}

			bool Next(HistoryItem& history_item)
			{
				std::vector<uint8_t> serialized;
				uint64_t meta;
				if (!archive->get_part(stream_id, serialized, meta))
					return false;

				history_item.load(serialized);
				return true;
			}

			bool Last(HistoryItem& history_item)
			{
				return Ith(history_item, static_cast<int>(GetNumHistoryItems()) - 1);
			}

			void Reset()
			{
				archive->rewind(stream_id);
			}

		};

		template<typename VALUE_T>
		class ReaderBase
		{
		protected:
			const Metadata& metadata;
			archive_t* archive;
		public:
			ReaderBase(MetadataReader& metadata_reader) :
				metadata(get_metadata_from_metadata_reader(metadata_reader)),
				archive(get_archive_from_metadata_reader(metadata_reader))
			{
				//mkokot_TODO: different exception type?
				if (!AreSame<VALUE_T>(metadata.value_types))
					throw std::runtime_error("Wrong template instantiation for " + detail::to_string(metadata.value_types));

				detail::IterateValues(VALUE_T{}, [this]<typename T>(auto idx, T& /*val*/)
				{
					(void)idx; //suppress warning unused (idx is used in assert, so not always depending on NDEBUG)
					if constexpr (!std::is_integral_v<T>)
						assert(metadata.config.num_bytes_single_value[idx] == sizeof(T));
				});
			}

			void GetSampleNames(std::vector<std::string>& sample_names) const
			{
				sample_names.clear();
				const auto names_sample_stream_id = archive->get_stream_id(stream_names::SAMPLE_NAMES);
				//meaning it was not stored in the archive
				if (names_sample_stream_id == -1)
					return;

				std::vector<uint8_t> serialized;
				size_t meta;
				if (!archive->get_part(names_sample_stream_id, serialized, meta))
					throw std::runtime_error("Stream" + stream_names::SAMPLE_NAMES + " is empty");

				sample_names.reserve(metadata.config.num_samples);

				size_t pos{};
				for (size_t i = 0 ; i < metadata.config.num_samples ; ++i)
				{
					sample_names.emplace_back();
					refresh::serialization::load_string(sample_names.back(), serialized, pos);
				}

				if (archive->get_part(names_sample_stream_id, serialized, meta))
					throw std::runtime_error("Unexpected data in stream" + stream_names::SAMPLE_NAMES);
			}

			std::vector<std::string> GetSampleNames() const
			{
				std::vector<std::string> res;
				GetSampleNames(res);
				return res;
			}
		};

		template<typename VALUE_T>
		class ReaderSortedPlainBase : public ReaderBase<VALUE_T>
		{
		protected:
			using ReaderBase<VALUE_T>::metadata;
			using ReaderBase<VALUE_T>::archive;
		public:
			ReaderSortedPlainBase(MetadataReader& metadata_reader) :
				ReaderBase<VALUE_T>(metadata_reader)
			{
				if (metadata.kmers_representation != detail::KmersRepresentation::SortedPlain)
					throw std::runtime_error("This k-mer representation (" + detail::to_string(metadata.kmers_representation) + ") is not supported by ReaderSortedPlainBase");
			}
		};

		template<typename VALUE_T>
		class ReaderSortedWithLUTBase : public ReaderBase<VALUE_T>
		{
		protected:
			using ReaderBase<VALUE_T>::metadata;
			using ReaderBase<VALUE_T>::archive;
		public:
			ReaderSortedWithLUTBase(MetadataReader& metadata_reader) :
				ReaderBase<VALUE_T>(metadata_reader)
			{
				if (metadata.kmers_representation != detail::KmersRepresentation::SortedWithLUT)
					throw std::runtime_error("This k-mer representation (" + detail::to_string(metadata.kmers_representation) + ") is not supported by ReaderSortedWithLUTBase");
			}
		};
	}


	template<typename VALUE_T>
	class ReaderSortedPlainForListing : public detail::ReaderSortedPlainBase<VALUE_T>
	{
		using detail::ReaderBase<VALUE_T>::metadata;
		using detail::ReaderBase<VALUE_T>::archive;
		std::vector<std::unique_ptr<BinReaderSortedPlainForListing<VALUE_T>>> bins{}; //mkokot_TODO: maybe template base class on a bin type, and move "bins" member to base?
	public:
		ReaderSortedPlainForListing(MetadataReader& metadata_reader,
			size_t max_single_read_size = 1ull << 20) : //mkokot_TODO: is this a good default?
			detail::ReaderSortedPlainBase<VALUE_T>(metadata_reader)
		{
			//mkokot_TODO: in each reader for listing it may be small issue that I create all the bins
			//the problem is that it allocates memory for internal buffer
			//when we only want to iterate bins, but in sequential order (or even in parallel we probably will not read all at once)
			//so it should work such that it allocates memory at first read (at bin level) and deallocates after the last read
			//probably this will result in having all in memory only if we want final in sorted order
			//maybe the issue is only to deallocate after last read, i.e. when NextKmer returns false
			bins.reserve(metadata.config.num_bins);
			for (uint64_t i = 0; i < metadata.config.num_bins; ++i)
				bins.emplace_back(std::make_unique<BinReaderSortedPlainForListing<VALUE_T>>(
					i,
					archive,
					metadata.config.kmer_len,
					metadata.config.num_samples,
					metadata.config.num_bytes_single_value,
					max_single_read_size));
		}

		BinReaderSortedPlainForListing<VALUE_T>* GetBin(uint32_t id)
		{
			return bins[id].get();
		}
	};

	template<typename VALUE_T>
	class ReaderSortedPlainForRandomAccess : public detail::ReaderSortedPlainBase<VALUE_T>
	{
		using detail::ReaderBase<VALUE_T>::metadata;
		using detail::ReaderBase<VALUE_T>::archive;
		std::vector<std::unique_ptr<BinReaderSortedPlainForRandomAccess<VALUE_T>>> bins{};
	public:
		ReaderSortedPlainForRandomAccess(MetadataReader& metadata_reader) :
			detail::ReaderSortedPlainBase<VALUE_T>(metadata_reader)
		{
			bins.reserve(metadata.config.num_bins);
			for (uint64_t i = 0; i < metadata.config.num_bins; ++i)
				bins.emplace_back(std::make_unique<BinReaderSortedPlainForRandomAccess<VALUE_T>>(
					i,
					archive,
					metadata.config.kmer_len,
					metadata.config.num_samples,
					metadata.config.num_bytes_single_value));
		}

		// mkokot_TODO: this may be not needed, but there should be a method to search
		// for the whole database, but bin id determination is not ready yet
		BinReaderSortedPlainForRandomAccess<VALUE_T>* GetBin(uint32_t id)
		{
			return bins[id].get();
		}

		template<unsigned SIZE>
		bool CheckKmer(const CKmer<SIZE>& kmer, VALUE_T* values) const
		{
			//auto s = std::chrono::high_resolution_clock::now();
			auto bin_id = get_bin_id(kmer, metadata.config.kmer_len,
				metadata.config.signature_len, metadata.config.num_bins,
				metadata.config.signature_selection_scheme,
				metadata.config.signature_to_bin_mapping);

			return bins[bin_id]->CheckKmer(kmer, values);
		}
	};

	template<typename VALUE_T>
	class ReaderSortedWithLUTForListing : public detail::ReaderSortedWithLUTBase<VALUE_T>
	{
		using detail::ReaderBase<VALUE_T>::metadata;
		using detail::ReaderBase<VALUE_T>::archive;
		std::vector<std::unique_ptr<BinReaderSortedWithLUTForListing<VALUE_T>>> bins{};
	public:
		ReaderSortedWithLUTForListing(MetadataReader& metadata_reader,
			size_t max_single_read_size = 1ull << 20) : //mkokot_TODO: is this a good default?
			detail::ReaderSortedWithLUTBase<VALUE_T>(metadata_reader)
		{
			bins.reserve(metadata.config.num_bins);
			for (uint64_t i = 0; i < metadata.config.num_bins; ++i)
				bins.emplace_back(std::make_unique<BinReaderSortedWithLUTForListing<VALUE_T>>(
					i,
					archive,
					metadata.config.kmer_len,
					metadata.config.num_samples,
					std::get<ConfigSortedWithLUT>(metadata.representation_config).lut_prefix_len,
					metadata.config.num_bytes_single_value,
					max_single_read_size));
		}
		BinReaderSortedWithLUTForListing<VALUE_T>* GetBin(uint32_t id)
		{
			return bins[id].get();
		}
	};

	template<typename VALUE_T>
	class ReaderSortedWithLUTForRandomAccess : public detail::ReaderSortedWithLUTBase<VALUE_T>
	{
		using detail::ReaderBase<VALUE_T>::metadata;
		using detail::ReaderBase<VALUE_T>::archive;
		std::vector<std::unique_ptr<BinReaderSortedWithLUTForRandomAccess<VALUE_T>>> bins{};
	public:
		ReaderSortedWithLUTForRandomAccess(MetadataReader& metadata_reader) :
			detail::ReaderSortedWithLUTBase<VALUE_T>(metadata_reader)
		{
			bins.reserve(metadata.config.num_bins);
			for (uint64_t i = 0; i < metadata.config.num_bins; ++i)
				bins.emplace_back(std::make_unique<BinReaderSortedWithLUTForRandomAccess<VALUE_T>>(
					i,
					archive,
					metadata.config.kmer_len,
					metadata.config.num_samples,
					std::get<ConfigSortedWithLUT>(metadata.representation_config).lut_prefix_len,
					metadata.config.num_bytes_single_value));
		}

		BinReaderSortedWithLUTForRandomAccess<VALUE_T>* GetBin(uint32_t id)
		{
			return bins[id].get();
		}

		template<unsigned SIZE>
		bool CheckKmer(const CKmer<SIZE>& kmer, VALUE_T* values) const
		{
			auto bin_id = get_bin_id(kmer, metadata.config.kmer_len,
				metadata.config.signature_len, metadata.config.num_bins,
				metadata.config.signature_selection_scheme,
				metadata.config.signature_to_bin_mapping);

			return bins[bin_id]->CheckKmer(kmer, values);
		}
	};
}

#endif// ! READERS_H_
