#ifndef METADATA_READER_H_
#define METADATA_READER_H_

#include "stream_names.h"
#include "utils.h"
#include "metadata.h"

namespace kmcdb
{
	class MetadataReader;
	namespace detail
	{
		//some trick because I don't want to expose archive as a public method for a client code
		//at least not explicitly, I will assume client code will not use the detail namespace
		inline archive_input_t* get_archive_from_metadata_reader(MetadataReader& metadata_reader);
		inline Metadata& get_metadata_from_metadata_reader(MetadataReader& metadata_reader);
	}

	//this class will just read the metadata which may be needed to decide which class to use further, i.e. which template concretization
	class MetadataReader
	{
		friend archive_input_t* detail::get_archive_from_metadata_reader(MetadataReader& metadata_reader);
		friend detail::Metadata& detail::get_metadata_from_metadata_reader(MetadataReader& metadata_reader);
		archive_input_t archive;
		detail::Metadata metadata;
		int stream_id;

		std::string get_captured_impl(const std::string& stream_name)
		{
			const auto captured_stream_id = archive.get_stream_id(stream_name);
			//if was not captured just return empty string
			if (captured_stream_id == -1)
				return  "";

			std::string result;
			result.reserve(archive.get_stream_data_size(captured_stream_id));
			const auto no_parts = archive.get_no_parts(captured_stream_id);
			for (size_t part_id = 0; part_id < no_parts; ++part_id)
			{
				//mkokot_TODO: for now I am reading to vec because there is no method to get via raw pointer while specifying part_id
				//in the future either use rewind or raw pointer version
				std::vector<uint8_t> data;
				uint64_t meta;
				auto res = archive.get_part(captured_stream_id, static_cast<int>(part_id), data, meta);
				if (!res)
					throw std::runtime_error("Cannot read expected part form stream id " + std::to_string(captured_stream_id));

				result.insert(result.end(), data.begin(), data.end());
			}
			return result;
		}
	public:
		MetadataReader(const std::string& path, const bool reopen_mode)
		{
			//mkokot_TODO: different exception type?
			//if (!archive.open_file_unbuffered(path, reopen_mode))
			if (!archive.open_file_memory_mapped(path, reopen_mode))
				throw std::runtime_error("Cannot open file " + path);

			stream_id = archive.get_stream_id(stream_names::METADATA);

			if (stream_id == -1)
				throw std::runtime_error("Cannot find " + stream_names::METADATA + " stream in file " + path);

			std::vector<uint8_t> vec;
			uint64_t tmp;
			if (!archive.get_part(stream_id, vec, tmp))
				throw std::runtime_error("Stream" + stream_names::METADATA + " is empty in file" + path);

			metadata.load(vec);

			if (archive.get_part(stream_id, vec, tmp))
				throw std::runtime_error("Unexpected data in stream" + stream_names::METADATA + " in file" + path);

			if (!detail::IsCompatible(metadata.config.kmcdb_version))
				throw std::runtime_error("Incompatible kmcdb API (" + KMCDB_VERSION.to_string() +
					") and kmcdb file (" + metadata.config.kmcdb_version.to_string() + ") versions.");
		}

		const Config& GetConfig() const
		{
			return metadata.config;
		}
	};
}
#endif // ! METADATA_READER_H_

