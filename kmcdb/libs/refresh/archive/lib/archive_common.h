#pragma once

#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <thread>
#include <mutex>
#include <functional>
#include <cinttypes>

namespace refresh
{
	const uint32_t REFRESH_BUILD_ARCHIVE = 2;

	class archive_output; //forward declaration for friendness
	class archive_common
	{
	public:
		struct part_t {
			size_t offset;
			size_t size;

			part_t() : offset(0), size(0)
			{};

			part_t(size_t _offset, size_t _size) : offset(_offset), size(_size)
			{};
		};

		struct stream_t {
			std::string stream_name;
			size_t cur_id;
			uint64_t metadata;
			size_t total_size;
			size_t data_size;
			std::vector<part_t> parts;
			size_t sub_part_id;
			size_t sub_part_offset;
			size_t sub_part_metadata_size;
			uint64_t sub_part_metadata;

			stream_t() :
				stream_name(""),
				cur_id(0),
				metadata(0),
				total_size(0),
				data_size(0),
				sub_part_id(0),
				sub_part_offset(0),
				sub_part_metadata_size(0),
				sub_part_metadata(0)
			{}

			stream_t(const std::string& stream_name, size_t cur_id, uint64_t metadata, size_t total_size, size_t data_size) :
				stream_name(stream_name),
				cur_id(cur_id),
				metadata(metadata),
				total_size(total_size),
				data_size(data_size),
				sub_part_id(0),
				sub_part_offset(0),
				sub_part_metadata_size(0),
				sub_part_metadata(0)
			{}

			stream_t(const stream_t&) = default;
			stream_t(stream_t&&) = default;

			stream_t& operator=(const stream_t&) = default;
			stream_t& operator=(stream_t&&) = default;
		};

		struct params_t
		{
			friend class archive_output;

			uint64_t archive_version = 1;
			bool parts_metadata_fixed_size = false;
			bool parts_metadata_empty = false;

			static const inline std::string str_archive_version = "archive_version";
			static const inline std::string str_parts_metadata_fixed_size = "parts_metadata_fixed_size";
			static const inline std::string str_parts_metadata_empty = "parts_metadata_empty";

		protected:
			// Internal method (only for archive* classes) to check if the parameters are consistent
			bool verify() const
			{
				if (archive_version == 1)
				{
					return (parts_metadata_empty == false) && (parts_metadata_fixed_size == false);
				}
				else if (archive_version == 2)
				{
					return true;
				}
				else return false;
			}
		} params;

	protected:
		std::map<int, std::vector<std::pair<std::vector<uint8_t>, uint64_t>>> m_buffer;
		std::vector<stream_t> v_streams;
		std::unordered_map<std::string, size_t> rm_streams;

		std::mutex mtx;

		static constexpr uint64_t ec_ok = 0;
		static constexpr uint64_t ec_file_open = 101;
		static constexpr uint64_t ec_file_close = 102;
		static constexpr uint64_t ec_file_write = 103;
		static constexpr uint64_t ec_file_read = 104;
		static constexpr uint64_t ec_archive_params = 105;
		static constexpr uint64_t ec_wrong_part_desc = 201;
		static constexpr uint64_t ec_unknown_stream = 301;
		static constexpr uint64_t ec_unknown_stream_part = 302;

		std::atomic<uint64_t> err_code{ ec_ok };

		// *******************************************************************************************
		int get_stream_id_impl(const std::string& stream_name) const
		{
			auto p = rm_streams.find(stream_name);
			if (p != rm_streams.end())
				return (int)p->second;

			return -1;
		}

	public:
		archive_common() = default;

		// *******************************************************************************************
		int get_stream_id(const std::string& stream_name) 
		{
			std::lock_guard<std::mutex> lck(mtx);

			return get_stream_id_impl(stream_name);
		}

		// *******************************************************************************************
		std::vector<int> get_stream_ids(const std::vector<std::string>& stream_names) 
		{
			std::lock_guard<std::mutex> lck(mtx);

			std::vector<int> ret;

			ret.reserve(stream_names.size());

			for (const auto& sn : stream_names)
				ret.emplace_back(get_stream_id_impl(sn));

			return ret;
		}

		// *******************************************************************************************
		size_t get_no_streams() 
		{
			std::lock_guard<std::mutex> lck(mtx);

			return v_streams.size();
		}

		// *******************************************************************************************
		uint64_t err_status()	const
		{
			return err_code.load();
		}

		// *******************************************************************************************
		std::string err_message(const uint64_t code) const
		{
			switch (code)
			{
			case ec_ok:	return "";
			case ec_archive_params: return "archive params inconsistent";
			case ec_file_open: return "cannot open file";
			case ec_file_close: return "cannot close file";
			case ec_file_write: return "cannot write";
			case ec_file_read: return "cannot read";
			case ec_wrong_part_desc: return "wrong size of parameters";
			case ec_unknown_stream: return "unknown stream id";
			case ec_unknown_stream_part: return "uknown stream part id";
			}

			return "unknown error";
		}
	};
}