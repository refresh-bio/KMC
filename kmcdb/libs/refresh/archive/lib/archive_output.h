#pragma once

#include <atomic>

#include "archive_common.h"
#include "output_stdout_memory.h"
#include "output_file_buffered.h"
#include "output_file_unbuffered.h"
#include "output_file_memory_mapped.h"
#include "output_file_low_level.h"

namespace refresh
{
	class archive_output : public archive_common
	{
	private:
		std::unique_ptr<refresh::io::output_common> output;
		size_t f_offset{};

		// *******************************************************************************************
		size_t serialize_params()
		{
			size_t params_size = 0;
			size_t nb = 0;

			std::map<std::string, uint64_t> dict;

			dict.emplace(params.str_archive_version, (uint64_t)params.archive_version);
			dict.emplace(params.str_parts_metadata_fixed_size, (uint64_t)params.parts_metadata_fixed_size);
			dict.emplace(params.str_parts_metadata_empty, (uint64_t)params.parts_metadata_empty);

			params_size += (nb = write(dict.size()));
			if (nb == 0)	return 0;

			for (const auto& par : dict)
			{
				params_size += (nb = write(par.first));
				if (nb == 0)	return 0;

				params_size += (nb = write(par.second));
				if (nb == 0)	return 0;
			}

			return params_size;
		}

		// *******************************************************************************************
		bool serialize()
		{
			size_t footer_size = 0;
			size_t nb = 0;

			// Store stream part offsets
			footer_size += (nb = write(v_streams.size()));
			if (nb == 0)	return 0;

			for (auto& stream : v_streams)
			{
				size_t p = footer_size;

				footer_size += (nb = write(stream.stream_name));
				if (nb == 0)	return 0;
				footer_size += (nb = write(stream.parts.size()));
				if (nb == 0)	return 0;
				footer_size += (nb = write(stream.metadata));
				if (nb == 0)	return 0;

				for (auto& part : stream.parts)
				{
					footer_size += (nb = write(part.offset));
					if (nb == 0)	return 0;
					footer_size += (nb = write(part.size));
					if (nb == 0)	return 0;
				}

				stream.total_size += footer_size - p;
			}

			// Optional params dictionary
			footer_size += (nb = serialize_params());
			if (nb == 0)	return 0;

			nb = write_fixed(footer_size);
			if (nb == 0)	return 0;

			return true;
		}

		// *******************************************************************************************
		template<typename T>
		size_t write_fixed(const T x)
		{
			if (output->write_uint(static_cast<uint64_t>(x), 8) != 8)
				return 0;

			return 8;
		}

		// *******************************************************************************************
		template<typename T>
		size_t write(const T _x)
		{
			uint8_t no_bytes = 0;
			uint64_t x = static_cast<uint64_t>(_x);

			for (size_t tmp = x; tmp; tmp >>= 8)
				++no_bytes;

			if (!output->put(no_bytes))
				return 0;

			for (int i = no_bytes; i; --i)
				if (!output->put((x >> ((i - 1) * 8)) & 0xff))
					return 0;

			return no_bytes + 1;
		}

		// *******************************************************************************************
		size_t write(const std::string& s)
		{
			if (output->write((uint8_t*)s.c_str(), s.size()) != s.size())
				return 0;
			if (!output->put(0))
				return 0;

			return s.size() + 1;
		}

		// *******************************************************************************************
		int register_stream_impl(const std::string& stream_name)
		{
			// Before adding new stream check if stream_name is already registered
			auto p = rm_streams.find(stream_name);
			if (p != rm_streams.end())
				return (int)p->second;

			int id = (int)v_streams.size();

			v_streams.emplace_back(stream_t());

			v_streams[id].cur_id = 0;
			v_streams[id].stream_name = stream_name;
			v_streams[id].metadata = 0;
			v_streams[id].total_size = 0;
			v_streams[id].data_size = 0;

			rm_streams[stream_name] = id;

			return id;
		}

		// *******************************************************************************************
		bool add_part_impl(const int stream_id, const uint8_t* data, const size_t data_size, const uint64_t metadata)
		{
			size_t nb = 0;

			v_streams[stream_id].parts.push_back(part_t(f_offset, data_size));

			if (params.parts_metadata_empty)
				;
			else if (params.parts_metadata_fixed_size)
			{
				f_offset += (nb = write_fixed(metadata));
				if (nb == 0)	return false;
			}
			else
			{
				f_offset += (nb = write(metadata));
				if (nb == 0)	return false;
			}

			if (output->write(data, data_size) != data_size)
				return false;

			f_offset += data_size;

			v_streams[stream_id].total_size += f_offset - v_streams[stream_id].parts.back().offset;
			v_streams[stream_id].data_size += data_size;

			return true;
		}

		// *******************************************************************************************
		int add_part_prepare_impl(const int stream_id)
		{
			v_streams[stream_id].parts.push_back(part_t(0, 0));

			return static_cast<int>(v_streams[stream_id].parts.size()) - 1;
		}

		// *******************************************************************************************
		bool add_part_complete_impl(const int stream_id, const int part_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
		{
			size_t nb = 0;

			v_streams[stream_id].parts[part_id] = part_t(f_offset, v_data.size());

			if (params.parts_metadata_empty)
				;
			else if (params.parts_metadata_fixed_size)
			{
				f_offset += (nb = write_fixed(metadata));
				if (nb == 0)	return false;
			}
			else
			{
				f_offset += (nb = write(metadata));
				if (nb == 0)	return false;
			}

			if(output->write(v_data.data(), v_data.size()) != v_data.size())
				return false;

			f_offset += v_data.size();

			v_streams[stream_id].total_size += f_offset - v_streams[stream_id].parts[part_id].offset;
			v_streams[stream_id].data_size += v_data.size();

			return true;
		}

		// *******************************************************************************************
		bool flush_out_buffers_imp()
		{
			for (auto& x : m_buffer)
				for (auto& y : x.second)
					if (!add_part_impl(x.first, y.first.data(), y.first.size(), y.second))
						return false;

			m_buffer.clear();

			return true;
		}

		// *******************************************************************************************
		bool open_check()
		{
			if (!output->opened())
			{
				err_code = ec_file_open;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		bool archive_params_check(const params_t& _params)
		{
			if (!_params.verify())
			{
				err_code = ec_archive_params;
				return false;
			}

			return true;
		}

	public:
		// *******************************************************************************************
		archive_output() : archive_common()
		{}

		// *******************************************************************************************
		bool open_memory(std::vector<uint8_t>& data, const params_t& _params = params_t())
		{
			if (!archive_params_check(_params))
				return false;

			output = std::make_unique<refresh::io::output_memory>(data);
			params = _params;

			return open_check();
		}

		// *******************************************************************************************
		bool open_stdout_unbuffered(const params_t& _params = params_t())
		{
			if (!archive_params_check(_params))
				return false;

			output = std::make_unique<refresh::io::output_stdout>();
			params = _params;

			return open_check();
		}

		// *******************************************************************************************
		bool open_stdout_buffered(const size_t buffer_size = 1 << 20, const params_t& _params = params_t())
		{
			if (!archive_params_check(_params))
				return false;

			output = std::make_unique<refresh::io::output_stdout_buffered>(buffer_size);
			params = _params;

			return open_check();
		}

		// *******************************************************************************************
		bool open_file_unbuffered(const std::string& file_name, const bool reopen_mode, const params_t& _params = params_t())
		{
			if (!archive_params_check(_params))
				return false;

			if(reopen_mode)
				output = std::make_unique<refresh::io::output_file_unbuffered_reopen>(file_name);
			else
				output = std::make_unique<refresh::io::output_file_unbuffered>(file_name);
			params = _params;

			return open_check();
		}

		// *******************************************************************************************
		bool open_file_buffered(const std::string& file_name, const bool reopen_mode, const size_t buffer_size = 32 << 20, const params_t& _params = params_t())
		{
			if (!archive_params_check(_params))
				return false;

			if (reopen_mode)
				output = std::make_unique<refresh::io::output_file_buffered_reopen>(file_name, buffer_size);
			else
				output = std::make_unique<refresh::io::output_file_buffered>(file_name, buffer_size);
			params = _params;

			return open_check();
		}

		// *******************************************************************************************
		bool open_file_low_level(const std::string& file_name, const size_t buffer_size = 32 << 20, const params_t& _params = params_t())
		{
			if (!archive_params_check(_params))
				return false;

			output = std::make_unique<refresh::io::output_file_low_level>(file_name, buffer_size);
			params = _params;

			return open_check();
		}

		// *******************************************************************************************
		bool close()
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (output)
			{
				if (!flush_out_buffers_imp())
				{
					err_code = ec_file_write;
					return false;
				}
				if (!output->start_transaction())
				{
					err_code = ec_file_write;
					return false;
				}
				if (!serialize())
				{
					err_code = ec_file_write;
					return false;
				}
				if (!output->stop_transaction())
				{
					err_code = ec_file_write;
					return false;
				}
				if (!output->close())
				{
					err_code = ec_file_close;
					return false;
				}

				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_open;
				return false;
			}
		}

		// *******************************************************************************************
		int register_stream(const std::string& stream_name)
		{
			std::lock_guard<std::mutex> lck(mtx);

			return register_stream_impl(stream_name);
		}

		// *******************************************************************************************
		std::vector<int> register_streams(const std::vector<std::string>& stream_names)
		{
			std::lock_guard<std::mutex> lck(mtx);

			std::vector<int> ret;

			ret.reserve(stream_names.size());

			for (const auto& sn : stream_names)
				ret.emplace_back(register_stream_impl(sn));

			return ret;
		}

		// *******************************************************************************************
		bool add_part(const int stream_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!output->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			if (!add_part_impl(stream_id, v_data.data(), v_data.size(), metadata))
			{
				err_code = ec_file_write;
				return false;
			}

			if (!output->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		bool add_part(const int stream_id, const uint8_t* data, const size_t data_size, const uint64_t metadata = 0)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!output->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			if (!add_part_impl(stream_id, data, data_size, metadata))
			{
				err_code = ec_file_write;
				return false;
			}

			if (!output->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		// metadatas can be empty!
		bool add_parts(const std::vector<int>& stream_ids, const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, const std::vector<uint64_t>& metadatas)
		{
			std::lock_guard<std::mutex> lck(mtx);

			auto n_items = stream_ids.size();

			if (v_datas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (!metadatas.empty() && metadatas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (!output->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			for (size_t i = 0; i < n_items; ++i)
				if (!add_part_impl(stream_ids[i], v_datas.at(i).get().data(), v_datas.at(i).get().size(), metadatas.empty() ? 0 : metadatas[i]))
				{
					err_code = ec_file_write;
					return false;
				}

			if (!output->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		// metadatas can be empty!
		bool add_parts(const std::vector<int>& stream_ids, const std::vector<uint8_t*>& datas, const std::vector<size_t> data_sizes, const std::vector<uint64_t>& metadatas)
		{
			std::lock_guard<std::mutex> lck(mtx);

			auto n_items = stream_ids.size();

			if (datas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (data_sizes.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (!metadatas.empty() && metadatas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (!output->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			for (size_t i = 0; i < n_items; ++i)
				if (!add_part_impl(stream_ids[i], datas[i], data_sizes[i], metadatas.empty() ? 0 : metadatas[i]))
				{
					err_code = ec_file_write;
					return false;
				}

			if (!output->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		int add_part_prepare(const int stream_id)
		{
			std::lock_guard<std::mutex> lck(mtx);

			return add_part_prepare_impl(stream_id);
		}

		// *******************************************************************************************
		std::vector<int> add_parts_prepare(const std::vector<int>& stream_ids)
		{
			std::lock_guard<std::mutex> lck(mtx);

			std::vector<int> ret;

			ret.reserve(stream_ids.size());

			for (const auto sid : stream_ids)
				ret.emplace_back(add_part_prepare_impl(sid));

			return ret;
		}

		// *******************************************************************************************
		bool add_part_complete(const int stream_id, const int part_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!output->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			if (!add_part_complete_impl(stream_id, part_id, v_data, metadata))
			{
				err_code = ec_file_write;
				return false;
			}

			if (!output->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		// metadatas can be empty
		bool add_parts_complete(const std::vector<int>& stream_ids, const std::vector<int>& part_ids, const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, const std::vector<uint64_t>& metadatas)
		{
			std::lock_guard<std::mutex> lck(mtx);

			auto n_items = stream_ids.size();

			if (part_ids.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (v_datas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (!metadatas.empty() && metadatas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (!output->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			for (size_t i = 0; i < n_items; ++i)
				if (!add_part_complete_impl(stream_ids[i], part_ids[i], v_datas.at(i), metadatas.empty() ? 0 : metadatas[i]))
				{
					err_code = ec_file_write;
					return false;
				}

			if (!output->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		bool add_part_buffered(const int stream_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
		{
			std::lock_guard<std::mutex> lck(mtx);

			m_buffer[stream_id].emplace_back(v_data, metadata);

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		bool add_parts_buffered(const std::vector<int>& stream_ids, const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, const std::vector<uint64_t>& metadatas)
		{
			std::lock_guard<std::mutex> lck(mtx);

			auto n_items = stream_ids.size();

			if (v_datas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			if (!metadatas.empty() && metadatas.size() != n_items)
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			for (size_t i = 0; i < n_items; ++i)
				m_buffer[stream_ids[i]].emplace_back(v_datas.at(i), metadatas.empty() ? 0 : metadatas[i]);

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		bool flush_out_buffers()
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!output->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			if (!flush_out_buffers_imp())
			{
				err_code = ec_file_write;
				return false;
			}

			if (!output->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

		// *******************************************************************************************
		void set_stream_metadata(const int stream_id, const uint64_t metadata)
		{
			std::lock_guard<std::mutex> lck(mtx);

			v_streams[stream_id].metadata = metadata;
		}
	};
}