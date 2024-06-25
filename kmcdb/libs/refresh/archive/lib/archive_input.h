#pragma once

#include <memory>
#include <filesystem>

#include "archive_common.h"
#include "input_stdin_memory.h"
#include "input_file_buffered.h"
#include "input_file_unbuffered.h"

namespace refresh
{
	class archive_input : public archive_common
	{
		std::unique_ptr<refresh::io::input_common> input;

		// *******************************************************************************************
		void unget_part(const int id)
		{
			if (v_streams[id].cur_id)
				--v_streams[id].cur_id;
		}

		// *******************************************************************************************
		bool deserialize_params()
		{
			std::map<std::string, uint64_t> dict;

			size_t dict_size;
			read(dict_size);
			std::string par_name;
			uint64_t par_val;

			// Load whole dictionary
			for (size_t i = 0; i < dict_size; ++i)
			{
				if (!read(par_name))		return false;
				if (!read(par_val))			return false;

				dict[par_name] = par_val;
			}

			// Set known parameters
			auto p = dict.find(params_t::str_archive_version);
			if (p != dict.end())
				params.archive_version = p->second;

			p = dict.find(params_t::str_parts_metadata_fixed_size);
			if (p != dict.end())
				params.parts_metadata_fixed_size = (bool)p->second;

			p = dict.find(params_t::str_parts_metadata_empty);
			if (p != dict.end())
				params.parts_metadata_empty = (bool)p->second;

			return true;
		}

		// *******************************************************************************************
		bool deserialize()
		{
			size_t nb = 0;
			size_t footer_size;
			size_t file_size_ = input->file_size();

			input->seek(file_size_ - 8ull);
			read_fixed(footer_size);

			input->seek(file_size_ - (size_t)(8 + footer_size));

			size_t readed = 0;

			// Read stream part offsets
			size_t n_streams;
			readed += (nb = read(n_streams));
			if (nb == 0)		return false;

			v_streams.resize(n_streams, stream_t());

			for (size_t i = 0; i < n_streams; ++i)
			{
				auto& stream_second = v_streams[i];

				readed += (nb = read(stream_second.stream_name));
				if (nb == 0)		return false;
				readed += (nb = read(stream_second.cur_id));
				if (nb == 0)		return false;
				readed += (nb = read(stream_second.metadata));
				if (nb == 0)		return false;

				stream_second.parts.resize(stream_second.cur_id);
				for (size_t j = 0; j < stream_second.cur_id; ++j)
				{
					readed += (nb = read(stream_second.parts[j].offset));
					if (nb == 0)		return false;
					readed += (nb = read(stream_second.parts[j].size));
					if (nb == 0)		return false;

					stream_second.total_size += stream_second.parts[j].size + 8;		// approximated (cannot check part metadata sizes w/o reading the archive)
					stream_second.data_size += stream_second.parts[j].size;
				}

				stream_second.cur_id = 0;

				rm_streams[stream_second.stream_name] = i;
			}

			// Optional params dictionary 
			if (readed < footer_size)
				if (!deserialize_params())
					return false;

			if (!input->seek(0))
				return false;

			return true;
		}

		// *******************************************************************************************
		template<typename T>
		size_t read_fixed(T& x)
		{
			uint64_t tmp;

			if (!input->read_uint(sizeof(T), tmp))
				return 0;

			x = static_cast<T>(tmp);

			return sizeof(T);
		}

		// *******************************************************************************************
		size_t read(std::string& s)
		{
			s.clear();

			while (true)
			{
				uint8_t c;
				if (!input->get(c))
					return 0;

				if (c == 0)
					return s.size() + 1;

				s.push_back((char)c);
			}

			return 0;			// Never be here
		}

		// *******************************************************************************************
		template<typename T>
		size_t read(T& x)
		{
			uint8_t no_bytes;
			uint8_t c;

			if (!input->get(no_bytes))
				return 0;

			x = 0;

			for (int i = 0; i < no_bytes; ++i)
			{
				x <<= 8;
				if (!input->get(c))
					return 0;
				x += static_cast<T>(c);
			}

			return no_bytes + 1;
		}

		// *******************************************************************************************
		size_t get_stream_total_size_impl(const int stream_id)
		{
			if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
				return 0;

			return v_streams[stream_id].total_size;
		}

		// *******************************************************************************************
		size_t get_stream_data_size_impl(const int stream_id)
		{
			if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
				return 0;

			return v_streams[stream_id].data_size;
		}

		// *******************************************************************************************
		bool get_part_impl(const int stream_id, uint8_t* data, uint64_t& metadata)
		{
			auto& p = v_streams[stream_id];

			input->seek(p.parts[p.cur_id].offset);

			if (p.parts[p.cur_id].size != 0)
			{
				if (params.parts_metadata_empty)
					metadata = 0;
				else if (params.parts_metadata_fixed_size)
				{
					if (read_fixed(metadata) == 0)
						return false;
				}
				else
				{
					if (read(metadata) == 0)
						return false;
				}
			}
			else
			{
				metadata = 0;
				p.cur_id++;
				return true;
			}

			if (input->read(data, p.parts[p.cur_id].size) != p.parts[p.cur_id].size)
				return false;

			p.cur_id++;

			p.sub_part_id = 0;
			p.sub_part_offset = 0;
			p.sub_part_metadata = 0;
			p.sub_part_metadata_size = 0;

			return true;
		}

		// *******************************************************************************************
		bool get_part_impl(const int stream_id, const int part_id, uint8_t* data, uint64_t& metadata)
		{
			auto& p = v_streams[stream_id];

			if ((size_t)part_id >= p.parts.size())
				return false;

			if (!input->seek(p.parts[part_id].offset))
				return false;

			if (p.parts[part_id].size != 0)
			{
				if (params.parts_metadata_empty)
					metadata = 0;
				else if (params.parts_metadata_fixed_size)
				{
					if (read_fixed(metadata) == 0)
						return false;
				}
				else
				{
					if (read(metadata) == 0)
						return false;
				}
			}
			else
			{
				metadata = 0;
				return true;
			}

			if (input->read(data, p.parts[part_id].size) != p.parts[part_id].size)
				return false;

			return true;
		}

		// *******************************************************************************************
		bool get_part_size_impl(const int stream_id, size_t& size)
		{
			auto& p = v_streams[stream_id];

			if (p.cur_id >= p.parts.size())
				return false;

			size = p.parts[p.cur_id].size;

			return true;
		}

		// *******************************************************************************************
		bool get_part_size_impl(const int stream_id, const int part_id, size_t& size)
		{
			if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
				return false;

			auto& p = v_streams[stream_id];

			if (part_id >= static_cast<int>(p.parts.size()))
				return false;

			size = p.parts[part_id].size;

			return true;
		}

		// *******************************************************************************************
		bool get_parts_impl(const std::vector<int>& stream_ids, std::vector<uint8_t*> datas, std::vector<uint64_t>& metadatas)
		{
			size_t no_streams = stream_ids.size();

			for (size_t i = 0; i < no_streams; ++i)
				if (!get_part_impl(stream_ids[i], datas[i], metadatas.at(i)))
				{
					for (int j = 0; j < (int) i; ++j)
						unget_part(j);

					return false;
				}

			return true;
		}

		// *******************************************************************************************
		bool get_sub_part_impl(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, bool& last_sub_part, uint8_t* data, size_t& sub_part_size, uint64_t& metadata)
		{
			auto& p = v_streams[stream_id];

			if (p.cur_id >= p.parts.size())
				return false;

			part_id = static_cast<int>(p.cur_id);
			sub_part_id = static_cast<int>(p.sub_part_id++);

			if (sub_part_id == 0)
			{
				input->seek(p.parts[p.cur_id].offset);
				p.sub_part_offset = 0;

				if (p.parts[p.cur_id].size != 0)
				{
					if (params.parts_metadata_empty)
						p.sub_part_metadata_size = 0;
					else if (params.parts_metadata_fixed_size)
					{
						p.sub_part_metadata_size = read_fixed(p.sub_part_metadata);
						if (p.sub_part_metadata_size == 0)
							return false;
					}
					else
					{
						p.sub_part_metadata_size = read(p.sub_part_metadata);
						if (p.sub_part_metadata_size == 0)
							return false;
					}
				}
				else
				{
					metadata = 0;
					sub_part_size = 0;
					p.cur_id++;
					last_sub_part = true;
					return true;
				}
			}
			else
			{
				if (!input->seek(p.parts[p.cur_id].offset + p.sub_part_metadata_size + p.sub_part_offset))
					return false;
			}

			sub_part_size = std::min(max_sub_part_size, p.parts[p.cur_id].size - p.sub_part_offset);

			input->read(data, sub_part_size);
			p.sub_part_offset += sub_part_size;

			metadata = p.sub_part_metadata;

			if (p.parts[p.cur_id].size == p.sub_part_offset)
			{
				p.cur_id++;
				p.sub_part_offset = 0;
				p.sub_part_id = 0;
				p.sub_part_metadata_size = 0;
				last_sub_part = true;
			}
			else
				last_sub_part = false;

			return true;
		}

		// *******************************************************************************************
		bool get_excerpt_impl(const int stream_id, const size_t max_excerpt_size, const size_t offset, uint8_t* v_data, size_t& excerpt_size)
		{
			size_t acc_part_size = 0;

			excerpt_size = 0;

			if (stream_id < 0 || stream_id >= v_streams.size())
				return false;

			for (const auto& part : v_streams[stream_id].parts)
			{
				if (offset < acc_part_size + part.size)
				{
					size_t metadata_len;

					if (!params.parts_metadata_fixed_size)
					{
						input->seek(part.offset);
						uint8_t c;
						input->get(c);
						metadata_len = c + 1;
					}
					else
						metadata_len = 8;

					size_t local_offset = offset - acc_part_size;
					excerpt_size = std::min<size_t>(max_excerpt_size, part.size - local_offset);

					if (!input->seek(part.offset + metadata_len + local_offset))
						return false;
					input->read(v_data, excerpt_size);

					return excerpt_size == max_excerpt_size;
				}

				acc_part_size += part.size;
			}

			return false;
		}

		// *******************************************************************************************
		bool rewind_impl(const int stream_id)
		{
			if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
				return false;

			auto& s = v_streams[stream_id];

			s.cur_id = 0;
			s.sub_part_id = 0;
			s.sub_part_offset = 0;
			s.sub_part_metadata_size = 0;
			s.sub_part_metadata = 0;

			return true;
		}

		// *******************************************************************************************
		bool prepare()
		{
			if (!input->opened())
			{
				err_code = ec_file_open;
				return false;
			}

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			if (!deserialize())
			{
				err_code = ec_file_read;
				return false;
			}

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			err_code = ec_ok;
			return true;
		}

	public:
		// *******************************************************************************************
		archive_input() : archive_common()
		{}

		// *******************************************************************************************
		bool open_memory(const std::vector<uint8_t>& data)
		{
			input = std::make_unique<refresh::io::input_memory>(data);

			return prepare();
		}

		// *******************************************************************************************
		bool open_memory(std::vector<uint8_t>&& data)
		{
			input = std::make_unique<refresh::io::input_memory>(data);

			return prepare();
		}

		// *******************************************************************************************
		bool open_stdin()
		{
			input = std::make_unique<refresh::io::input_stdin>();

			return prepare();
		}
		
		// *******************************************************************************************
		bool open_file_unbuffered(const std::string& file_name, const bool reopen_mode)
		{
			if (reopen_mode)
				input = std::make_unique<refresh::io::input_file_unbuffered_reopen>(file_name);
			else
				input = std::make_unique<refresh::io::input_file_unbuffered>(file_name);

			return prepare();
		}

		// *******************************************************************************************
		bool open_file_all_buffered(const std::string& file_name)
		{
			input = std::make_unique<refresh::io::input_file_all_buffered>(file_name);

			return prepare();
		}

		// *******************************************************************************************
		bool open_file_buffered(const std::string& file_name, const bool reopen_mode, const size_t buffer_size = 32 << 20)
		{
			std::filesystem::path fp(file_name);

			if (!std::filesystem::exists(fp))
			{
				err_code = ec_file_open;
				return false;
			}
			
			if (std::filesystem::file_size(fp) <= buffer_size)
				return open_file_all_buffered(file_name);

			if (reopen_mode)
				input = std::make_unique<refresh::io::input_file_buffered_reopen>(file_name, buffer_size);
			else
				input = std::make_unique<refresh::io::input_file_buffered>(file_name, buffer_size);

			return prepare();
		}

		// *******************************************************************************************
		bool close()
		{
			if (input)
			{
				bool r = input->close();
				input.reset();

				return r;
			}

			return true;
		}

		// *******************************************************************************************
		size_t get_stream_total_size(const int stream_id)
		{
			std::lock_guard<std::mutex> lck(mtx);

			return get_stream_total_size_impl(stream_id);
		}

		// *******************************************************************************************
		std::vector<size_t> get_stream_total_sizes(const std::vector<int>& stream_ids)
		{
			std::lock_guard<std::mutex> lck(mtx);

			std::vector<size_t> ret;

			ret.reserve(stream_ids.size());

			for (const int sid : stream_ids)
				ret.emplace_back(get_stream_total_size(sid));

			return ret;
		}

		// *******************************************************************************************
		size_t get_stream_data_size(const int stream_id)
		{
			std::lock_guard<std::mutex> lck(mtx);

			return get_stream_data_size_impl(stream_id);
		}

		// *******************************************************************************************
		std::vector<size_t> get_stream_data_sizes(const std::vector<int>& stream_ids)
		{
			std::lock_guard<std::mutex> lck(mtx);

			std::vector<size_t> ret;

			ret.reserve(stream_ids.size());

			for (const int sid : stream_ids)
				ret.emplace_back(get_stream_data_size_impl(sid));

			return ret;
		}

		// *******************************************************************************************
		bool get_part_size(const int stream_id, size_t& size)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (get_part_size_impl(stream_id, size))
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_unknown_stream;
				return false;
			}
		}
		// *******************************************************************************************
		bool get_part_size(const int stream_id, const int part_id, size_t& size)
		{
			std::lock_guard<std::mutex> lck(mtx);

			err_code = ec_ok;

			if(get_part_size_impl(stream_id, part_id, size))
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_unknown_stream_part;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_part_sizes(const std::vector<int>& stream_ids, std::vector<size_t>& sizes)
		{
			std::lock_guard<std::mutex> lck(mtx);

			auto n_items = stream_ids.size();

			sizes.resize(n_items);

			bool ret = true;

			for (size_t i = 0; i < n_items; ++i)
				if (!get_part_size_impl(stream_ids[i], sizes[i]))
				{
					err_code = ec_unknown_stream_part;
					return false;
				}

			err_code = ec_ok;
			return ret;
		}

		// *******************************************************************************************
		bool get_part(const int stream_id, std::vector<uint8_t>& v_data, uint64_t& metadata)
		{
			std::lock_guard<std::mutex> lck(mtx);

			size_t size;
			if (!get_part_size_impl(stream_id, size))
			{
				err_code = ec_unknown_stream_part;
				return false;
			}

			v_data.resize(size);

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return true;
			}

			bool ret = get_part_impl(stream_id, v_data.data(), metadata);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_parts(const std::vector<int>& stream_ids, std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, std::vector<uint64_t>& metadatas)
		{
			std::lock_guard<std::mutex> lck(mtx);

			size_t no_streams = stream_ids.size();
			if (no_streams != v_datas.size() || no_streams != metadatas.size())
			{
				err_code = ec_wrong_part_desc;
				return false;
			}

			size_t size;

			std::vector<uint8_t*> datas;
			std::vector<size_t> data_sizes;

			datas.reserve(no_streams);
			data_sizes.reserve(no_streams);

			for (size_t i = 0; i < no_streams; ++i)
				if (get_part_size_impl(stream_ids[i], size))
				{
					v_datas.at(i).get().resize(size);
					datas.emplace_back(static_cast<uint8_t*>(v_datas.at(i).get().data()));
					data_sizes.emplace_back(size);
				}
				else
				{
					for (size_t j = 0; j < i; ++j)
						v_datas.at(j).get().clear();

					err_code = ec_unknown_stream_part;
					return false;
				}

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_parts_impl(stream_ids, datas, metadatas);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_part(const int stream_id, uint8_t* data, size_t& data_size, uint64_t& metadata)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!get_part_size_impl(stream_id, data_size))
			{
				err_code = ec_unknown_stream_part;
				return false;
			}

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_part_impl(stream_id, data, metadata);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_parts(const std::vector<int>& stream_ids, std::vector<uint8_t*> datas, std::vector<size_t>& data_sizes, std::vector<uint64_t>& metadatas)
		{
			std::lock_guard<std::mutex> lck(mtx);

			size_t no_streams = stream_ids.size();
			if (no_streams != datas.size() || no_streams != metadatas.size())
			{
				err_code = ec_unknown_stream_part;
				return false;
			}

			for (size_t i = 0; i < no_streams; ++i)
				if (get_part_size_impl(stream_ids[i], data_sizes[i]))
				{
					err_code = ec_unknown_stream_part;
					return false;
				}

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_parts_impl(stream_ids, datas, metadatas);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_part(const int stream_id, const int part_id, std::vector<uint8_t>& v_data, uint64_t& metadata)
		{
			std::lock_guard<std::mutex> lck(mtx);

			size_t size;
			if (!get_part_size_impl(stream_id, part_id, size))
			{
				err_code = ec_unknown_stream_part;
				return false;
			}

			v_data.resize(size);

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_part_impl(stream_id, part_id, v_data.data(), metadata);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_part(const int stream_id, const int part_id, uint8_t* data, size_t& data_size, uint64_t& metadata)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!get_part_size_impl(stream_id, part_id, data_size))
			{
				err_code = ec_unknown_stream_part;
				return false;
			}

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_part_impl(stream_id, part_id, data, metadata);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_sub_part(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, bool& last_sub_part, std::vector<uint8_t>& v_data, uint64_t& metadata)
		{
			std::lock_guard<std::mutex> lck(mtx);

			v_data.resize(max_sub_part_size);
			size_t sub_part_size;

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_sub_part_impl(stream_id, max_sub_part_size, part_id, sub_part_id, last_sub_part, v_data.data(), sub_part_size, metadata);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				v_data.resize(sub_part_size);
				err_code = ec_ok;
			}
			else
			{
				err_code = ec_file_read;
			}

			return ret;
		}

		// *******************************************************************************************
		bool get_sub_part(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, std::vector<uint8_t>& v_data, uint64_t& metadata)
		{
			bool to_ignore;

			return get_sub_part(stream_id, max_sub_part_size, part_id, sub_part_id, to_ignore, v_data, metadata);
		}

		// *******************************************************************************************
		bool get_sub_part(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, bool& last_part, uint8_t* v_data, size_t& sub_part_size, uint64_t& metadata)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_sub_part_impl(stream_id, max_sub_part_size, part_id, sub_part_id, last_part, v_data, sub_part_size, metadata);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_sub_part(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, uint8_t* v_data, size_t& sub_part_size, uint64_t& metadata)
		{
			bool to_ignore;

			return get_sub_part(stream_id, max_sub_part_size, part_id, sub_part_id, to_ignore, v_data, sub_part_size, metadata);
		}

		// *******************************************************************************************
		bool get_excerpt(const int stream_id, const size_t max_excerpt_size, const size_t offset, std::vector<uint8_t>& v_data)
		{
			std::lock_guard<std::mutex> lck(mtx);

			v_data.resize(max_excerpt_size);
			size_t excerpt_size;

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_excerpt_impl(stream_id, max_excerpt_size, offset, v_data.data(), excerpt_size);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			v_data.resize(excerpt_size);

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool get_excerpt(const int stream_id, const size_t max_excerpt_size, const size_t offset, uint8_t* v_data, size_t& excerpt_size)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!input->start_transaction())
			{
				err_code = ec_file_open;
				return false;
			}

			bool ret = get_excerpt_impl(stream_id, max_excerpt_size, offset, v_data, excerpt_size);

			if (!input->stop_transaction())
			{
				err_code = ec_file_close;
				return false;
			}

			if (ret)
			{
				err_code = ec_ok;
				return true;
			}
			else
			{
				err_code = ec_file_read;
				return false;
			}
		}

		// *******************************************************************************************
		bool rewind(const int stream_id)
		{
			std::lock_guard<std::mutex> lck(mtx);

			return rewind_impl(stream_id);
		}

		// *******************************************************************************************
		bool rewind(const std::vector<int>& stream_ids)
		{
			std::lock_guard<std::mutex> lck(mtx);

			for (auto x : stream_ids)
				if (x < 0 || x >= v_streams.size())
					return false;

			bool r = true;

			for (auto stream_id : stream_ids)
				r &= rewind_impl(stream_id);

			return r;
		}

		// *******************************************************************************************
		bool rewind_all()
		{
			std::lock_guard<std::mutex> lck(mtx);

			for (size_t i = 0; i < v_streams.size(); ++i)
				rewind_impl((int)i);

			return true;
		}

		// *******************************************************************************************
		size_t get_stream_metadata(const int stream_id)
		{
			std::lock_guard<std::mutex> lck(mtx);

			return v_streams[stream_id].metadata;
		}

		// *******************************************************************************************
		size_t get_no_parts(const int stream_id)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (stream_id < 0 || (size_t)stream_id >= v_streams.size())
				return 0;

			return v_streams[stream_id].parts.size();
		}

		// *******************************************************************************************
		void list_streams(std::vector<stream_t>& _v_streams)
		{
			std::lock_guard<std::mutex> lck(mtx);

			_v_streams = v_streams;
		}

		// *******************************************************************************************
		std::vector<std::string> get_params_string()
		{
			std::vector<std::string> info;

			info.push_back(params.str_archive_version + ": " + std::to_string(params.archive_version));
			info.push_back(params.str_parts_metadata_empty + ": " + (params.parts_metadata_empty ? "true" : "false"));
			info.push_back(params.str_parts_metadata_fixed_size + ": " + (params.parts_metadata_fixed_size ? "true" : "false"));

			return info;
		}
	};
}