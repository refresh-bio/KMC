#ifndef _ARCHIVE_H
#define _ARCHIVE_H

#include <cstdio>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <string>
#include <thread>
#include <mutex>
#include <functional>
#include <cinttypes>
#include <stdarg.h>

#include "buffered_io.h"
#include "unbuffered_io.h"

namespace refresh
{
	const uint32_t REFRESH_BUILD_ARCHIVE = 2;

template<class T_INPUT, class T_OUTPUT>
class archive_basic
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

		stream_t(const std::string &stream_name, size_t cur_id, uint64_t metadata, size_t total_size, size_t data_size) :
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
		uint64_t archive_version = 1;
		bool parts_metadata_fixed_size = false;
		bool parts_metadata_empty = false;

		static const inline std::string str_archive_version = "archive_version";
		static const inline std::string str_parts_metadata_fixed_size = "parts_metadata_fixed_size";
		static const inline std::string str_parts_metadata_empty = "parts_metadata_empty";
	} params;
	
private:
	bool input_mode;
	T_INPUT f_in;
	T_OUTPUT f_out;
	size_t io_buffer_size;

	size_t f_offset{};

	std::map<int, std::vector<std::pair<std::vector<uint8_t>, uint64_t>>> m_buffer;
	std::vector<stream_t> v_streams;
	std::unordered_map<std::string, size_t> rm_streams;

	std::mutex mtx;

	// *******************************************************************************************
	void unget_part(const int id)
	{
		if (v_streams[id].cur_id)
			--v_streams[id].cur_id;
	}

	// *******************************************************************************************
	size_t serialize_params()
	{
		size_t params_size = 0;

		std::map<std::string, uint64_t> dict;

		dict.emplace(params.str_archive_version, (uint64_t)params.archive_version);
		dict.emplace(params.str_parts_metadata_fixed_size, (uint64_t)params.parts_metadata_fixed_size);
		dict.emplace(params.str_parts_metadata_empty, (uint64_t)params.parts_metadata_empty);

		params_size += write(dict.size());
		for (const auto& par : dict)
		{
			params_size += write(par.first);
			params_size += write(par.second);
		}

		return params_size;
	}

	// *******************************************************************************************
	bool serialize()
	{
		size_t footer_size = 0;

		// Store stream part offsets
		footer_size += write(v_streams.size());

		for (auto& stream : v_streams)
		{
			size_t p = footer_size;

			footer_size += write(stream.stream_name);
			footer_size += write(stream.parts.size());	
			footer_size += write(stream.metadata);

			for (auto& part : stream.parts)
			{
				footer_size += write(part.offset);
				footer_size += write(part.size);
			}

			stream.total_size += footer_size - p;
		}

		// Optional params dictionary
		footer_size += serialize_params();

		write_fixed(footer_size);

		return true;
	}

	// *******************************************************************************************
	void deserialize_params()
	{
		std::map<std::string, uint64_t> dict;

		size_t dict_size;
		read(dict_size);
		std::string par_name;
		uint64_t par_val;

		// Load whole dictionary
		for (size_t i = 0; i < dict_size; ++i)
		{
			read(par_name);
			read(par_val);

			dict[par_name] = par_val;
		}

		// Set known parameters
		auto p = dict.find(params_t::str_archive_version);
		if (p != dict.end())
			params.archive_version = p->second;

		p = dict.find(params_t::str_parts_metadata_fixed_size);
		if (p != dict.end())
			params.parts_metadata_fixed_size = (bool) p->second;

		p = dict.find(params_t::str_parts_metadata_empty);
		if (p != dict.end())
			params.parts_metadata_empty = (bool) p->second;
	}

	// *******************************************************************************************
	bool deserialize()
	{
		size_t footer_size;
		size_t file_size_ = f_in.file_size();

		f_in.seek(file_size_ - 8ull);
		read_fixed(footer_size);

		f_in.seek(file_size_ - (size_t)(8 + footer_size));

		size_t readed = 0;

		// Read stream part offsets
		size_t n_streams;
		readed += read(n_streams);

		v_streams.resize(n_streams, stream_t());

		for (size_t i = 0; i < n_streams; ++i)
		{
			auto& stream_second = v_streams[i];

			readed += read(stream_second.stream_name);
			readed += read(stream_second.cur_id);
			readed += read(stream_second.metadata);

			stream_second.parts.resize(stream_second.cur_id);
			for (size_t j = 0; j < stream_second.cur_id; ++j)
			{
				readed += read(stream_second.parts[j].offset);
				readed += read(stream_second.parts[j].size);

				stream_second.total_size += stream_second.parts[j].size + 8;		// approximated (cannot check part metadata sizes w/o reading the archive)
				stream_second.data_size += stream_second.parts[j].size;
			}

			stream_second.cur_id = 0;

			rm_streams[stream_second.stream_name] = i;
		}

		// Optional params dictionary 
		if (readed < footer_size)
			deserialize_params();

		f_in.seek(0);

		return true;
	}

	// *******************************************************************************************
	template<typename T>
	size_t write_fixed(const T x)
	{
		f_out.write_uint(static_cast<uint64_t>(x), 8);

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

		f_out.put(no_bytes);

		for (int i = no_bytes; i; --i)
			f_out.put((x >> ((i - 1) * 8)) & 0xff);

		return no_bytes + 1;
	}

	// *******************************************************************************************
	size_t write(const std::string& s)
	{
//		f_out.write(s);
		f_out.write((uint8_t*)s.c_str(), s.size());
		f_out.put(0);

		return s.size() + 1;
	}

	// *******************************************************************************************
	template<typename T>
	size_t read_fixed(T& x)
	{
		x = static_cast<T>(f_in.read_uint(8));

		return 8;
	}

	// *******************************************************************************************
	size_t read(std::string& s)
	{
		s.clear();

		while (true)
		{
			int c = f_in.get();
			if (c == EOF)
				return 0;

			if (c == 0)
				return s.size() + 1;

			s.push_back((char)c);
		}

		return 0;
	}

	// *******************************************************************************************
	template<typename T>
	size_t read(T& x)
	{
		int no_bytes = f_in.get();

		x = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			x <<= 8;
			x += static_cast<T>(f_in.get());
		}

		return no_bytes + 1;
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
	int get_stream_id_impl(const std::string& stream_name)
	{
		auto p = rm_streams.find(stream_name);
		if (p != rm_streams.end())
			return (int)p->second;

		return -1;
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
	bool add_part_impl(const int stream_id, const uint8_t* data, const size_t data_size, const uint64_t metadata)
	{
		v_streams[stream_id].parts.push_back(part_t(f_offset, data_size));

		if (params.parts_metadata_empty)
			;
		else if(params.parts_metadata_fixed_size)
			f_offset += write_fixed(metadata);
		else
			f_offset += write(metadata);

		f_out.write(data, data_size);

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
		v_streams[stream_id].parts[part_id] = part_t(f_offset, v_data.size());

		if (params.parts_metadata_empty)
			;
		else if (params.parts_metadata_fixed_size)
			f_offset += write_fixed(metadata);
		else
			f_offset += write(metadata);

		f_out.write(v_data.data(), v_data.size());

		f_offset += v_data.size();

		v_streams[stream_id].total_size += f_offset - v_streams[stream_id].parts[part_id].offset;
		v_streams[stream_id].data_size += v_data.size();

		return true;
	}

	// *******************************************************************************************
	bool get_part_impl(const int stream_id, uint8_t* data, uint64_t& metadata)
	{
		auto& p = v_streams[stream_id];

		f_in.seek(p.parts[p.cur_id].offset);

		if (p.parts[p.cur_id].size != 0)
		{
			if (params.parts_metadata_empty)
				metadata = 0;
			else if(params.parts_metadata_fixed_size)
				read_fixed(metadata);
			else
				read(metadata);
		}
		else
		{
			metadata = 0;
			p.cur_id++;
			return true;
		}

		f_in.read(data, p.parts[p.cur_id].size);

		p.cur_id++;

		p.sub_part_id = 0;
		p.sub_part_offset = 0;
		p.sub_part_metadata = 0;
		p.sub_part_metadata_size = 0;

		//if (r != p.parts[p.cur_id-1].size)
			//return false;

		//return r == p.parts[p.cur_id-1].size;
		return true;
	}

	// *******************************************************************************************
	bool get_part_impl(const int stream_id, const int part_id, uint8_t *data, uint64_t& metadata)
	{
		auto& p = v_streams[stream_id];

		if ((size_t)part_id >= p.parts.size())
			return false;

		f_in.seek(p.parts[part_id].offset);

		if (p.parts[part_id].size != 0)
		{
			if (params.parts_metadata_empty)
				metadata = 0;
			else if (params.parts_metadata_fixed_size)
				read_fixed(metadata);
			else
				read(metadata);
		}
		else
		{
			metadata = 0;
			return true;
		}

		f_in.read(data, p.parts[part_id].size);

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
	bool get_parts_impl(const std::vector<int>& stream_ids, std::vector<uint8_t*> datas, std::vector<size_t>& data_sizes, std::vector<uint64_t>& metadatas)
	{
		size_t no_streams = stream_ids.size();

		for (size_t i = 0; i < no_streams; ++i)
			if (!get_part_impl(stream_ids[i], datas[i], metadatas.at(i)))
			{
				for (size_t j = 0; j < i; ++j)
					unget_part(j);

				return false;
			}

		return true;
	}

	// *******************************************************************************************
	bool get_sub_part_impl(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, bool &last_sub_part, uint8_t* data, size_t& sub_part_size, uint64_t& metadata)
	{
		auto& p = v_streams[stream_id];

		if (p.cur_id >= p.parts.size())
			return false;

		part_id = static_cast<int>(p.cur_id);
		sub_part_id = static_cast<int>(p.sub_part_id++);

		if (sub_part_id == 0)
		{
			f_in.seek(p.parts[p.cur_id].offset);
			p.sub_part_offset = 0;

			if (p.parts[p.cur_id].size != 0)
			{
				if (params.parts_metadata_empty)
					p.sub_part_metadata_size = 0;
				else if(params.parts_metadata_fixed_size)
					p.sub_part_metadata_size = read_fixed(p.sub_part_metadata);
				else
					p.sub_part_metadata_size = read(p.sub_part_metadata);
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
			f_in.seek(p.parts[p.cur_id].offset + p.sub_part_metadata_size + p.sub_part_offset);
		}

		sub_part_size = std::min(max_sub_part_size, p.parts[p.cur_id].size - p.sub_part_offset);

		f_in.read(data, sub_part_size);
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
					f_in.seek(part.offset);
					metadata_len = f_in.get() + 1;
				}
				else
					metadata_len = 8;

				size_t local_offset = offset - acc_part_size;
				excerpt_size = std::min<size_t>(max_excerpt_size, part.size - local_offset);

				f_in.seek(part.offset + metadata_len + local_offset);
				f_in.read(v_data, excerpt_size);

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
	bool flush_out_buffers_imp()
	{
		for (auto& x : m_buffer)
			for (auto& y : x.second)
				add_part(x.first, y.first, y.second);

		m_buffer.clear();

		return true;
	}

public:
	// *******************************************************************************************
	archive_basic(const bool _input_mode, const size_t _io_buffer_size = 64 << 20)
	{
		input_mode = _input_mode;
		io_buffer_size = _io_buffer_size;
	}

	// *******************************************************************************************
	~archive_basic()
	{
		close();
	}

	// *******************************************************************************************
	bool open(const std::string& file_name, const params_t &_params = params_t()) 
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (f_in.opened())
			f_in.close();
		if (f_out.opened())
			f_out.close();

		if (input_mode)
			f_in.open(file_name, io_buffer_size);
		else
			f_out.open(file_name, io_buffer_size);

		if (!f_in.opened() && !f_out.opened())
			return false;

		if (input_mode)
			deserialize();
		else
			params = _params;

		f_offset = 0;

		return true;
	}

	// *******************************************************************************************
	bool close()
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (!f_in.opened() && !f_out.opened())
			return false;

		if (input_mode)
			f_in.close();
		else
		{
			flush_out_buffers_imp();
			serialize();
			f_out.close();
		}

		return true;
	}

	// *******************************************************************************************
	int register_stream(const std::string& stream_name)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return register_stream_impl(stream_name);
	}

	// *******************************************************************************************
	std::vector<int> register_streams(const std::vector<std::string> &stream_names)
	{
		std::lock_guard<std::mutex> lck(mtx);

		std::vector<int> ret;

		ret.reserve(stream_names.size());

		for (const auto& sn : stream_names)
			ret.emplace_back(register_stream_impl(sn));

		return ret;
	}

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
	size_t get_stream_total_size(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return get_stream_total_size_impl(stream_id);
	}

	// *******************************************************************************************
	std::vector<size_t> get_stream_total_sizes(const std::vector<int> &stream_ids)
	{
		std::lock_guard<std::mutex> lck(mtx);

		std::vector<size_t> ret;

		ret.reserve(stream_ids.size());

		for(const int sid : stream_ids)
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
	std::vector<size_t> get_stream_data_sizes(const std::vector<int> &stream_ids)
	{
		std::lock_guard<std::mutex> lck(mtx);

		std::vector<size_t> ret;

		ret.reserve(stream_ids.size());

		for (const int sid : stream_ids)
			ret.emplace_back(get_stream_data_size_impl(sid));

		return ret;
	}

	// *******************************************************************************************
	bool add_part(const int stream_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return add_part_impl(stream_id, v_data.data(), v_data.size(), metadata);
	}

	// *******************************************************************************************
	bool add_part(const int stream_id, const uint8_t* data, const size_t data_size, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return add_part_impl(stream_id, data, data_size, metadata);
	}

	// *******************************************************************************************
	// metadatas can be empty!
	bool add_parts(const std::vector<int> &stream_ids, const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, const std::vector<uint64_t> &metadatas)
	{
		std::lock_guard<std::mutex> lck(mtx);

		auto n_items = stream_ids.size();

		if (v_datas.size() != n_items)
			return false;

		if (!metadatas.empty() && metadatas.size() != n_items)
			return false;

		bool ret = true;

		for(size_t i = 0; i < n_items; ++i)
			ret &= add_part_impl(stream_ids[i], v_datas.at(i).get().data(), v_datas.at(i).get().size(), metadatas.empty() ? 0 : metadatas[i]);

		return ret;
	}

	// *******************************************************************************************
	// metadatas can be empty!
	bool add_parts(const std::vector<int> &stream_ids, const std::vector<uint8_t *> & datas, const std::vector<size_t> data_sizes, const std::vector<uint64_t> & metadatas)
	{
		std::lock_guard<std::mutex> lck(mtx);

		auto n_items = stream_ids.size();

		if (datas.size() != n_items)
			return false;

		if (data_sizes.size() != n_items)
			return false;

		if (!metadatas.empty() && metadatas.size() != n_items)
			return false;

		bool ret = true;

		for (size_t i = 0; i < n_items; ++i)
			ret &= add_part_impl(stream_ids[i], datas[i], data_sizes[i], metadatas.empty() ? 0 : metadatas[i]);

		return ret;
	}

	// *******************************************************************************************
	int add_part_prepare(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return add_part_prepare_impl(stream_id);
	}

	// *******************************************************************************************
	std::vector<int> add_parts_prepare(const std::vector<int> &stream_ids)
	{
		std::lock_guard<std::mutex> lck(mtx);

		std::vector<int> ret;

		ret.reserve(stream_ids.size());

		for(const auto sid : stream_ids)
			ret.emplace_back(add_part_prepare_impl(sid));

		return ret;
	}

	// *******************************************************************************************
	bool add_part_complete(const int stream_id, const int part_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return add_part_complete_impl(stream_id, part_id, v_data, metadata);
	}

	// *******************************************************************************************
	// metadatas can be empty
	bool add_parts_complete(const std::vector<int> &stream_ids, const std::vector<int> &part_ids, const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, const std::vector<uint64_t> &metadatas)
	{
		std::lock_guard<std::mutex> lck(mtx);

		auto n_items = stream_ids.size();

		if (part_ids.size() != n_items)
			return false;

		if (v_datas.size() != n_items)
			return false;

		if (!metadatas.empty() && metadatas.size() != n_items)
			return false;

		bool ret = true;

		for(size_t i = 0; i < n_items; ++i)
			ret &= add_part_complete_impl(stream_ids[i], part_ids[i], v_datas.at(i), metadatas.empty() ? 0 : metadatas[i]);

		return ret;
	}

	// *******************************************************************************************
	bool add_part_buffered(const int stream_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		m_buffer[stream_id].emplace_back(v_data, metadata);

		return true;
	}

	// *******************************************************************************************
	bool add_parts_buffered(const std::vector<int> &stream_ids, const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, const std::vector<uint64_t> &metadatas)
	{
		std::lock_guard<std::mutex> lck(mtx);

		auto n_items = stream_ids.size();

		if (v_datas.size() != n_items)
			return false;

		if (!metadatas.empty() && metadatas.size() != n_items)
			return false;

		for(size_t i = 0; i < n_items; ++i)
			m_buffer[stream_ids[i]].emplace_back(v_datas.at(i), metadatas.empty() ? 0 : metadatas[i]);

		return true;
	}

	// *******************************************************************************************
	bool flush_out_buffers()
	{
		std::lock_guard<std::mutex> lck(mtx);

		return flush_out_buffers_imp();
	}

	// *******************************************************************************************
	bool get_part_size(const int stream_id, size_t& size)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return get_part_size_impl(stream_id, size);
	}
	// *******************************************************************************************
	bool get_part_size(const int stream_id, const int part_id, size_t& size)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return get_part_size_impl(stream_id, part_id, size);
	}

	// *******************************************************************************************
	bool get_part_sizes(const std::vector<int> &stream_ids, std::vector<size_t>& sizes)
	{
		std::lock_guard<std::mutex> lck(mtx);

		auto n_items = stream_ids.size();

		sizes.resize(n_items);

		bool ret = true;

		for (size_t i = 0; i < n_items; ++i)
			ret &= get_part_size_impl(stream_ids[i], sizes[i]);

		return ret;
	}

	// *******************************************************************************************
	bool get_part(const int stream_id, std::vector<uint8_t>& v_data, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		size_t size;
		if (!get_part_size_impl(stream_id, size))
			return false;

		v_data.resize(size);

		return get_part_impl(stream_id, v_data.data(), metadata);
	}

	// *******************************************************************************************
	bool get_parts(const std::vector<int> &stream_ids, std::vector<std::reference_wrapper<std::vector<uint8_t>>>& v_datas, std::vector<uint64_t>& metadatas)
	{
		std::lock_guard<std::mutex> lck(mtx);

		size_t no_streams = stream_ids.size();
		if (no_streams != v_datas.size() || no_streams != metadatas.size())
			return false;

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
				return false;
			}

		return get_parts_impl(stream_ids, datas, data_sizes, metadatas);
	}

	// *******************************************************************************************
	bool get_part(const int stream_id, uint8_t* data, size_t& data_size, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (!get_part_size_impl(stream_id, data_size))
			return false;

		return get_part_impl(stream_id, data, metadata);
	}

	// *******************************************************************************************
	bool get_parts(const std::vector<int>& stream_ids, std::vector<uint8_t*> datas, std::vector<size_t>& data_sizes, std::vector<uint64_t>& metadatas)
	{
		std::lock_guard<std::mutex> lck(mtx);

		size_t no_streams = stream_ids.size();
		if (no_streams != datas.size() || no_streams != metadatas.size())
			return false;

		return get_parts_impl(stream_ids, datas, data_sizes, metadatas);
	}

	// *******************************************************************************************
	bool get_part(const int stream_id, const int part_id, std::vector<uint8_t>& v_data, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		size_t size;
		if (!get_part_size_impl(stream_id, part_id, size))
			return false;

		v_data.resize(size);

		return get_part_impl(stream_id, part_id, v_data.data(), metadata);
	}

	// *******************************************************************************************
	bool get_part(const int stream_id, const int part_id, uint8_t* data, size_t& data_size, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (!get_part_size_impl(stream_id, part_id, data_size))
			return false;

		return get_part_impl(stream_id, part_id, data, metadata);
	}

	// *******************************************************************************************
	bool get_sub_part(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, bool &last_sub_part, std::vector<uint8_t>& v_data, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		v_data.resize(max_sub_part_size);
		size_t sub_part_size;

		bool r = get_sub_part_impl(stream_id, max_sub_part_size, part_id, sub_part_id, last_sub_part, v_data.data(), sub_part_size, metadata);

		if(r)
			v_data.resize(sub_part_size);

		return r;
	}

	// *******************************************************************************************
	bool get_sub_part(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, std::vector<uint8_t> &v_data, uint64_t& metadata)
	{
		bool to_ignore;
		return get_sub_part(stream_id, max_sub_part_size, part_id, sub_part_id, to_ignore, v_data, metadata);
	}

	// *******************************************************************************************
	bool get_sub_part(const int stream_id, const size_t max_sub_part_size, int& part_id, int& sub_part_id, bool &last_part, uint8_t* v_data, size_t& sub_part_size, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return get_sub_part_impl(stream_id, max_sub_part_size, part_id, sub_part_id, last_part, v_data, sub_part_size, metadata);
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

		bool r = get_excerpt_impl(stream_id, max_excerpt_size, offset, v_data.data(), excerpt_size);

		v_data.resize(excerpt_size);

		return r;
	}

	// *******************************************************************************************
	bool get_excerpt(const int stream_id, const size_t max_excerpt_size, const size_t offset, uint8_t* v_data, size_t &excerpt_size)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return get_excerpt_impl(stream_id, max_excerpt_size, offset, v_data, excerpt_size);
	}

	// *******************************************************************************************
	bool rewind(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return rewind_impl(stream_id);
	}

	// *******************************************************************************************
	bool rewind(const std::vector<int> &stream_ids)
	{
		std::lock_guard<std::mutex> lck(mtx);

		for (auto x : stream_ids)
			if (x < 0 || x >= v_streams.size())
				return false;

		bool r = true;

		for(auto stream_id : stream_ids)
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
	void set_stream_metadata(const int stream_id, const uint64_t metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		v_streams[stream_id].metadata = metadata;
	}

	// *******************************************************************************************
	// void set_raw_sizes(const int stream_id, const size_t raw_size)

	// *******************************************************************************************
	size_t get_stream_metadata(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return v_streams[stream_id].metadata;
	}

	// *******************************************************************************************
	// size_t get_raw_sizes(const int stream_id)

	// *******************************************************************************************
	size_t get_no_streams()
	{
		std::lock_guard<std::mutex> lck(mtx);

		return v_streams.size();
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
	// size_t get_no_parts(const int stream_id)

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

using archive_buffered_io = archive_basic<buffered_in_file, buffered_out_file>;
using archive_unbuffered_io = archive_basic<unbuffered_in_file, unbuffered_out_file>;

}
// EOF
#endif