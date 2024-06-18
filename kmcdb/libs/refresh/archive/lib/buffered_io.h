#ifndef _BUFFERED_IO_H
#define _BUFFERED_IO_H

#include <algorithm>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <cinttypes>
#include <iostream>

#include "archive_helper.h"

namespace refresh
{
// *******************************************************************************************
// Buffered input file
class buffered_in_file
{
	size_t BUFFER_SIZE = 0;

	FILE* f;
	uint8_t* buffer;
	size_t buffer_pos;
	size_t buffer_filled;

	size_t file_size_;
	size_t before_buffer_bytes;

	bool use_stdin;

public:
	// *******************************************************************************************
	buffered_in_file() : f(nullptr), buffer(nullptr), buffer_pos(0), buffer_filled(0), file_size_(0), before_buffer_bytes(0), use_stdin(false)
	{};

	// *******************************************************************************************
	~buffered_in_file()
	{
		close();
	}

	// *******************************************************************************************
	bool open(const std::string& file_name, const size_t _BUFFER_SIZE = 128 << 20)
	{
		if (f)
			return false;

		use_stdin = file_name.empty();
		if (use_stdin)
		{
			f = stdin;
#ifdef _WIN32
			int res = _setmode(_fileno(f), _O_BINARY);
			if (res == -1)
				return false;
#endif
		}
		else
		{
			f = fopen(file_name.c_str(), "rb");
			if (!f)
				return false;
		}

		refresh::io::fseek(f, 0, SEEK_END);
		file_size_ = refresh::io::ftell(f);
		refresh::io::fseek(f, 0, SEEK_SET);
		before_buffer_bytes = 0;

		if (_BUFFER_SIZE == ~0ull)
			BUFFER_SIZE = file_size_;
		else
			BUFFER_SIZE = _BUFFER_SIZE;

		buffer = new uint8_t[BUFFER_SIZE];

		buffer_pos = 0;

		buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);

		return true;
	}

	// *******************************************************************************************
	bool close()
	{
		if (f)
		{
			if(!use_stdin)
				fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return true;
	}

	// *******************************************************************************************
	bool opened()
	{
		return f != nullptr;
	}

	// *******************************************************************************************
	int get()
	{
		if (buffer_pos < buffer_filled)
			return buffer[buffer_pos++];

		if (feof(f))
			return EOF;

		before_buffer_bytes += buffer_filled;

		buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
		if (buffer_filled == 0)
			return EOF;

		buffer_pos = 0;

		return buffer[buffer_pos++];
	}

	// *******************************************************************************************
	bool un_get()
	{
		if (buffer_pos)
		{
			--buffer_pos;
			return true;
		}

		return false;
	}

	// *******************************************************************************************
	uint64_t read_uint(const int no_bytes)
	{
		uint64_t x = 0;
		uint64_t shift = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			uint64_t c = get();
			x += c << shift;
			shift += 8;
		}

		return x;
	}

	// *******************************************************************************************
	void read(uint8_t* ptr, size_t size)
	{
		if (before_buffer_bytes + buffer_pos + size > file_size_)
			size = file_size_ - (before_buffer_bytes + buffer_pos);

		size_t to_read = size;

		while (buffer_pos + to_read > BUFFER_SIZE)
		{
			memcpy(ptr, buffer + buffer_pos, BUFFER_SIZE - buffer_pos);
			ptr += BUFFER_SIZE - buffer_pos;
			to_read -= BUFFER_SIZE - buffer_pos;

			before_buffer_bytes += buffer_filled;
			buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
			buffer_pos = 0;
		}

		memcpy(ptr, buffer + buffer_pos, to_read);
		buffer_pos += to_read;
	}

	// *******************************************************************************************
	bool eof() const
	{
		return before_buffer_bytes + buffer_pos >= file_size_;
	}
	
	// *******************************************************************************************
	bool seek(const size_t requested_pos)
	{
		if (requested_pos >= before_buffer_bytes && requested_pos < before_buffer_bytes + buffer_filled)
			buffer_pos = requested_pos - before_buffer_bytes;
		else
		{
			before_buffer_bytes = requested_pos;
			refresh::io::fseek(f, static_cast<long long>(requested_pos), SEEK_SET);
			buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
			buffer_pos = 0;
		}

		return true;
	}
	
	// *******************************************************************************************
	size_t file_size() const
	{
		if (f)
			return file_size_;
		else
			return 0;
	}

	// *******************************************************************************************
	size_t get_pos() const
	{
		return before_buffer_bytes + buffer_pos;
	}
};

// *******************************************************************************************
// Buffered output file
class buffered_out_file
{
	size_t BUFFER_SIZE;

	FILE* f;
	uint8_t* buffer;
	size_t buffer_pos;
	bool success;
	bool use_stdout;

public:
	// *******************************************************************************************
	buffered_out_file() : BUFFER_SIZE(8u << 20), f(nullptr), buffer(nullptr), buffer_pos(0), success(false), use_stdout(false)
	{};

	// *******************************************************************************************
	~buffered_out_file()
	{
		close();
	}

	// *******************************************************************************************
	bool open(const std::string& file_name, const size_t _BUFFER_SIZE = 8 << 20)
	{
		if (f)
			return false;

		use_stdout = file_name.empty();

		if (use_stdout)
		{
			f = stdout;
#ifdef _WIN32
			int res = _setmode(_fileno(f), _O_BINARY);
			if (res == -1)
				return false;
#endif
		}
		else
		{
			f = fopen(file_name.c_str(), "wb");
			if (!f)
				return false;
		}

		BUFFER_SIZE = _BUFFER_SIZE;
		buffer = new uint8_t[BUFFER_SIZE];
		buffer_pos = 0;
		success = true;

		return true;
	}

	// *******************************************************************************************
	bool close()
	{
		if (buffer_pos)
		{
			success &= fwrite(buffer, 1, buffer_pos, f) == buffer_pos;
			buffer_pos = 0;
		}

		if (f)
		{
			fflush(f);

			if(!use_stdout)
				fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return success;
	}

	// *******************************************************************************************
	bool opened()
	{
		return f != nullptr;
	}

	// *******************************************************************************************
	void put(const uint8_t c)
	{
		if (buffer_pos == BUFFER_SIZE)
		{
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;
			buffer_pos = 0;
		}

		buffer[buffer_pos++] = c;
	}

	// *******************************************************************************************
	void write(const uint8_t* p, size_t n)
	{
		uint8_t* q = (uint8_t*)p;

		while (buffer_pos + n > BUFFER_SIZE)
		{
			size_t small_n = BUFFER_SIZE - buffer_pos;
			memcpy(buffer + buffer_pos, q, small_n);
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;

			buffer_pos = 0;
			n -= small_n;
			q += small_n;
		}

		memcpy(buffer + buffer_pos, q, n);
		buffer_pos += n;
	}

	// *******************************************************************************************
	void write_uint(const uint64_t _x, const int no_bytes)
	{
		uint64_t x = _x;

		for (int i = 0; i < no_bytes; ++i)
		{
			put(x & 0xff);
			x >>= 8;
		}
	}

	// *******************************************************************************************
/*	void write(const std::string& s)
	{
		write((uint8_t*)s.c_str(), s.size());
	}

	// *******************************************************************************************
	void write(const std::string& s, const size_t start_pos, const size_t len)
	{
		write((uint8_t*)s.c_str() + start_pos, len);
	}*/
};
}
// EOF
#endif