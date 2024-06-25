#pragma once

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <vector>

#include "input_common.h"
#include "io_utils.h"

namespace refresh
{
	namespace io
	{
		// *******************************************************************************************
		// Common class for inputs that are all buffered
		// *******************************************************************************************
		class input_all_buffered : public input_common
		{
		protected:
			std::vector<uint8_t> buffer;
		
		public:
			// *******************************************************************************************
			input_all_buffered() : input_common()
			{
			}
		
			// *******************************************************************************************
			virtual ~input_all_buffered()
			{}

			// *******************************************************************************************
			virtual bool close() final
			{
				active = false;

				buffer.clear();
				buffer.shrink_to_fit();

				return true;
			}

			// *******************************************************************************************
			virtual bool get(uint8_t& x) final
			{
				if (f_pos >= f_size)
					return false;

				x = buffer[f_pos++];

				return true;
			}

			// *******************************************************************************************
			virtual bool read_uint(const int no_bytes, uint64_t& x) final
			{
				if (f_pos + no_bytes > f_size)
				{
					f_pos = f_size;
					return false;
				}

				uint64_t shift = 0;
				x = 0;

				for (int i = 0; i < no_bytes; ++i, shift += 8)
					x += ((uint64_t)buffer[f_pos++]) << shift;

				return true;
			}

			// *******************************************************************************************
			virtual size_t read(uint8_t* ptr, size_t size) final
			{
				if (f_pos + size > f_size)
					size = f_size - f_pos;

				std::copy_n(buffer.data() + f_pos, size, ptr);
				f_pos += size;

				return size;
			}

			// *******************************************************************************************
			virtual bool eof() const final
			{
				return f_pos >= f_size;
			}

			// *******************************************************************************************
			virtual bool seek(const size_t requested_pos) final
			{
				if (requested_pos > f_size)
					return false;

				f_pos = requested_pos;

				return true;
			}
		};

		// *******************************************************************************************
		// 
		// *******************************************************************************************
		class input_stdin : public input_all_buffered
		{
			// *******************************************************************************************
			void fill_buffer()
			{
				const size_t init_size = 64ull << 10;
				const double mult_factor = 1.4;

				buffer.resize(init_size);
				size_t buffer_filled = 0;

				while (true)
				{
					size_t to_read = buffer.size() - buffer_filled;

					size_t readed = refresh::io::fread(buffer.data() + buffer_filled, 1, to_read, stdin);

					buffer_filled += readed;

					if (buffer_filled != buffer.size())
					{
						buffer.resize(buffer_filled);
						break;
					}

					to_read = (size_t)(to_read * mult_factor);
				}

				f_size = buffer.size();
				f_pos = 0;
			}

		public:
			// *******************************************************************************************
			input_stdin() : input_all_buffered()
			{
#ifdef _WIN32
				int res = _setmode(_fileno(stdin), _O_BINARY);
				active = res != -1;
#endif

				if(active)
					fill_buffer();
			}

			// *******************************************************************************************
			virtual ~input_stdin()
			{
				// Nothing to do
			}
		};

		// *******************************************************************************************
		// 
		// *******************************************************************************************
		class input_memory : public input_all_buffered
		{
		public:
			// *******************************************************************************************
			input_memory(const std::vector<uint8_t> &data) : input_all_buffered()
			{
				buffer = data;
				f_size = buffer.size();
				f_pos = 0;

				active = true;
			}

			// *******************************************************************************************
			input_memory(std::vector<uint8_t> &&data) : input_all_buffered()
			{
				buffer = move(data);
				f_size = buffer.size();
				f_pos = 0;

				active = true;
			}

			// *******************************************************************************************
			virtual ~input_memory()
			{
				// Nothing to do
			}
		};

		// *******************************************************************************************
		// 
		// *******************************************************************************************
		class input_file_all_buffered : public input_all_buffered
		{
			// *******************************************************************************************
			void fill_buffer(const std::string &file_name)
			{
				std::filesystem::path fp(file_name);

				if (!std::filesystem::exists(fp))
					return;

				buffer.resize(std::filesystem::file_size(fp));
				FILE* f = fopen(file_name.c_str(), "rb");

				if (!f)
					return;

				active = refresh::io::fread(buffer.data(), 1, buffer.size(), f) == buffer.size();

				f_pos = 0;
				f_size = buffer.size();

				fclose(f);
			}

		public:
			// *******************************************************************************************
			input_file_all_buffered(const std::string &file_name) : input_all_buffered()
			{
				fill_buffer(file_name);
			}

			// *******************************************************************************************
			virtual ~input_file_all_buffered()
			{
				// Nothing to do
			}
		};
	}
}