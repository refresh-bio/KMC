#pragma once

#include "input_common.h"

namespace refresh
{
	namespace io
	{
		class input_file_buffered : public input_common
		{
		protected:
			FILE* f{};

			size_t buffer_size{};
			std::vector<uint8_t> buffer;
			size_t buffer_pos{};
			size_t buffer_filled{};
			size_t before_buffer_bytes{};

		public:
			// *******************************************************************************************
			input_file_buffered(const std::string& file_name, const size_t _buffer_size = 32 << 20) : 
				input_common(),
				buffer_size(_buffer_size)
			{
				std::filesystem::path fp(file_name);

				if (!std::filesystem::exists(fp))
					return;

				f_size = std::filesystem::file_size(fp);
				f_pos = 0;

				if (f_size > buffer_size)
					buffer_size = f_size;

				f = fopen(file_name.c_str(), "rb");

				active = f != nullptr;

				buffer.resize(buffer_size);
				buffer_filled = refresh::io::fread(buffer.data(), 1, buffer.size(), f);
				buffer_pos = 0;
				before_buffer_bytes = 0;
			}

			// *******************************************************************************************
			virtual ~input_file_buffered()
			{
				close();
			}

			// *******************************************************************************************
			virtual bool close()
			{
				if (f)
				{
					if (std::fclose(f) != 0)
						return false;
					f = nullptr;
				}

				active = false;

				return true;
			}

			// *******************************************************************************************
			virtual bool get(uint8_t& x)
			{
				if (buffer_pos < buffer_filled)
				{
					x = buffer[buffer_pos++];
					++f_pos;
					return true;
				}
					
				if (std::feof(f))
					return false;

				before_buffer_bytes += buffer_filled;

				buffer_filled = refresh::io::fread(buffer.data(), 1, buffer.size(), f);
				if (buffer_filled == 0)
					return false;

				buffer_pos = 0;

				x = buffer[buffer_pos++];
				++f_pos;

				return true;
			}
			
			// *******************************************************************************************
			virtual bool read_uint(const int no_bytes, uint64_t& x)
			{
				uint64_t shift = 0;
				x = 0;

				for (int i = 0; i < no_bytes; ++i)
				{
					uint8_t c;
					if (!get(c))
						return false;

					x += ((uint64_t) c) << shift;
					shift += 8;
				}

				return true;
			}
			
			// *******************************************************************************************
			virtual size_t read(uint8_t* ptr, size_t size)
			{
				if (before_buffer_bytes + buffer_pos + size > f_size)
					size = f_size - (before_buffer_bytes + buffer_pos);

				size_t to_read = size;
				uint8_t *ptr0 = ptr;

				while (buffer_pos + to_read > buffer_size)
				{
					memcpy(ptr, buffer.data() + buffer_pos, buffer_size - buffer_pos);
					ptr += buffer_size - buffer_pos;
					to_read -= buffer_size - buffer_pos;

					before_buffer_bytes += buffer_filled;
					buffer_filled = refresh::io::fread(buffer.data(), 1, buffer.size(), f);

					if (buffer_filled == 0)
					{
						f_pos += ptr - ptr0;
						return ptr - ptr0;
					}

					buffer_pos = 0;
				}

				memcpy(ptr, buffer.data() + buffer_pos, to_read);
				buffer_pos += to_read;
				ptr += to_read;

				f_pos += ptr - ptr0;

				return ptr - ptr0;
			}
			
			// *******************************************************************************************
			virtual bool eof() const
			{
				return before_buffer_bytes + buffer_pos >= f_pos;
			}
			
			// *******************************************************************************************
			virtual bool seek(const size_t requested_pos)
			{
				if (requested_pos > f_size)
					return false;

				if (requested_pos >= before_buffer_bytes && requested_pos < before_buffer_bytes + buffer_filled)
					buffer_pos = requested_pos - before_buffer_bytes;
				else
				{
					before_buffer_bytes = requested_pos;
					if (refresh::io::fseek(f, static_cast<long long>(requested_pos), SEEK_SET) != 0)
						return false;
					buffer_filled = refresh::io::fread(buffer.data(), 1, buffer.size(), f);
					buffer_pos = 0;
				}

				f_pos = requested_pos;

				return true;
			}
		};

		// *******************************************************************************************
		//
		// *******************************************************************************************
		class input_file_buffered_reopen : public input_file_buffered
		{
			std::string file_name;
			bool in_transaction;

		public:
			// *******************************************************************************************
			input_file_buffered_reopen(const std::string& file_name, const size_t _buffer_size = 32 << 20) :
				input_file_buffered(file_name, _buffer_size),
				file_name(file_name),
				in_transaction(false)
			{
				if (f)		// In reopen-mode file is closed after the openning
				{
					fclose(f);
					f = nullptr;
				}
			}

			// *******************************************************************************************
			virtual ~input_file_buffered_reopen()
			{}

			// *******************************************************************************************
			virtual bool start_transaction()
			{
				if (in_transaction)
					return false;

				assert(f == nullptr);

				f = fopen(file_name.c_str(), "rb");
				if (!f)
					return false;

				in_transaction = true;

				seek(f_pos);

				return true;
			}

			// *******************************************************************************************
			virtual bool stop_transaction()
			{
				if (!in_transaction)
					return false;

				assert(f != nullptr);

				in_transaction = false;
				if (std::fclose(f) != 0)
					return false;
				f = nullptr;

				return true;
			}
		};
	}
}