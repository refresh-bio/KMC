#pragma once

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "output_common.h"
#include "io_utils.h"

#include <vector>

namespace refresh
{
	namespace io
	{
		// *******************************************************************************************
		//
		// *******************************************************************************************
		class output_memory : public output_common
		{
			const size_t init_buffer_size = 1 << 20;
			const double mult_factor = 1.4;

			std::vector<uint8_t> &buffer;

			// *******************************************************************************************
			void check_extra_space(const size_t n_bytes)
			{
				if (buffer.size() + n_bytes <= buffer.capacity())
					return;

				size_t new_size = buffer.size() + n_bytes;
				new_size = (size_t)((double) new_size * mult_factor);

				buffer.reserve(new_size);
			}

		public:
			// *******************************************************************************************
			output_memory(std::vector<uint8_t>& _buffer) :
				output_common(),
				buffer(_buffer)
			{
				buffer.clear();
				buffer.reserve(init_buffer_size);
				
				active = true;
			}

			// *******************************************************************************************
			virtual ~output_memory()
			{}

			// *******************************************************************************************
			virtual bool close() final
			{
				buffer.shrink_to_fit();

				return true;
			}

			// *******************************************************************************************
			virtual bool put(const uint8_t c) final
			{
				check_extra_space(1);

				buffer.emplace_back(c);
				++f_size;
				
				return true;
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) final
			{
				check_extra_space(n_bytes);

				buffer.insert(buffer.end(), ptr, ptr + n_bytes);
				f_size += n_bytes;

				return n_bytes;
			}

			// *******************************************************************************************
			virtual size_t write_uint(const uint64_t val, const size_t n_bytes) final
			{
				check_extra_space(n_bytes);

				uint64_t x = val;

				for (size_t i = 0; i < n_bytes; ++i)
				{
					put((uint8_t)(x & 0xffull));
					x >>= 8;
				}

				f_size += n_bytes;

				return n_bytes;
			}
		};

		// *******************************************************************************************
		//
		// *******************************************************************************************
		class output_stdout : public output_common
		{
		public:
			// *******************************************************************************************
			output_stdout() :
				output_common()
			{
#ifdef _WIN32
				int res = _setmode(_fileno(stdout), _O_BINARY);
				active = res != -1;
#endif
				active = true;
			}

			// *******************************************************************************************
			virtual ~output_stdout()
			{}

			// *******************************************************************************************
			virtual bool close() final
			{
				fflush(stdout);

				return true;
			}

			// *******************************************************************************************
			virtual bool put(const uint8_t c) final
			{
				if (std::putc((int)c, stdout) == EOF)
					return false;

				++f_size;

				return true;
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) final
			{
				auto written = refresh::io::fwrite(ptr, 1, n_bytes, stdout);
			
				f_size += written;

				return written;
			}

			// *******************************************************************************************
			virtual size_t write_uint(const uint64_t val, const size_t n_bytes) final
			{
				uint64_t x = val;

				for (size_t i = 0; i < n_bytes; ++i)
				{
					if (std::putc((int)(x & 0xffull), stdout) == EOF)
						return 0;

					x >>= 8;
					++f_size;
				}

				return n_bytes;
			}
		};

		// *******************************************************************************************
		//
		// *******************************************************************************************
		class output_stdout_buffered : public output_common
		{
			size_t buffer_size;
			std::vector<uint8_t> buffer;

			// *******************************************************************************************
			bool flush_buffer()
			{
				if (refresh::io::fwrite(buffer.data(), 1, buffer.size(), stdout) != buffer.size())
					return false;

				buffer.clear();

				return true;
			}

		public:
			// *******************************************************************************************
			output_stdout_buffered(const size_t buffer_size = 1 << 20) :
				output_common(),
				buffer_size(buffer_size)
			{
				buffer.reserve(buffer_size);

#ifdef _WIN32
				int res = _setmode(_fileno(stdout), _O_BINARY);
				active = res != -1;
#endif
				active = true;
			}

			// *******************************************************************************************
			virtual ~output_stdout_buffered()
			{}

			// *******************************************************************************************
			virtual bool close() final
			{
				if (!buffer.empty())
					if (!flush_buffer())
						return false;
				
				fflush(stdout);

				return true;
			}

			// *******************************************************************************************
			virtual bool put(const uint8_t c) final
			{
				buffer.emplace_back(c);

				++f_size;

				if (buffer.size() == buffer_size)
					if (!flush_buffer())
						return true;

				return true;
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) final
			{
				if (buffer.size() + n_bytes >= buffer_size)
				{
					if (!flush_buffer())
						return 0;
					if (refresh::io::fwrite(ptr, 1, n_bytes, stdout) != n_bytes)
						return 0;
				}
				else
					buffer.insert(buffer.end(), ptr, ptr + n_bytes);

				f_size += n_bytes;

				return n_bytes;
			}

			// *******************************************************************************************
			virtual size_t write_uint(const uint64_t val, const size_t n_bytes) final
			{
				uint64_t x = val;

				for (size_t i = 0; i < n_bytes; ++i)
				{
					if (!put((uint8_t)(x & 0xffull)))
						return 0;
					x >>= 8;
				}

				return n_bytes;
			}
		};
	}
}