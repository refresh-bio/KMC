#ifndef _ARCHIVE_HELPER_H
#define _ARCHIVE_HELPER_H

#include <cstdio>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace refresh
{
	namespace io
	{
		// *******************************************************************************************
		// Write in 2GB parts - due to MSVC bug when writing 4GB+ data at one call
		inline size_t fwrite(const void* buffer, size_t size, size_t count, FILE* stream)
		{
			const size_t part_size = 2ull << 30;

			size_t tot_written = 0;
			size_t to_write = size * count;
			const char* ptr = static_cast<const char*>(buffer);

			while (true)
			{
				if (to_write <= part_size)
				{
					tot_written += std::fwrite(ptr, 1, to_write, stream);
					break;
				}

				tot_written += std::fwrite(ptr, 1, part_size, stream);
				to_write -= part_size;
				ptr += part_size;
			}

			return tot_written;
		}

		// *******************************************************************************************
		inline size_t fread(void* buffer, size_t size, size_t count, FILE* stream)
		{
			const size_t part_size = 2ull << 30;

			size_t tot_readed = 0;
			size_t to_read = size * count;
			char* ptr = static_cast<char*>(buffer);

			while (true)
			{
				if (to_read <= part_size)
				{
					tot_readed += std::fread(ptr, 1, to_read, stream);
					break;
				}

				tot_readed += std::fread(ptr, 1, part_size, stream);
				to_read -= part_size;
				ptr += part_size;
			}

			return tot_readed;
		}

		// *******************************************************************************************
		inline int fseek(FILE* stream, long long offset, int origin)
		{
#ifdef _WIN32
			return _fseeki64(stream, offset, origin);
#else
			static_assert(sizeof(long) >= 8); //if it fails we need to figure something else than std::fseek
			return std::fseek(stream, offset, origin);
#endif
		}

		// *******************************************************************************************
		inline long long ftell(std::FILE* stream)
		{
#ifdef _WIN32
			return _ftelli64(stream);
#else
			return std::ftell(stream);
#endif
		}
	}
}

#endif