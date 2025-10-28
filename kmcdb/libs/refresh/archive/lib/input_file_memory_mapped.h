#pragma once

#include <cstring>
#include <filesystem>
#include <vector>
#include <string>
#include <algorithm>
#include <sys/types.h>

#ifdef _WIN32
// Windows
#include <windows.h>
#include <io.h>
#else
// Linux, MacOS
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif

#include "input_common.h"

namespace refresh
{
	namespace io
	{
		class input_file_memory_mapped : public input_common
		{
		protected:
#ifdef _WIN32
			using FileHandleType = HANDLE;
			const FileHandleType INVALID_HANDLE = INVALID_HANDLE_VALUE;
#else // Linux, MacOS
			using FileHandleType = int;
			const FileHandleType INVALID_HANDLE = -1;
#endif

			FileHandleType file_handle = INVALID_HANDLE;
			FileHandleType map_handle = INVALID_HANDLE;

			uint8_t* segment = nullptr;

			size_t segment_size{};
			size_t in_segment_pos{};
			size_t segment_filled{};
			size_t segment_start{};

			// *******************************************************************************************
			size_t determine_segment_size(const size_t preferred_size)
			{
#ifdef _WIN32
				SYSTEM_INFO sys_info;
				GetSystemInfo(&sys_info);
				size_t allocation_granularity = sys_info.dwAllocationGranularity;
#else
				auto page_size = sysconf(_SC_PAGESIZE);
				size_t allocation_granularity = page_size != -1 ? static_cast<size_t>(page_size) : 4096;
#endif

				if (preferred_size % allocation_granularity != 0)
					return ((preferred_size / allocation_granularity) + 1) * allocation_granularity;
				else
					return preferred_size;
			}

			// *******************************************************************************************
			bool map_file_section(const size_t offset, const size_t size)
			{
				segment_start = offset / segment_size * segment_size;
				segment_filled = std::min<size_t>(size, f_size - segment_start);

#ifdef _WIN32
				if(segment != nullptr)
					UnmapViewOfFile(segment);

				// Windows: MapViewOfFile
				void *addr = MapViewOfFile(
					map_handle, FILE_MAP_READ,
					(DWORD)(segment_start >> 32), (DWORD)(segment_start & 0xFFFFFFFF), // High/Low part of offset
					segment_filled);

				if (addr == NULL)
					return false;
#else
				if(segment != nullptr)
					munmap(segment, segment_filled);

				// POSIX: mmap 
				void *addr = mmap(NULL, segment_filled, PROT_READ, MAP_SHARED, file_handle, segment_start);

				if (addr == MAP_FAILED)
					return false;

				// Optional optimization 
				madvise(addr, segment_filled, MADV_SEQUENTIAL);
#endif

				segment = static_cast<uint8_t*>(addr);

				in_segment_pos = offset - segment_start;

				return true;
			}

			// *******************************************************************************************
			bool open(const std::string &file_name)
			{
#ifdef _WIN32
				file_handle = CreateFileA(
					file_name.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
				/*				if (hFile == INVALID_HANDLE_VALUE) {
				//					std::cerr << "CreateFile error: " << GetLastError() << std::endl;
									active = false;
									return;
								}*/

				map_handle = CreateFileMapping(
					file_handle, NULL, PAGE_READONLY, 0, 0, NULL);
				if (map_handle == NULL) {
					//					std::cerr << "CreateFileMapping error: " << GetLastError() << std::endl;
					CloseHandle(file_handle);
					//					active = false;
					//					return;
					file_handle = INVALID_HANDLE;
				}
#else
				// POSIX
				file_handle = ::open(file_name.c_str(), O_RDONLY);
				/*				if (fd < 0) {
				//					std::cerr << "open error: " << strerror(errno) << std::endl;
									active = false;
									return;
								}*/
#endif

				active = file_handle != INVALID_HANDLE;

				return active;
			}

		public:
			// *******************************************************************************************
			input_file_memory_mapped(const std::string& file_name, const size_t _segment_size = 32 << 20) :
				input_common()
			{
				std::filesystem::path fp(file_name);

				if (!std::filesystem::exists(fp))
					return;

				segment_size = determine_segment_size(_segment_size);

				f_size = std::filesystem::file_size(fp);
				f_pos = 0;

				if (f_size < segment_size)
					segment_size = f_size;

				(void) open(file_name);
				(void)map_file_section(0, segment_size);
			}

			// *******************************************************************************************
			virtual ~input_file_memory_mapped()
			{
				close();
			}

			// *******************************************************************************************
			virtual bool close()
			{
				if (file_handle != INVALID_HANDLE)
				{
#ifdef _WIN32
					if (segment != nullptr)
						UnmapViewOfFile(segment);

					CloseHandle(map_handle);
					CloseHandle(file_handle);
					map_handle = INVALID_HANDLE;
					file_handle = INVALID_HANDLE;
#else
					if (segment != nullptr)
						munmap(segment, segment_filled);

					::close(file_handle);
					file_handle = INVALID_HANDLE;
#endif				
				}

				segment = nullptr;
				in_segment_pos = 0;
				segment_filled = 0;

				active = false;

				return true;
			}

			// *******************************************************************************************
			virtual bool get(uint8_t& x)
			{
				if (in_segment_pos < segment_filled)
				{
					x = segment[in_segment_pos++];
					++f_pos;
					return true;
				}
					
				if (segment_filled == 0)
					return false;

				if(!map_file_section(segment_start + segment_size, segment_size))
					return false;

				if (in_segment_pos >= segment_filled)
					return false;

				x = segment[in_segment_pos++];
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
				if (segment_start + in_segment_pos + size > f_size)
					size = f_size - (segment_start + in_segment_pos);

				size_t to_read = size;
				uint8_t *ptr0 = ptr;

				while (in_segment_pos + to_read > segment_filled)
				{
					auto chunk_size = segment_filled - in_segment_pos;
					memcpy(ptr, segment + in_segment_pos, chunk_size);
					ptr += chunk_size;
					to_read -= chunk_size;

					if(!map_file_section(segment_start + segment_size, segment_size) || segment_filled == 0)
					{
						f_pos += ptr - ptr0;
						return ptr - ptr0;
					}

					in_segment_pos = 0;
				}

				memcpy(ptr, segment + in_segment_pos, to_read);
				in_segment_pos += to_read;
				ptr += to_read;

				f_pos += ptr - ptr0;

				return ptr - ptr0;
			}
			
			// *******************************************************************************************
			virtual bool eof() const
			{
//				return before_buffer_bytes + buffer_pos >= f_pos;			// !!! Was it a bug ?
				return segment_start + in_segment_pos >= f_size;
			}
			
			// *******************************************************************************************
			virtual bool seek(const size_t requested_pos)
			{
				if (requested_pos > f_size)
					return false;

				if (requested_pos >= segment_start && requested_pos < segment_start + segment_filled)
					in_segment_pos = requested_pos - segment_start;
				else
					map_file_section(requested_pos, segment_size);

				f_pos = requested_pos;

				return true;
			}
		};

		// *******************************************************************************************
		//
		// *******************************************************************************************
		class input_file_memory_mapped_reopen : public input_file_memory_mapped
		{
			std::string file_name;
			bool in_transaction;

		public:
			// *******************************************************************************************
			input_file_memory_mapped_reopen(const std::string& file_name, const size_t _buffer_size = 32 << 20) :
				input_file_memory_mapped(file_name, _buffer_size),
				file_name(file_name),
				in_transaction(false)
			{
				if (active)		// In reopen-mode file is closed after the openning
				{
					close();
					active = true;		// close() sets active=false but we want to stay active in reopen mode
				}
			}

			// *******************************************************************************************
			virtual ~input_file_memory_mapped_reopen()
			{}

			// *******************************************************************************************
			virtual bool start_transaction()
			{
				if (in_transaction)
					return false;

				assert(file_handle == INVALID_HANDLE);

				if(!open(file_name))
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

				assert(file_handle != INVALID_HANDLE);

				in_transaction = false;
				if (!close() != 0)
					return false;

				active = true;		// close() sets active=false but we want to stay active in reopen mode
				
				return true;
			}
		};
	}
}