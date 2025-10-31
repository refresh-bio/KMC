#pragma once

#include "output_common.h"
#include "io_utils.h"

#include <string>
#include <cstring>
#include <future>
#include <vector>
#include <list>
#include <stack>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>

#ifdef _WIN32
#include <windows.h>
#include <malloc.h> 
#else // Linux, MacOS
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <cstdlib>
#endif

namespace refresh
{
	namespace io
	{
#ifdef _WIN32
		using FileHandleType = HANDLE;
		const FileHandleType INVALID_HANDLE = INVALID_HANDLE_VALUE;
#else // Linux, MacOS
		using FileHandleType = int;
		const FileHandleType INVALID_HANDLE = -1;
#endif

		// *******************************************************************************************
		// 
		// *******************************************************************************************
		class output_file_low_level : public output_common
		{
		protected:
			std::string file_name;
			FileHandleType file_handle = INVALID_HANDLE;
			size_t file_size_on_disk{};
			size_t prealocated_size{};
			uint8_t* buffer;
			size_t buffer_offset{};
			const size_t ALIGNMENT = 4096;
			size_t buffer_size;

			// *******************************************************************************************
			void allocate_buffer()
			{
#ifdef _WIN32
				// Windows
				buffer = static_cast<uint8_t*>(_aligned_malloc(buffer_size, ALIGNMENT));
				if (buffer == nullptr) 
					throw std::bad_alloc();
#else
				// Linux, macOS
				if (posix_memalign((void **) & buffer, ALIGNMENT, buffer_size) != 0)
					throw std::bad_alloc();
#endif
			}

			// *******************************************************************************************
			void free_buffer()
			{
				if (buffer == nullptr) 
					return;

#ifdef _WIN32
				_aligned_free(buffer);
#else
				free(buffer);
#endif
			}

			// *******************************************************************************************
			bool open()
			{
#ifdef _WIN32
				// Windows
				file_handle = CreateFileA(
					file_name.c_str(),
					GENERIC_WRITE,
					0, // No sharing
					NULL,
					CREATE_ALWAYS,
					FILE_ATTRIBUTE_NORMAL | FILE_FLAG_NO_BUFFERING, // No buffering
					NULL
				);

				if (file_handle == INVALID_HANDLE) 
					return false;
#else 
				// Linux, macOS
#ifdef __APPLE__
				// macOS - F_NOCACHE must be set after open
				file_handle = ::open(file_name.c_str(), O_WRONLY | O_CREAT, 0644);
/*				if (file_handle != INVALID_HANDLE) {
					if (fcntl(file_handle, F_NOCACHE, 1) == -1) {
//						Warning: Impossible to set F_NOCACHE
					}
				}*/
#else 
				// Linux - O_DIRECT
//				file_handle = ::open(file_name.c_str(), O_WRONLY | O_CREAT | O_DIRECT, 0644);
				file_handle = ::open(file_name.c_str(), O_WRONLY | O_CREAT, 0644);
#endif
				if (file_handle == INVALID_HANDLE) {
					perror("open");
					return false;
				}
#endif

				buffer_offset = 0;

				return true;
			}

			// *******************************************************************************************
			bool write_block_to_file(const uint8_t* ptr, const size_t n_bytes)
			{
#ifdef _WIN32
				// Windows
				DWORD bytes_written;
				BOOL success = WriteFile(
					file_handle,
					ptr,
					(DWORD)n_bytes,
					&bytes_written,
					NULL
				);
				
				if (!success || bytes_written != n_bytes)
					return false;

#else 
				// Linux, macOS
#ifdef __linux__
/*				if (file_size_on_disk + n_bytes > prealocated_size)
				{
					prealocated_size = file_size_on_disk + n_bytes + 16 * buffer_size;
					if (fallocate(file_handle, 0, 0, prealocated_size) != 0)
						if (ftruncate(file_handle, prealocated_size) != 0) 
							perror("ftruncate");
				}*/

				if (posix_fallocate(file_handle, 0, file_size_on_disk + n_bytes) != 0)
					if (ftruncate(file_handle, file_size_on_disk + n_bytes) != 0) 
						perror("ftruncate");
#else
				if (ftruncate(file_handle, file_size_on_disk + n_bytes) != 0)
					perror("ftruncate");
#endif			

				ssize_t bytes_written = pwrite(file_handle, ptr, n_bytes, file_size_on_disk);

				if(bytes_written == -1 || (size_t)bytes_written != n_bytes)
				{
					perror("pwrite");
					return false;
				}
#endif

//				file_size_on_disk += bytes_written;

				return true;
			}

		public:
			// *******************************************************************************************
			output_file_low_level(const std::string& file_name, size_t buffer_size) :
				output_common(),
				file_name(file_name),
				buffer_size(buffer_size)
			{
				if(buffer_size % ALIGNMENT != 0)
					this->buffer_size = ((buffer_size / ALIGNMENT) + 1) * ALIGNMENT;

				allocate_buffer();
				open();
				active = file_handle != INVALID_HANDLE;
			}

			// *******************************************************************************************
			~output_file_low_level()
			{
				close();
				free_buffer();
			}

			// *******************************************************************************************
			void prealocate_file_space(size_t TOTAL_SIZE)
			{
#ifdef _WIN32
				// Windows
				LONG high_size = (LONG)(TOTAL_SIZE >> 32);
				DWORD low_size = (DWORD)TOTAL_SIZE;
				if (SetFilePointer(file_handle, low_size, &high_size, FILE_BEGIN) != INVALID_SET_FILE_POINTER && GetLastError() == NO_ERROR) 
					SetEndOfFile(file_handle);
#else
#ifdef __linux__
				// Linux, macOS
				if (posix_fallocate(file_handle, 0, TOTAL_SIZE) != 0)
					if (ftruncate(file_handle, TOTAL_SIZE) != 0)
						perror("ftruncate");
#else
				// macOS/Inne Unixy: używamy ftruncate do ustawienia rozmiaru
				if (ftruncate(file_handle, TOTAL_SIZE) != 0) {
					perror("ftruncate");
				}
#endif			
#endif			
			}

			// *******************************************************************************************
			virtual bool close() 
			{
				if (!active)
					return false;

				size_t bytes_to_flush = (buffer_offset + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;

				if (bytes_to_flush != 0)
				{
					write_block_to_file(buffer, bytes_to_flush);
					file_size_on_disk += buffer_offset;
				}

#ifdef _WIN32
				CloseHandle(file_handle);
#else
				::close(file_handle);
#endif

				// Truncate file to actual size
				if (bytes_to_flush != 0)
				{
#ifdef _WIN32
					// Windows
					HANDLE hFile_trim = CreateFileA(
						file_name.c_str(),
						GENERIC_WRITE,
						0,
						NULL,
						OPEN_EXISTING, 
						FILE_ATTRIBUTE_NORMAL, // Buffered mode
						NULL
					);

					if (hFile_trim != INVALID_HANDLE_VALUE) {
						LONG high_size = (LONG)(file_size_on_disk >> 32);
						DWORD low_size = (DWORD)file_size_on_disk;

						if (SetFilePointer(hFile_trim, low_size, &high_size, FILE_BEGIN) != INVALID_SET_FILE_POINTER || GetLastError() == NO_ERROR)
							SetEndOfFile(hFile_trim);
						CloseHandle(hFile_trim);
					}
#else 
					// Linux, macOS
					int fd_trim = ::open(file_name.c_str(), O_WRONLY);
					if (fd_trim != INVALID_HANDLE) 
					{
						if(ftruncate(fd_trim, file_size_on_disk) != 0)
							perror("ftruncate");
						::close(fd_trim);
					}
#endif
				}

				active = false;

				return true;
			}

			// *******************************************************************************************
			virtual bool put(const uint8_t c) 
			{
				buffer[buffer_offset++] = c;

				if (buffer_offset == buffer_size)
				{
					if(!write_block_to_file(buffer, buffer_size))
						return false;

					file_size_on_disk += buffer_size;
					buffer_offset = 0;
				}

				return true;
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) 
			{
				uint8_t* local_ptr = const_cast<uint8_t*>(ptr);
				size_t bytes_to_write = 0;

				if (buffer_offset == 0 && uintptr_t(local_ptr) % ALIGNMENT == 0)
				{
					size_t full_blocks = n_bytes / ALIGNMENT;

					bytes_to_write = full_blocks * ALIGNMENT;
					if (bytes_to_write != 0 && !write_block_to_file(local_ptr, bytes_to_write))
						return 0;

					file_size_on_disk += bytes_to_write;
					local_ptr += bytes_to_write;
					bytes_to_write = n_bytes - bytes_to_write;
				}
				else
					bytes_to_write = n_bytes;

				while(bytes_to_write != 0)
				{
					size_t space_in_buffer = buffer_size - buffer_offset;
					size_t to_copy = (bytes_to_write < space_in_buffer) ? bytes_to_write : space_in_buffer;
					std::memcpy(buffer + buffer_offset, local_ptr, to_copy);
					buffer_offset += to_copy;
					local_ptr += to_copy;
					bytes_to_write -= to_copy;

					if (buffer_offset == buffer_size)
					{
						if (!write_block_to_file(buffer, buffer_size))
							return 0;
						buffer_offset = 0;
						file_size_on_disk += buffer_size;
					}
				}

				return n_bytes;
			}

			// *******************************************************************************************
			virtual size_t write_uint(const uint64_t val, const size_t n_bytes) 
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

		// *******************************************************************************************
		// 
		// *******************************************************************************************
		class output_file_low_level_threaded : public output_common
		{
			class simple_queue
			{
			private:
				size_t capacity;
				std::list<std::pair<uint8_t*, size_t>> data;

				std::mutex mtx;
				std::condition_variable cv_not_empty;
				std::condition_variable cv_not_full;
				bool stopped{ false };

			public:
				simple_queue(size_t capacity) : capacity(capacity)
				{
				}

				bool push(std::pair<uint8_t*, size_t> item)
				{
					std::unique_lock<std::mutex> lck(mtx);
					cv_not_full.wait(lck, [this]() { return data.size() < capacity; });

					data.push_back(item);

					cv_not_empty.notify_one();

					return true;
				}

				bool pop(std::pair<uint8_t*, size_t>& item)
				{
					std::unique_lock<std::mutex> lck(mtx);
					cv_not_empty.wait(lck, [this]() { return !data.empty() || stopped; });

					if(data.empty() && stopped)
					{
						item = std::make_pair<uint8_t*, size_t>(nullptr, 0);
						return false;
					}

					item = data.front();
					data.pop_front();

					cv_not_full.notify_one();
					
					return true;
				}

				void stop()
				{
					std::unique_lock<std::mutex> lck(mtx);
					stopped = true;
					cv_not_empty.notify_all();
				}
			};

		protected:
			std::string file_name;
			FileHandleType file_handle = INVALID_HANDLE;
			size_t file_size_on_disk{};
			size_t prealocated_size{};

			uint8_t* buffer;
			size_t buffer_offset{};
			const size_t ALIGNMENT = 4096;
			size_t buffer_size;
			size_t queue_capacity;

			simple_queue write_queue;
			std::stack<uint8_t*> free_buffers;
			std::mutex free_buffers_mtx;

			std::thread write_thread;
			std::atomic<bool> writing_ok{ true };

			// *******************************************************************************************
			uint8_t *allocate_buffer()
			{
				uint8_t* ptr;

#ifdef _WIN32
				// Windows
				ptr = static_cast<uint8_t*>(_aligned_malloc(buffer_size, ALIGNMENT));
				if (ptr == nullptr)
					throw std::bad_alloc();
#else
				// Linux, macOS
				if (posix_memalign((void**)&ptr, ALIGNMENT, buffer_size) != 0)
					throw std::bad_alloc();
#endif

				return ptr;
			}

			// *******************************************************************************************
			void free_buffer(uint8_t *ptr)
			{
				if (ptr == nullptr)
					return;

#ifdef _WIN32
				_aligned_free(ptr);
#else
				free(ptr);
#endif
			}

			// *******************************************************************************************
			bool open()
			{
#ifdef _WIN32
				// Windows
				file_handle = CreateFileA(
					file_name.c_str(),
					GENERIC_WRITE,
					0, // No sharing
					NULL,
					CREATE_ALWAYS,
					FILE_ATTRIBUTE_NORMAL | FILE_FLAG_NO_BUFFERING, // No buffering
					NULL
				);

				if (file_handle == INVALID_HANDLE)
					return false;
#else 
				// Linux, macOS
#ifdef __APPLE__
				// macOS - F_NOCACHE must be set after open
				file_handle = ::open(file_name.c_str(), O_WRONLY | O_CREAT, 0644);
				/*				if (file_handle != INVALID_HANDLE) {
									if (fcntl(file_handle, F_NOCACHE, 1) == -1) {
				//						Warning: Impossible to set F_NOCACHE
									}
								}*/
#else 
				// Linux - O_DIRECT
//				file_handle = ::open(file_name.c_str(), O_WRONLY | O_CREAT | O_DIRECT, 0644);
				file_handle = ::open(file_name.c_str(), O_WRONLY | O_CREAT, 0644);
#endif
				if (file_handle == INVALID_HANDLE) {
					perror("open");
					return false;
				}
#endif

				buffer_offset = 0;

				return true;
			}

			// *******************************************************************************************
			bool write_block_to_file(const uint8_t* ptr, const size_t n_bytes)
			{
#ifdef _WIN32
				// Windows
				DWORD bytes_written;
				BOOL success = WriteFile(
					file_handle,
					ptr,
					(DWORD)n_bytes,
					&bytes_written,
					NULL
				);

				if (!success || bytes_written != n_bytes)
					return false;

#else 
				// Linux, macOS
#ifdef __linux__
/*				if (file_size_on_disk + n_bytes > prealocated_size)
				{
					prealocated_size = file_size_on_disk + n_bytes + 16 * buffer_size;
					if (fallocate(file_handle, 0, 0, prealocated_size) != 0)
						if (ftruncate(file_handle, prealocated_size) != 0)
							perror("ftruncate");
				}*/

				if (posix_fallocate(file_handle, 0, file_size_on_disk + n_bytes) != 0)
					if (ftruncate(file_handle, file_size_on_disk + n_bytes) != 0)
						perror("ftruncate");
#else
				if (ftruncate(file_handle, file_size_on_disk + n_bytes) != 0)
					perror("ftruncate");
#endif			

				ssize_t bytes_written = pwrite(file_handle, ptr, n_bytes, file_size_on_disk);

				if (bytes_written == -1 || (size_t)bytes_written != n_bytes)
				{
					perror("pwrite");
					return false;
				}
#endif

				//				file_size_on_disk += bytes_written;

				return true;
			}

			// *******************************************************************************************
			void put_buffer_to_queue()
			{
				if (writing_ok.load() == false)
				{
					buffer_offset = 0;
					return;
				}

				free_buffers_mtx.lock();
				uint8_t* ptr;
				if (free_buffers.empty())
					ptr = allocate_buffer();
				else
				{
					ptr = free_buffers.top();
					free_buffers.pop();
				}

				free_buffers_mtx.unlock();

				std::swap(buffer, ptr);
				write_queue.push(std::make_pair(ptr, buffer_offset));
				buffer_offset = 0;
			}

			// *******************************************************************************************
			size_t align_size(size_t size)
			{
				if (size % ALIGNMENT == 0)
					return size;
				else
					return ((size / ALIGNMENT) + 1) * ALIGNMENT;
			}

		public:
			// *******************************************************************************************
			output_file_low_level_threaded(const std::string& file_name, size_t buffer_size, size_t queue_capacity) :
				output_common(),
				file_name(file_name),
				buffer_size(buffer_size),
				queue_capacity(queue_capacity),
				write_queue(queue_capacity)

			{
				if (buffer_size % ALIGNMENT != 0)
					this->buffer_size = align_size(buffer_size);

				buffer = allocate_buffer();

//				for(size_t i = 0; i < queue_capacity + 2; ++i)
//					free_buffers.push(allocate_buffer());

				open();
				active = file_handle != INVALID_HANDLE;

				if(active == false)
					return;

				write_thread = std::thread([this]() {
					std::pair<uint8_t*, size_t> block;

					bool last_block_not_aligned = false;

					while (write_queue.pop(block))
					{
						if(last_block_not_aligned)
						{
							// Previous block was not aligned, so the file is broken
							writing_ok = false;
							continue;
						}

						size_t aligned_size = align_size(block.second);
						last_block_not_aligned = (aligned_size != block.second);

						if (!write_block_to_file(block.first, aligned_size))
							writing_ok = false;

						file_size_on_disk += block.second;

						free_buffers_mtx.lock();
						free_buffers.push(block.first);
						free_buffers_mtx.unlock();
					}
				});
			}

			// *******************************************************************************************
			~output_file_low_level_threaded()
			{
				try
				{
					(void) close();
				}
				catch (...) {}

				free_buffer(buffer);

				while (!free_buffers.empty())
				{
					uint8_t *ptr = free_buffers.top();
					free_buffers.pop();
					free_buffer(ptr);
				}
			}

			// *******************************************************************************************
			virtual bool close()
			{
				if (!active)
					return false;

				bool need_to_flush = (buffer_offset != 0);

				if (need_to_flush)
					put_buffer_to_queue();

				write_queue.stop();
				write_thread.join();

#ifdef _WIN32
				CloseHandle(file_handle);
#else
				::close(file_handle);
#endif

				// Truncate file to actual size
				if (writing_ok.load() && need_to_flush)
				{
#ifdef _WIN32
					// Windows
					HANDLE hFile_trim = CreateFileA(
						file_name.c_str(),
						GENERIC_WRITE,
						0,
						NULL,
						OPEN_EXISTING,
						FILE_ATTRIBUTE_NORMAL, // Buffered mode
						NULL
					);

					if (hFile_trim != INVALID_HANDLE_VALUE) {
						LONG high_size = (LONG)(file_size_on_disk >> 32);
						DWORD low_size = (DWORD)file_size_on_disk;

						if (SetFilePointer(hFile_trim, low_size, &high_size, FILE_BEGIN) != INVALID_SET_FILE_POINTER || GetLastError() == NO_ERROR)
							SetEndOfFile(hFile_trim);
						CloseHandle(hFile_trim);
					}
#else 
					// Linux, macOS
					int fd_trim = ::open(file_name.c_str(), O_WRONLY);
					if (fd_trim != INVALID_HANDLE)
					{
						if (ftruncate(fd_trim, file_size_on_disk) != 0)
							perror("ftruncate");
						::close(fd_trim);
					}
#endif
				}

				active = false;

				return writing_ok.load();			// To tell if writing was successful
			}

			// *******************************************************************************************
			virtual bool put(const uint8_t c)
			{
				if (writing_ok.load() == false)
				{
					write_queue.stop();
					return false;
				}

				buffer[buffer_offset++] = c;

				if (buffer_offset == buffer_size)
					put_buffer_to_queue();

				return writing_ok.load();
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes)
			{
				uint8_t* local_ptr = const_cast<uint8_t*>(ptr);
				size_t bytes_to_write = n_bytes;

				while (bytes_to_write != 0 && writing_ok.load())
				{
					size_t space_in_buffer = buffer_size - buffer_offset;
					size_t to_copy = (bytes_to_write < space_in_buffer) ? bytes_to_write : space_in_buffer;
					std::memcpy(buffer + buffer_offset, local_ptr, to_copy);
					buffer_offset += to_copy;
					local_ptr += to_copy;
					bytes_to_write -= to_copy;

					if (buffer_offset == buffer_size)
						put_buffer_to_queue();
				}

				if(!writing_ok.load())
					write_queue.stop();

				return writing_ok ? n_bytes : 0;
			}

			// *******************************************************************************************
			virtual size_t write_uint(const uint64_t val, const size_t n_bytes)
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


		// *******************************************************************************************
		// 
		// *******************************************************************************************
#if 0
		class output_file_memory_mapped_reopen : public output_file_memory_mapped
		{
			std::string file_name;
			bool in_transaction{};

		public:
			// *******************************************************************************************
			output_file_memory_mapped_reopen(const std::string& file_name, const size_t buffer_size) :
				output_file_memory_mapped(file_name, buffer_size),
				file_name(file_name)
			{
				// After creating the file it is closed and we are out-of-transaction
/*				if (f)
				{
					fclose(f);
					f = nullptr;
				}*/
			}

			// *******************************************************************************************
			~output_file_memory_mapped_reopen()
			{
				if (active && !buffer.empty())
					close();
			}

			// *******************************************************************************************
			virtual bool close() final
			{
//				if (!f)
				{
/*					f = fopen(file_name.c_str(), "ab");
					if (!f)
						return false;
						*/
					if (!flush_buffer())
						return false;

/*					if (fclose(f) != 0)
						return false;
					f = nullptr;*/
				}

				active = false;

				return true;
			}

			// *******************************************************************************************
			virtual bool start_transaction()
			{
				if (in_transaction)
					return false;

				// In out-of-transaction state file must be closed
//				assert(f == nullptr);

				// Delayed file open - no need to open file until need to store buffer contents

				in_transaction = true;

				return true;
			};

			// *******************************************************************************************
			virtual bool stop_transaction()
			{
				if (!in_transaction)
					return false;

				in_transaction = false;
				
/*				if (f)					// file can be not opened due to delayed file open in start_transaction
				{
					fclose(f);
					f = nullptr;
				}
				*/
				return true;
			};

			// *******************************************************************************************
			virtual bool put(const uint8_t c) final
			{
				buffer.emplace_back(c);

				++f_size;

				if (buffer.size() == buffer_size)
				{
/*					if (!f && in_transaction)
					{
						f = fopen(file_name.c_str(), "ab");
						if (!f)
							return false;
					}
					*/
					if (!flush_buffer())
						return false;
				}

				return true;
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) final
			{
				if (buffer.size() + n_bytes >= buffer_size)
				{
/*					if (!f && in_transaction)
					{
						f = fopen(file_name.c_str(), "ab");
						if (!f)
							return false;
					}*/

					if (!flush_buffer())
						return false;
					if (!flush_block(ptr, n_bytes))
						return false;
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
#endif
	}
}