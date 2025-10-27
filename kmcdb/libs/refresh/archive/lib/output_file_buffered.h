#pragma once

#include "output_common.h"
#include "io_utils.h"

namespace refresh
{
	namespace io
	{
		// *******************************************************************************************
		// 
		// *******************************************************************************************
		class output_file_buffered : public output_common
		{
		protected:
			FILE* f{};
			std::vector<uint8_t> buffer;
			const size_t buffer_size;

			// *******************************************************************************************
			bool flush_buffer()
			{
				if (refresh::io::fwrite(buffer.data(), 1, buffer.size(), f) != buffer.size())
					return false;

				buffer.clear();
				return true;
			}

		public:
			// *******************************************************************************************
			output_file_buffered(const std::string& file_name, size_t buffer_size) :
				output_common(),
				buffer_size(buffer_size)
			{
				f = fopen(file_name.c_str(), "wb");

				buffer.reserve(buffer_size);

				active = f != nullptr;
			}

			// *******************************************************************************************
			~output_file_buffered()
			{
				close();
			}

			// *******************************************************************************************
			virtual bool close() 
			{
				if (f)
				{
					if (!buffer.empty())
						if (!flush_buffer())
							return false;

					if (fclose(f) != 0)
						return false;
					f = nullptr;
				}

				active = false;

				return true;
			}

			// *******************************************************************************************
			virtual bool put(const uint8_t c) 
			{
				buffer.emplace_back(c);

				++f_size;

				if (buffer.size() == buffer_size)
					if (!flush_buffer())
						return false;

				return true;
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) 
			{
				if (buffer.size() + n_bytes >= buffer_size)
				{
					if (!flush_buffer())
						return 0;
					if (refresh::io::fwrite(ptr, 1, n_bytes, f) != n_bytes)
						return 0;
				}
				else
					buffer.insert(buffer.end(), ptr, ptr + n_bytes);

				f_size += n_bytes;

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
		class output_file_buffered_reopen : public output_file_buffered
		{
			std::string file_name;
			bool in_transaction{};

		public:
			// *******************************************************************************************
			output_file_buffered_reopen(const std::string& file_name, const size_t buffer_size) :
				output_file_buffered(file_name, buffer_size),
				file_name(file_name)
			{
				// After creating the file it is closed and we are out-of-transaction
				if (f)
				{
					fclose(f);
					f = nullptr;
				}
			}

			// *******************************************************************************************
			~output_file_buffered_reopen()
			{
				if (active && !buffer.empty())
					close();
			}

			// *******************************************************************************************
			virtual bool close() final
			{
				if (!f)
				{
					f = fopen(file_name.c_str(), "ab");
					if (!f)
						return false;

					if (!flush_buffer())
						return false;

					if (fclose(f) != 0)
						return false;
					f = nullptr;
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
				assert(f == nullptr);

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
				
				if (f)					// file can be not opened due to delayed file open in start_transaction
				{
					fclose(f);
					f = nullptr;
				}

				return true;
			};

			// *******************************************************************************************
			virtual bool put(const uint8_t c) final
			{
				buffer.emplace_back(c);

				++f_size;

				if (buffer.size() == buffer_size)
				{
					if (!f && in_transaction)
					{
						f = fopen(file_name.c_str(), "ab");
						if (!f)
							return false;
					}

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
					if (!f && in_transaction)
					{
						f = fopen(file_name.c_str(), "ab");
						if (!f)
							return false;
					}

					if (!flush_buffer())
						return false;
					if (refresh::io::fwrite(ptr, 1, n_bytes, f) != n_bytes)
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
	}
}