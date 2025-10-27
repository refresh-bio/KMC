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
		class output_file_unbuffered : public output_common
		{
		protected:
			FILE* f{};

		public:
			// *******************************************************************************************
			output_file_unbuffered(const std::string &file_name) :
				output_common()
			{
				f = fopen(file_name.c_str(), "wb");

				active = f != nullptr;
			}

			// *******************************************************************************************
			~output_file_unbuffered()
			{
				if (f)
					fclose(f);
			}

			// *******************************************************************************************
			virtual bool close()
			{
				if (f)
				{
					if (fclose(f) != 0)
						return false;
					f = nullptr;
				}

				active = false;

				return true;
			}

			// *******************************************************************************************
			virtual bool put(const uint8_t c) final
			{
				if (std::putc(c, f) == EOF)
					return false;

				++f_size;

				return true;
			}

			// *******************************************************************************************
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) final
			{
				auto written = refresh::io::fwrite(ptr, 1, n_bytes, f);

				f_size += written;

				return written;
			}

			// *******************************************************************************************
			virtual size_t write_uint(const uint64_t val, const size_t n_bytes) final
			{
				uint64_t x = val;

				for (size_t i = 0; i < n_bytes; ++i)
				{
					int c = (int) (x & 0xffull);

					if (std::putc(c, f) == EOF)
						return 0;

					x >>= 8;
				}

				return n_bytes;
			}
		};

		// *******************************************************************************************
		// 
		// *******************************************************************************************
		class output_file_unbuffered_reopen : public output_file_unbuffered
		{
			std::string file_name;
			bool in_transaction{};

		public:
			// *******************************************************************************************
			output_file_unbuffered_reopen(const std::string& file_name) :
				output_file_unbuffered(file_name),
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
			~output_file_unbuffered_reopen() = default;

			// *******************************************************************************************
			virtual bool start_transaction() final
			{
				if (in_transaction)
					return false;

				// In out-of-transaction state file must be closed
				assert(f == nullptr);

				f = fopen(file_name.c_str(), "ab");
				if (!f)
					return false;

				in_transaction = true;

				return true;
			};

			// *******************************************************************************************
			virtual bool stop_transaction() final
			{
				if (!in_transaction)
					return false;

				// In in-transaction state file must be opened
				assert(f != nullptr);

				in_transaction = false;
				if (fclose(f) != 0)
					return false;
				f = nullptr;

				return true;
			};
		};
	}
}