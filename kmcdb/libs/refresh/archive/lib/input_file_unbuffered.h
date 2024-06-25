#pragma once

#include "input_common.h"

namespace refresh
{
	namespace io
	{
		// *******************************************************************************************
		//
		// *******************************************************************************************
		class input_file_unbuffered : public input_common
		{
		protected:
			FILE* f;

		public:
			// *******************************************************************************************
			input_file_unbuffered(const std::string& file_name)
			{
				std::filesystem::path fp(file_name);

				if (!std::filesystem::exists(fp))
					return;

				f_size = std::filesystem::file_size(fp);
				f_pos = 0;

				f = fopen(file_name.c_str(), "rb");

				active = f != nullptr;
			}

			// *******************************************************************************************
			virtual ~input_file_unbuffered()
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
				int c = std::getc(f);
				++f_pos;

				if (c < 0)
					return false;

				x = (uint8_t)c;

				return true;
			}

			// *******************************************************************************************
			virtual bool read_uint(const int no_bytes, uint64_t& x)
			{
				uint64_t shift = 0;
				x = 0;

				for (int i = 0; i < no_bytes; ++i)
				{
					int c = std::getc(f);
					++f_pos;

					if (c < 0)
						return false;

					x += ((uint64_t) c) << shift;
					shift += 8;
				}

				return true;
			}

			// *******************************************************************************************
			virtual size_t read(uint8_t* ptr, size_t size)
			{
				size_t readed = refresh::io::fread(ptr, 1, size, f);

				f_pos += readed;

				return readed;
			}

			// *******************************************************************************************
			virtual bool eof() const
			{
				return std::feof(f);
			}

			// *******************************************************************************************
			virtual bool seek(const size_t requested_pos) 
			{
				auto r = refresh::io::fseek(f, requested_pos, SEEK_SET) == 0;

				if (r)
					f_pos = requested_pos;

				return r;
			}
		};

		// *******************************************************************************************
		//
		// *******************************************************************************************
		class input_file_unbuffered_reopen : public input_file_unbuffered
		{
			std::string file_name;
			bool in_transaction;

		public:
			// *******************************************************************************************
			input_file_unbuffered_reopen(const std::string& file_name) :
				input_file_unbuffered(file_name),
				file_name(file_name),
				in_transaction(false)
			{
				if (f)
				{
					fclose(f);
					f = nullptr;
				}
			}

			// *******************************************************************************************
			virtual ~input_file_unbuffered_reopen()
			{
				close();
			}

			// *******************************************************************************************
			virtual bool start_transaction()
			{
				if (in_transaction)
					return false;

				assert(f == nullptr);

				f = fopen(file_name.c_str(), "rb");
				if (!f)
					return false;

				if (refresh::io::fseek(f, f_pos, SEEK_SET) != 0)
				{
					fclose(f);
					f = nullptr;

					return false;
				}

				in_transaction = true;

				return true;
			}

			// *******************************************************************************************
			virtual bool stop_transaction()
			{
				if (!in_transaction)
					return false;

				assert(f != nullptr);

				in_transaction = false;
				fclose(f);
				f = nullptr;

				return true;
			}
		};
	}
}