#pragma once

#include <cinttypes>
#include <cassert>

namespace refresh
{
	namespace io
	{
		class input_common
		{
		protected:
			bool active{false};
			size_t f_size{};
			size_t f_pos{};

		public:
			// *******************************************************************************************
			input_common() = default;

			// *******************************************************************************************
			virtual ~input_common() = default;

			// *******************************************************************************************
			virtual bool start_transaction() 
			{ 
				return true;						// For most types do nothing
			};
			
			// *******************************************************************************************
			virtual bool stop_transaction() 
			{ 
				return true;						// For most types do nothing
			};		

			// *******************************************************************************************
			size_t file_size() const 
			{
				return f_size;
			}

			// *******************************************************************************************
			size_t get_pos() const 
			{
				return f_pos;
			}

			// *******************************************************************************************
			bool opened() const
			{
				return active;
			}

			// *******************************************************************************************
			virtual bool close() = 0;
			virtual bool get(uint8_t& x) = 0;
			virtual bool read_uint(const int no_bytes, uint64_t& x) = 0;
			virtual size_t read(uint8_t* ptr, size_t size) = 0;
			virtual bool eof() const = 0;
			virtual bool seek(const size_t requested_pos) = 0;
		};
	}
}
