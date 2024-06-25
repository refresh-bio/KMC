#pragma once

#include <cinttypes>
#include <cassert>

namespace refresh
{
	namespace io
	{
		class output_common
		{
		protected:
			bool active{false};
			size_t f_size{};

		public:
			// *******************************************************************************************
			output_common() = default;

			// *******************************************************************************************
			virtual ~output_common() = default;

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
			bool opened() const
			{
				return active;
			}

			// *******************************************************************************************
			virtual bool close() = 0;
			virtual bool put(const uint8_t c) = 0;
			virtual size_t write(const uint8_t* ptr, const size_t n_bytes) = 0;
			virtual size_t write_uint(const uint64_t val, const size_t n_bytes) = 0;
		};
	}
}