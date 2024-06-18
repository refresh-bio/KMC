#ifndef CAPTURE_OSTREAM_H_
#define CAPTURE_OSTREAM_H_
#include <ios>

namespace kmcdb::detail
{
	class CaptureOstream : public std::streambuf
	{
		std::ostream* src;

		std::streambuf* org;

		std::string& capture_to;

	public:
		CaptureOstream(std::ostream& src, std::string& capture_to) :
			src(&src),
			org(src.rdbuf()),
			capture_to(capture_to)
		{
			src.rdbuf(this);
		}

		int overflow(int_type c) override
		{
			if (c == EOF)
				return c;

			capture_to += static_cast<char>(c);
			return org->sputc(static_cast<char>(c));
		}
		std::streamsize xsputn(const char* ptr, std::streamsize cnt) override
		{
			capture_to.insert(capture_to.end(), ptr, ptr + cnt);
			return org->sputn(ptr, cnt);
		}

		~CaptureOstream() override
		{
			try
			{
				src->rdbuf(org);
			}
			catch (...)
			{

			}
		}
	};
}
#endif // ! CAPTURE_OSTREAM_H_
