#ifndef CAPTURE_OSTREAM_H_
#define CAPTURE_OSTREAM_H_
#include <ios>

namespace kmcdb::detail
{
	//mkokot_TODO: z tym jest duzy problem jak sie odpala np. wiecej instancji na raz
	//chodzi o cos takiego ze mamy
	//oryginalny streambuf z std::cout (przykladowo): A
	//zastepujemy go pierwszym wywolaniem i wsadzamy B, a pamietamy A
	//zastepujamy drugim wywoalniem i wsadzamy C, a zapamietujemy B
	//powierdzmy ze najpierw jest niszczony ten, ktory pamietal A, on ustawie std::coutowi oryginalne A
	//ale potem poleci destruktor tego co pamietal B
	//on wpisze znow B, a B juz nie istnieje
	//problem tak na prawde jest taki, ze std::cout/cerr to sa globalne zmienne
	//no i na jakims poziomie nalezaloby to kontrolowac
	//np miec globalne nasluchiwanie dla std::cout i sobie tam rejestrowac i wyrejestrowywac
	//ale to nieco komplikuje sprawy

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
