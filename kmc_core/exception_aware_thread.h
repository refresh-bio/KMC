#pragma once

#include <thread>
#include <functional>
#include <cassert>


class CExceptionAwareThread
{
	struct Details
	{
		std::exception_ptr exc_ptr;
		std::function<void()> fun; //TODO: it should be possible to implement without std::function, some information avaiable here: https://stackoverflow.com/questions/47496358/c-lambdas-how-to-capture-variadic-parameter-pack-from-the-upper-scope
		std::thread thread;

		template<typename Callable, typename... Args>
		explicit Details(Callable&& f, Args&&... args) :
			fun(std::bind(std::forward<Callable>(f), std::forward<Args>(args)...)),
			thread([this] {
			try
			{
				//	std::cerr << "running: " << this << "\n";
				fun();
			}
			catch (...)
			{
				auto ex = std::current_exception();
				exc_ptr = ex;
				int a = 5;
			}
		})
		{

		}
	};
	//I'm wraping everythin into unique_ptr because if one creates a vector of CExceptionAwareThread, push_back may cause moves which may invalidate this pointer, which is needed to bound thread with state (exc_ptr)
	std::unique_ptr<Details> details;
public:
	CExceptionAwareThread() = default;
	CExceptionAwareThread(CExceptionAwareThread&& rhs) = default;
	CExceptionAwareThread& operator=(CExceptionAwareThread&& rhs) = default;
	CExceptionAwareThread(const CExceptionAwareThread& rhs) = delete;
	CExceptionAwareThread& operator=(const CExceptionAwareThread& rhs) = delete;

	template<typename Callable, typename... Args>
	explicit CExceptionAwareThread(Callable&& f, Args&&... args) :
		details(std::make_unique<Details>(std::forward<Callable>(f), std::forward<Args>(args)...))
	{

	}

	void join()
	{
		details->thread.join();
	}

	void RethrowIfException()
	{
		assert(details->thread.joinable() == false); //must be joined
		if (details->exc_ptr)
			std::rethrow_exception(details->exc_ptr);
	}
};