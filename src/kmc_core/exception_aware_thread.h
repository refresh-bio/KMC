#pragma once

#include <thread>
#include <functional>
#include <cassert>
#include "thread_cancellation_exception.h"

class CThreadExceptionCollector
{
	std::mutex mtx;
	std::vector<std::exception_ptr> collected_exceptions;
public:
	static CThreadExceptionCollector& Inst()
	{
		static CThreadExceptionCollector inst;
		return inst;
	}

	void CollectException(std::exception_ptr&& exc_ptr)
	{
		std::lock_guard<std::mutex> lck(mtx);
		collected_exceptions.push_back(std::move(exc_ptr));
	}

	void RethrowIfAnyException()
	{
		if (!collected_exceptions.empty())
		{
			auto first_exception = std::move(collected_exceptions.front());
			collected_exceptions.clear();
			std::rethrow_exception(first_exception);
		}
	}
};

class CExceptionAwareThread
{
	struct Details
	{
		std::function<void()> fun; //TODO: it should be possible to implement without std::function, some information avaiable here: https://stackoverflow.com/questions/47496358/c-lambdas-how-to-capture-variadic-parameter-pack-from-the-upper-scope
		std::thread thread;

		template<typename Callable, typename... Args>
		explicit Details(Callable&& f, Args&&... args) :
			fun(std::bind(std::forward<Callable>(f), std::forward<Args>(args)...)),
			thread([this] {
			try
			{
				fun();
			}
			catch (const CThreadCancellationException& ex) //threads is exiting because it was cancelled by critical error in other thread
			{

			}
			catch (...)
			{
				CThreadExceptionCollector::Inst().CollectException(std::current_exception());
			}
		})
		{

		}
	};
	//I'm wraping everythin into unique_ptr because if one creates a vector of CExceptionAwareThread, push_back may cause moves which may invalidate this pointer, which is needed to bound thread object with state
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

	~CExceptionAwareThread()
	{
		if(details)
			if (details->thread.joinable())
				details->thread.join();
	}
};
