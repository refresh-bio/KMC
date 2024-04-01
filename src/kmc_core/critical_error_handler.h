#pragma once

#include <set>
#include <mutex>
#include "thread_cancellation_exception.h" //TODO: moze ten wyjatek zdefiniowac tutaj?
#include <condition_variable>
#include <iostream>

class CThrowingOnCancelConditionVariable
{
	std::condition_variable cv;
	bool canceled = false;
public:
	CThrowingOnCancelConditionVariable();
	~CThrowingOnCancelConditionVariable();
	template<typename Predicate>
	void wait(std::unique_lock<std::mutex>& lock, Predicate predicate)
	{
		cv.wait(lock, [this, &predicate] {
			if (canceled)
				throw CThreadCancellationException{};
			return predicate();
		});
	}

	void notify_all()
	{
		cv.notify_all();
	}

	void notify_one()
	{
		cv.notify_one();
	}

	void CancelAllThreads()
	{
		canceled = true;
		notify_all();
	}
};

class CCriticalErrorHandler
{
	std::set<CThrowingOnCancelConditionVariable*> observed_condition_variables;
	std::mutex mtx;

	void cancelAllThreads()
	{
		std::unique_lock<std::mutex> lck(mtx);
		for (auto cv : observed_condition_variables)
			cv->CancelAllThreads();
	}
public:
	static CCriticalErrorHandler& Inst()
	{
		static CCriticalErrorHandler inst;
		return inst;
	}

	void RegisterConditionVariable(CThrowingOnCancelConditionVariable* cv)
	{
		std::unique_lock<std::mutex> lck(mtx);
		observed_condition_variables.insert(cv);
	}

	void UnregisterConditionVariable(CThrowingOnCancelConditionVariable* cv)
	{
		std::unique_lock<std::mutex> lck(mtx);
		observed_condition_variables.erase(cv);
	}

	void HandleCriticalError(const std::string& msg)
	{
		//std::cerr << msg << "\n";
		//exit(1);
		cancelAllThreads();
		throw std::runtime_error(msg);
	}
};

inline CThrowingOnCancelConditionVariable::CThrowingOnCancelConditionVariable()
{
	CCriticalErrorHandler::Inst().RegisterConditionVariable(this);
}

inline CThrowingOnCancelConditionVariable::~CThrowingOnCancelConditionVariable()
{
	CCriticalErrorHandler::Inst().UnregisterConditionVariable(this);
}
