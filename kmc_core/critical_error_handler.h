#pragma once

#include <set>
#include <mutex>

class IFinishableQueue
{
public:
	virtual void ForceToFinish() = 0;
	virtual ~IFinishableQueue() = default;
};

class CCriticalErrorHandler
{
	std::set<IFinishableQueue*> observed_queues;
	std::mutex mtx;

	void finishAllQueues()
	{
		std::unique_lock<std::mutex> lck(mtx);
		for (auto queue : observed_queues)
			queue->ForceToFinish();
	}
public:
	static CCriticalErrorHandler& Inst()
	{
		static CCriticalErrorHandler inst;
		return inst;
	}

	void RegisterQueue(IFinishableQueue* queue)
	{
		std::unique_lock<std::mutex> lck(mtx);
		observed_queues.insert(queue);
	}

	void UnregisterQueue(IFinishableQueue* queue)
	{
		std::unique_lock<std::mutex> lck(mtx);
		observed_queues.erase(queue);
	}

	void HandleCriticalError(const std::string& msg)
	{
		std::cerr << msg << "\n";
		exit(1);
		//finishAllQueues();
		//throw std::runtime_error(msg);
	}
};
