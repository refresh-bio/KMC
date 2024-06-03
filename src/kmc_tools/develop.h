#ifndef _DEVELOP_H
#define _DEVELOP_H

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include "defs.h"
#include "timer.h"

#ifdef ENABLE_LOGGER




class CLoger
{
	std::map<std::string, std::map<void*, std::vector<double>>> operations;

public:
	static CLoger& GetLogger()
	{
		static CLoger logger;
		return logger;
	}

	void add_to_operators_times(void* addr, double time)
	{
		operations["operators_times"][addr].push_back(time);
	}

	void add_to_wait_for_inputdata(void* addr, double time)
	{
		operations["wait_for_input"][addr].push_back(time);
	}

	void log_operation(const std::string& name, void* addr, double time)
	{
		operations[name][addr].push_back(time);
	}

	void print_stats()
	{
		std::cout << "--operators_times--\n";

		for (auto& oper : operations)
		{
			std::cout << "operation : -----------" << oper.first << "-----------\n";
			for (const auto& _operator : oper.second)
			{
				std::cout << "operator at " << _operator.first << "\n";
				double sum = 0;
				for (const auto& time : _operator.second)
				{
					sum += time;
				}
				std::cout << "time: " << sum << "\nmean: " << (sum / _operator.second.size()) << ", no of entries: " << _operator.second.size() << "\n";
			}
		}

	}
};
#endif


#endif