#pragma once

#include <thread>

class ComputationThread
{
protected:
	std::thread mThread;
	virtual ~ComputationThread();

	virtual void Calculate() = 0;
public:
	void Start();
	void join();
};

