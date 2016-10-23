#pragma once

#include <thread>

class ComputationThread
{
public:
	void Start();
protected:
	std::thread mThread;
	virtual ~ComputationThread();

	virtual void Calculate() = 0;
public:
	void join();
};

