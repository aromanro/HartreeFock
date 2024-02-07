#pragma once

#include <thread>

class ComputationThread
{
public:
	void Start();
	void join();

protected:
	std::thread mThread;
	virtual ~ComputationThread() = default;

	virtual void Calculate() = 0;
};

