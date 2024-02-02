#include "stdafx.h"
#include "ComputationThread.h"


void ComputationThread::Start()
{
	mThread = std::thread([this]() {
		Calculate();		
	});
}

void ComputationThread::join()
{
	if (mThread.joinable()) mThread.join();
}
