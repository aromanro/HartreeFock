#pragma once
#include "ComputationThread.h"

#include "Options.h"


#include "Molecule.h"

#include <atomic>
#include <vector>

namespace HartreeFock {
	class HartreeFockAlgorithm;
};

class CHartreeFockDoc;

class HartreeFockThread :
	public ComputationThread
{
protected:
	std::atomic_bool terminate;
	HartreeFock::HartreeFockAlgorithm *algorithm;
	CHartreeFockDoc* m_Doc;

	Systems::Molecule molecule;

	double m_start;
	double m_end;
	double m_step;

	double angle;

	bool converged;
public:
	HartreeFockThread(const Options& options, CHartreeFockDoc* doc, double start, double end, double step);
	virtual ~HartreeFockThread();

	std::vector<std::tuple<double, double, double>> results;

	virtual void Calculate();
	void Terminate();
	bool Converged() const;
};

