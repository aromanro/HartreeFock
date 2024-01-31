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

	Systems::AtomWithShells atom1, atom2;
	Options opt;

public:
	HartreeFockThread(const Options& options, CHartreeFockDoc* doc, const double start, const double end, const double step);
	~HartreeFockThread() override;

	std::vector<std::tuple<double, double, double>> results;

	void Calculate() override;
	void Terminate();
	bool Converged() const;

	bool computeFirstAtom;
	bool computeSecondAtom;

	double firstAtomEnergy;
	double secondAtomEnergy;

private:
	void ComputeAtoms();
	double ComputeAtom(const Systems::AtomWithShells& atom);
};

