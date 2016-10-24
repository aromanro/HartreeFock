#include "stdafx.h"

#include "RestrictedHartreeFock.h"
#include "UnrestrictedHartreeFock.h"

#include "HartreeFockThread.h"
#include "HartreeFockDoc.h"
#include "ChemUtils.h"

HartreeFockThread::HartreeFockThread(const Options& options, CHartreeFockDoc* doc, double start, double end, double step)
	: m_Doc(doc), m_start(start), m_end(end), m_step(step), terminate(false), converged(true)
{
	if (options.restricted && options.alphaElectrons == options.betaElectrons) {
		algorithm = new HartreeFock::RestrictedHartreeFock(options.iterations);		
	}
	else {
		algorithm = new HartreeFock::UnrestrictedHartreeFock(options.iterations);
	}

	algorithm->alpha = options.alpha;
	algorithm->initGuess = options.initialGuess;

	algorithm->integralsRepository.useLotsOfMemory = options.useLotsOfMemory;

	
	CT2CA psz1(options.m_atom1);
	std::string str1(psz1);
	unsigned int Z1 = Chemistry::ChemUtils::GetZForAtom(str1);
	CT2CA psz2(options.m_atom2);
	std::string str2(psz2);
	unsigned int Z2 = Chemistry::ChemUtils::GetZForAtom(str2);


	angle = options.bondAngle * M_PI / 180.;

	// construct the molecule
	Systems::AtomWithShells atom1, atom2;

	if (options.basis)
	{
		for (const auto &atom : doc->basisSTO6G.atoms)
		{
			if (Z1 == atom.Z) atom1 = atom;
			if (Z2 == atom.Z) atom2 = atom;
		}
	}
	else
	{
		for (const auto &atom : doc->basisSTO3G.atoms)
		{
			if (Z1 == atom.Z) atom1 = atom;
			if (Z2 == atom.Z) atom2 = atom;
		}
	}

	molecule.atoms.push_back(atom1);
	molecule.atoms.push_back(atom2);
	molecule.atoms[1].SetCenterForShells();
	if (options.twoAtom1) molecule.atoms.push_back(atom1);

	molecule.alphaElectrons = options.alphaElectrons;
	molecule.betaElectrons = options.betaElectrons;

	molecule.SetIDs();
}


HartreeFockThread::~HartreeFockThread()
{
	if (mThread.joinable()) mThread.join();

	delete algorithm;
}

void HartreeFockThread::Calculate()
{
	const double bohr = 0.5291772106712;

	for (double pos = m_start; pos < m_end; pos += m_step)
	{
		double dist = pos / bohr;

		// adjust the molecule coordinates
		if (molecule.atoms.size() <= 2)
		{
			molecule.atoms[0].position.X = -dist / 2.;
			molecule.atoms[0].SetCenterForShells();
			molecule.atoms[1].position.X = dist / 2.;
			molecule.atoms[1].SetCenterForShells();
		}
		else
		{
			molecule.atoms[0].position.X = dist * cos(angle / 2.);
			molecule.atoms[0].position.Y = dist * sin(angle / 2.);

			molecule.atoms[2].position.X = molecule.atoms[0].position.X;
			molecule.atoms[2].position.Y = -molecule.atoms[0].position.Y;

			molecule.atoms[0].SetCenterForShells();
			molecule.atoms[2].SetCenterForShells();
		}


		algorithm->Init(&molecule);

		

		double result = algorithm->Calculate();

		if (!algorithm->converged) converged = false;

		results.push_back(std::make_pair(pos, result * 27.211385056));
		if (terminate) break;
	}

	--m_Doc->runningThreads;
}


void HartreeFockThread::Terminate()
{
	algorithm->terminate = true;
	terminate = true;
}


bool HartreeFockThread::Converged() const
{
	return converged;
}
