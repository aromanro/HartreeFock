#include "stdafx.h"



#include "RestrictedHartreeFock.h"
#include "UnrestrictedHartreeFock.h"

#include "HartreeFockThread.h"
#include "HartreeFockDoc.h"
#include "ChemUtils.h"

HartreeFockThread::HartreeFockThread(const Options& options, CHartreeFockDoc* doc, double start, double end, double step)
	: m_Doc(doc), m_start(start), m_end(end), m_step(step), terminate(false), converged(true),
	computeFirstAtom(false), computeSecondAtom(false), firstAtomEnergy(0), secondAtomEnergy(0)
{
	if (options.restricted && options.alphaElectrons == options.betaElectrons) {
		algorithm = new HartreeFock::RestrictedHartreeFock(options.iterations);		
	}
	else {
		HartreeFock::UnrestrictedHartreeFock *alg = new HartreeFock::UnrestrictedHartreeFock(options.iterations);
		alg->addAsymmetry = options.addAsymmetry;
		alg->asymmetry = options.asymmetry;
		algorithm = alg;
	}

	algorithm->alpha = options.alpha;
	algorithm->initGuess = options.initialGuess;

	algorithm->integralsRepository.useLotsOfMemory = options.useLotsOfMemory;

	
	CT2CA psz1(options.m_atom1);
	std::string str1(psz1);
	const unsigned int Z1 = Chemistry::ChemUtils::GetZForAtom(str1);
	CT2CA psz2(options.m_atom2);
	std::string str2(psz2);
	const unsigned int Z2 = Chemistry::ChemUtils::GetZForAtom(str2);


	angle = options.bondAngle * M_PI / 180.;

	// construct the molecule
	

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
	opt = options;
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

		results.push_back(std::make_tuple(pos, result * 27.211385056, algorithm->HOMOEnergy * 27.211385056));
		if (terminate) break;
	}

	if (!terminate && computeFirstAtom)
	{
		delete algorithm;

		Systems::Molecule atom;
		atom.atoms.push_back(atom1);
		atom.alphaElectrons = static_cast<int>(atom1.Z / 2);
		atom.betaElectrons = atom1.Z - atom.alphaElectrons;
		atom.Init();

		if (opt.restricted && atom1.Z % 2 == 0) {		
			algorithm = new HartreeFock::RestrictedHartreeFock(opt.iterations);		
		}
		else {
			HartreeFock::UnrestrictedHartreeFock *alg = new HartreeFock::UnrestrictedHartreeFock(opt.iterations);
			alg->addAsymmetry = opt.addAsymmetry;
			alg->asymmetry = opt.asymmetry;
			algorithm = alg;
		}

		algorithm->alpha = opt.alpha;
		algorithm->initGuess = opt.initialGuess;

		algorithm->integralsRepository.useLotsOfMemory = opt.useLotsOfMemory;

		algorithm->Init(&atom);

		firstAtomEnergy = algorithm->Calculate() * 27.211385056;
	}

	if (!terminate && computeSecondAtom)
	{
		delete algorithm;

		Systems::Molecule atom;
		atom.atoms.push_back(atom2);
		atom.alphaElectrons = static_cast<int>(atom2.Z / 2);
		atom.betaElectrons = atom2.Z - atom.alphaElectrons;
		atom.Init();

		if (opt.restricted && atom2.Z % 2 == 0) {		
			algorithm = new HartreeFock::RestrictedHartreeFock(opt.iterations);		
		}
		else {
			HartreeFock::UnrestrictedHartreeFock *alg = new HartreeFock::UnrestrictedHartreeFock(opt.iterations);
			alg->addAsymmetry = opt.addAsymmetry;
			alg->asymmetry = opt.asymmetry;
			algorithm = alg;
		}

		algorithm->alpha = opt.alpha;
		algorithm->initGuess = opt.initialGuess;

		algorithm->integralsRepository.useLotsOfMemory = opt.useLotsOfMemory;

		algorithm->Init(&atom);

		secondAtomEnergy = algorithm->Calculate() * 27.211385056;
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
