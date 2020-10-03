#include "stdafx.h"

#include "MathUtils.h"

#include "GaussianOrbital.h"
#include "Molecule.h"
#include "IntegralsRepository.h"
#include "QuantumMatrix.h"

#include "BoysFunction.h"

#include "Test.h"

#include "RestrictedHartreeFock.h"
#include "UnrestrictedHartreeFock.h"

#include "Basis.h"
#include "ChemUtils.h"

#include <sstream>

#include <fstream>
#include <iostream>
#include <iomanip>


// will test values against examples from here:
// https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2303
// uses STO3G, so that's what is loaded

Test::Test()
{
	basisSTO3G.Load("sto3g.txt");
}

void Test::OutputMatricesForAtom(const std::string& atomName, const std::string& basisSetName, const std::string& fileName)
{
	Chemistry::Basis basisCustom;
	basisCustom.Load(basisSetName);

	Systems::AtomWithShells theAtom;

	const int Z = Chemistry::ChemUtils::GetZForAtom(atomName);
	for (auto& atom : basisCustom.atoms)
	{
		if (Z == atom.Z)
		{
			theAtom = atom;
			break;
		}
	}

	Systems::Molecule atom;

	atom.atoms.push_back(theAtom);
	atom.Init();

	std::ofstream file(fileName);

	OutputMatrices(atom, file);
}


void Test::OutputMatrices(Systems::Molecule& molecule, std::ofstream& file)
{
	HartreeFock::HartreeFockAlgorithm* hartreeFock;
	if (molecule.alphaElectrons % 2 == 0 && molecule.betaElectrons % 2 == 0)
		hartreeFock = new HartreeFock::RestrictedHartreeFock();
	else
		hartreeFock = new HartreeFock::UnrestrictedHartreeFock();

	hartreeFock->Init(&molecule);

	// now we have the molecule

	//hartreeFock->integralsRepository.useLotsOfMemory = false;

	// now check the overlap matrix

	file << "Overlap Matrix:\n";

	for (int i = 0; i < hartreeFock->overlapMatrix.matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < hartreeFock->overlapMatrix.matrix.cols(); ++j)
			line << hartreeFock->overlapMatrix.matrix(i, j) << "\t\t\t";

		line << std::endl;

		file << line.str();
	}

	file << "\nKinetic Matrix:\n";

	for (int i = 0; i < hartreeFock->kineticMatrix.matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < hartreeFock->kineticMatrix.matrix.cols(); ++j)
			line << hartreeFock->kineticMatrix.matrix(i, j) << "\t\t\t";

		line << std::endl;

		file << line.str();
	}

	file << "\nNuclear Matrix:\n";

	for (int i = 0; i < hartreeFock->nuclearMatrix.matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < hartreeFock->nuclearMatrix.matrix.cols(); ++j)
			line << hartreeFock->nuclearMatrix.matrix(i, j) << "\t\t\t";

		line << std::endl;

		file << line.str();
	}

	file << std::endl;


	/*
	for (Orbitals::QuantumNumbers::QuantumNumbers QN(0, 0, 0); QN <= 8; ++QN)
	{
	TRACE("(%d, %d, %d), Canonical Index: %d, Total Index: %d\n", QN.l, QN.m, QN.n, QN.GetCanonicalIndex(), QN.GetTotalCanonicalIndex());
	}
	*/


	// initialize


	const Eigen::MatrixXd h = hartreeFock->kineticMatrix.matrix + hartreeFock->nuclearMatrix.matrix;

	file << "\nh Matrix:\n";

	for (int i = 0; i < h.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < h.cols(); ++j)
			line << h(i, j) << "\t\t\t";

		line << std::endl;

		file << line.str();
	}

	file << "\nV Matrix:\n";

	for (int i = 0; i < hartreeFock->V.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < hartreeFock->V.cols(); ++j)
			line << hartreeFock->V(i, j) << "\t\t\t";

		line << std::endl;

		file << line.str();
	}


	file << "\ne-e Matrix:\n";

	file.precision(5);
	file << std::setprecision(5);


	const unsigned int nrorbs = hartreeFock->integralsRepository.getMolecule()->CountNumberOfContractedGaussians();

	for (unsigned int i = 0; i < nrorbs; i++)
		for (unsigned int j = 0; j < nrorbs; j++)
			for (unsigned int k = 0; k < nrorbs; k++)
				for (unsigned int l = 0; l < nrorbs; l++) {
					double val = hartreeFock->integralsRepository.getElectronElectron(i, j, k, l);

					file << "(" << i + 1 << "," << j + 1 << "," << k + 1 << "," << l + 1 << ") = " << std::setprecision(5) << std::setiosflags(std::ios::fixed) << std::setw(6) << val << std::endl;
				}


	delete hartreeFock;
}


// TODO: try to load the provided file at https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2303 and do comparisons in the code
// even the geometry of the molecule could be loaded from the file

void Test::TestWater(const std::string& fileName)
{
	Systems::AtomWithShells H1,H2,O;

	for (auto& atom : basisSTO3G.atoms)
	{
		if (1 == atom.Z)
		{
			H1 = H2 = atom;
		}
		else if (8 == atom.Z)
			O = atom;
	}

	Systems::Molecule molecule;

	O.position.X = 0;
	O.position.Y = -0.143225816552;
	O.position.Z = 0;

	H1.position.X = 1.638036840407;
	H1.position.Y = 1.136548822547;
	H1.position.Z = 0;

	H2.position.X = -1.638036840407;
	H2.position.Y = 1.136548822547;
	H2.position.Z = 0;

	molecule.atoms.push_back(O);
	molecule.atoms.push_back(H1);
	molecule.atoms.push_back(H2);
	molecule.Init();

	std::ofstream file(fileName);

	OutputMatrices(molecule, file);
}
