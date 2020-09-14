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
	// now we have the molecule
	GaussianIntegrals::IntegralsRepository integralsRepository(&molecule);

	//integralsRepository.useLotsOfMemory = false;

	// now check the overlap matrix
	Matrices::OverlapMatrix overlapMatrix(&integralsRepository);

	TRACE("Overlap Matrix:\n");
	file << "Overlap Matrix:\n";

	for (int i = 0; i < overlapMatrix.matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < overlapMatrix.matrix.cols(); ++j)
			line << overlapMatrix.matrix(i, j) << "\t\t\t";

		line << std::endl;
		TRACE(line.str().c_str());

		file << line.str();
	}

	Matrices::KineticMatrix kineticMatrix(&integralsRepository);

	TRACE("\nKinetic Matrix:\n");
	file << "\nKinetic Matrix:\n";

	for (int i = 0; i < kineticMatrix.matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < kineticMatrix.matrix.cols(); ++j)
			line << kineticMatrix.matrix(i, j) << "\t\t\t";

		line << std::endl;
		TRACE(line.str().c_str());

		file << line.str();
	}

	TRACE("\nNuclear Matrix:\n");
	file << "\nNuclear Matrix:\n";

	Matrices::NuclearMatrix nuclearMatrix(&integralsRepository);

	for (int i = 0; i < nuclearMatrix.matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < nuclearMatrix.matrix.cols(); ++j)
			line << nuclearMatrix.matrix(i, j) << "\t\t\t";

		line << std::endl;
		TRACE(line.str().c_str());

		file << line.str();
	}

	TRACE("\n");
	file << std::endl;


	/*
	for (Orbitals::QuantumNumbers::QuantumNumbers QN(0, 0, 0); QN <= 8; ++QN)
	{
	TRACE("(%d, %d, %d), Canonical Index: %d, Total Index: %d\n", QN.l, QN.m, QN.n, QN.GetCanonicalIndex(), QN.GetTotalCanonicalIndex());
	}
	*/


	integralsRepository.ClearElectronElectronIntermediaries();
	integralsRepository.ClearMatricesMaps();

	// initialize


	Eigen::MatrixXd h = kineticMatrix.matrix + nuclearMatrix.matrix;

	TRACE("\nh Matrix:\n");
	file << "\nh Matrix:\n";

	for (int i = 0; i < kineticMatrix.matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < kineticMatrix.matrix.cols(); ++j)
			line << h(i, j) << "\t\t\t";

		line << std::endl;
		TRACE(line.str().c_str());

		file << line.str();
	}


	TRACE("\ne-e Matrix:\n");

	file << "\ne-e Matrix:\n";

	file.precision(5);
	file << std::setprecision(5);
	integralsRepository.CalculateElectronElectronIntegrals();
	integralsRepository.ClearAllMaps();


	unsigned int nrorbs = integralsRepository.getMolecule()->CountNumberOfContractedGaussians();

	for (unsigned int i = 0; i < nrorbs; i++)
		for (unsigned int j = 0; j < nrorbs; j++)
			for (unsigned int k = 0; k < nrorbs; k++)
				for (unsigned int l = 0; l < nrorbs; l++) {
					double val = integralsRepository.getElectronElectron(i, j, k, l);

					TRACE("(%d, %d, %d, %d) = %f\n", i + 1, j + 1, k + 1, l + 1, val);

					file << "(" << i + 1 << "," << j + 1 << "," << k + 1 << "," << l + 1 << ") = " << std::setprecision(5) << std::setiosflags(std::ios::fixed) << std::setw(6) << val << std::endl;
				}
}
