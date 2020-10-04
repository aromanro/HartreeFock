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


void Test::OutputMatrix(const Eigen::MatrixXd& matrix, std::ofstream& file)
{
	for (int i = 0; i < matrix.rows(); ++i)
	{
		std::stringstream line;
		line.precision(6);
		line.setf(std::ios_base::fixed);

		for (int j = 0; j < matrix.cols(); ++j)
			line << matrix(i, j) << "\t\t\t";

		line << std::endl;

		file << line.str();
	}

	file << std::endl;
}



void Test::OutputMatrices(Systems::Molecule& molecule, std::ofstream& file, const std::string& sfileName, const std::string& tfileName, const std::string& vfileName, const std::string& erifileName, const double expectedNucEnergy)
{
	HartreeFock::HartreeFockAlgorithm* hartreeFock;
	const bool restricted = molecule.alphaElectrons == molecule.betaElectrons;
	if (restricted)
		hartreeFock = new HartreeFock::RestrictedHartreeFock();
	else
		hartreeFock = new HartreeFock::UnrestrictedHartreeFock();

	hartreeFock->UseDIIS = false;
	hartreeFock->alpha = 0.5;
	hartreeFock->initGuess = 0;

	hartreeFock->Init(&molecule);


	file.precision(15);
	file << "Nuclear repulsion energy: " << hartreeFock->nuclearRepulsionEnergy;

	if (expectedNucEnergy > 1E-15)
	{
		file << " Expected: " << expectedNucEnergy;
	}

	file << std::endl << std::endl;

	// now we have the molecule

	//hartreeFock->integralsRepository.useLotsOfMemory = false;

	// now check the overlap matrix

	file << "Overlap Matrix:\n";

	OutputMatrix(hartreeFock->overlapMatrix.matrix, file);
	CheckDifferences(hartreeFock->overlapMatrix.matrix, sfileName, file);

	file << "\nKinetic Matrix:\n";

	OutputMatrix(hartreeFock->kineticMatrix.matrix, file);
	CheckDifferences(hartreeFock->kineticMatrix.matrix, tfileName, file);

	file << "\nNuclear Matrix:\n";

	OutputMatrix(hartreeFock->nuclearMatrix.matrix, file);
	CheckDifferences(hartreeFock->nuclearMatrix.matrix, vfileName, file);

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

	OutputMatrix(h, file);

	file << "\nV Matrix:\n";

	OutputMatrix(hartreeFock->V, file);


	file << "\ne-e Matrix:\n";

	file.precision(5);
	file << std::setprecision(5);


	const unsigned int nrorbs = hartreeFock->integralsRepository.getMolecule()->CountNumberOfContractedGaussians();

	for (unsigned int i = 0; i < nrorbs; i++)
		for (unsigned int j = 0; j < nrorbs; j++)
			for (unsigned int k = 0; k < nrorbs; k++)
				for (unsigned int l = 0; l < nrorbs; l++) {
					const double val = hartreeFock->integralsRepository.getElectronElectron(i, j, k, l);

					file << "(" << i + 1 << "," << j + 1 << "," << k + 1 << "," << l + 1 << ") = " << std::setprecision(5) << std::setiosflags(std::ios::fixed) << std::setw(6) << val << std::endl;
				}


	CheckDifferences(hartreeFock->integralsRepository, erifileName, file);


	file << "\nChecking if orthogonalization really works, transformed overlap matrix:" << std::endl;

	const Eigen::MatrixXd Stransf = hartreeFock->Vt * hartreeFock->overlapMatrix.matrix * hartreeFock->V;

	OutputMatrix(Stransf, file);


	file << "\nV:" << std::endl;
	OutputMatrix(hartreeFock->V, file);

	file << "\nVt:" << std::endl;
	OutputMatrix(hartreeFock->Vt, file);

	// now check 'step 0' (initialization) 

	// check only the restricted algo for now
	if (restricted)
	{
		Eigen::MatrixXd FockMatrix;

		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->InitFockMatrix(0, FockMatrix);

		Eigen::MatrixXd FockTransformed = hartreeFock->Vt * FockMatrix * hartreeFock->V; // orthogonalize

		file << "\nInitial transformed Fock Matrix:\n";

		OutputMatrix(FockTransformed, file);

		const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(FockTransformed);
		const Eigen::MatrixXd& Cprime = es.eigenvectors();

		file << "\nInitial MO Coefficients:\n";

		OutputMatrix(Cprime, file);

		Eigen::MatrixXd Ce = hartreeFock->V * Cprime; // transform back the eigenvectors into the original non-orthogonalized AO basis

		Eigen::MatrixXd DensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				for (unsigned int vec = 0; vec < ((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied.size(); ++vec) // only eigenstates that are occupied 
					if (((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied[vec]) DensityMatrix(i, j) += Ce(i, vec) * Ce(j, vec); // 2 is for the number of electrons in the eigenstate, it's the restricted Hartree-Fock


		file << "\nInitial Density Matrix (no occupation multiplication):\n";

		OutputMatrix(DensityMatrix, file);

		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->eigenvals = es.eigenvalues();

		DensityMatrix *= 2;
		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->CalculateEnergy(((HartreeFock::RestrictedHartreeFock*)hartreeFock)->eigenvals, DensityMatrix/*, FockMatrix*/);

		file.precision(14);
		file << "\nInitial Hartree-Fock electronic energy: " << hartreeFock->GetTotalEnergy() - hartreeFock->nuclearRepulsionEnergy << " Hartrees" << std::endl;

		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->DensityMatrix = DensityMatrix;


		double rmsD = 0;
		double deltaE = 0;
		for (int step = 1; step < 100; ++step)
		{
			((HartreeFock::RestrictedHartreeFock*)hartreeFock)->InitFockMatrix(step, FockMatrix);

			if (1 == step)
			{
				file.precision(5);
				file << "\nThe first iteration Fock matrix in the AO basis: \n";

				OutputMatrix(FockMatrix, file);
			}

			// just copy/paste of the code from Restricted Hartree Fock class, with slight changes to work here, too:

			FockTransformed = hartreeFock->Vt * FockMatrix * hartreeFock->V; // orthogonalize
			const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esl(FockTransformed);
			const Eigen::MatrixXd& Cprime = esl.eigenvectors();
			Ce = hartreeFock->V * Cprime; // transform back the eigenvectors into the original non-orthogonalized AO basis

			DensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

			for (int i = 0; i < h.rows(); ++i)
				for (int j = 0; j < h.cols(); ++j)
					for (unsigned int vec = 0; vec < ((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied.size(); ++vec) // only eigenstates that are occupied 
						if (((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied[vec]) DensityMatrix(i, j) += 2. * Ce(i, vec) * Ce(j, vec); // 2 is for the number of electrons in the eigenstate, it's the restricted Hartree-Fock


			//**************************************************************************************************************

			((HartreeFock::RestrictedHartreeFock*)hartreeFock)->eigenvals = esl.eigenvalues();

			const double oldEnergy = hartreeFock->GetTotalEnergy();

			((HartreeFock::RestrictedHartreeFock*)hartreeFock)->CalculateEnergy(((HartreeFock::RestrictedHartreeFock*)hartreeFock)->eigenvals, DensityMatrix/*, FockMatrix*/);

			const double totalEnergy = hartreeFock->GetTotalEnergy();
			deltaE = totalEnergy - oldEnergy;

			// calculate rms for differences between new and old density matrices, it can be used to check for convergence, too
			Eigen::MatrixXd rmsDensityMatricesDif = DensityMatrix - ((HartreeFock::RestrictedHartreeFock*)hartreeFock)->DensityMatrix;

			rmsDensityMatricesDif *= 0.5; // there is a difference compared with the one from the reference, they don't multiply it with 2 as I do, so this is done for comparison purposes

			rmsD = sqrt(rmsDensityMatricesDif.cwiseProduct(rmsDensityMatricesDif).sum());

			file.precision(12);
			file << step << "\t" << totalEnergy - hartreeFock->nuclearRepulsionEnergy << "\t" << totalEnergy << "\t" << deltaE << "\t" << rmsD << std::endl;


			if (rmsD < 1E-12 && abs(deltaE)) break;

			// ***************************************************************************************************
			// go to the next density matrix

			((HartreeFock::RestrictedHartreeFock*)hartreeFock)->DensityMatrix = DensityMatrix;
		}
	}

	delete hartreeFock;
}



void Test::CheckDifferences(const Eigen::MatrixXd& matrix, const std::string& matrixFileName, std::ofstream& file)
{
	if (matrixFileName.empty() || !file) return;
	
	std::ifstream mfile(matrixFileName);
	if (!mfile) return;

	file << "\nDifferences, if there are any:\n";

	int countDifferences = 0;


	std::string line;
	while (std::getline(mfile, line))
	{
		std::istringstream lineStream(line);

		int i, j;
		double val;

		lineStream >> i >> j >> val;
		
		// indexed from 1 in the files
		--i;
		--j;


		if (i > matrix.rows() || j > matrix.cols())
			file << "Out of bonds: " << i << " " << j << std::endl;
		else
		{
			file.precision(6);
			if (abs(val - matrix(i, j)) > 1E-5)
			{
				file << "Differences, matrix(" << i << "," << j << ")=" << matrix(i, j) << " Expected: " << val << std::endl;
				++countDifferences;
			}
			//else file << "NO Differences, matrix(" << i << "," << j << ")=" << matrix(i, j) << " Expected: " << val << std::endl;
		}
	}

	if (0 == countDifferences) file << "No differences!" << std::endl;
}


void Test::CheckDifferences(const GaussianIntegrals::IntegralsRepository& repo, const std::string& eriFileName, std::ofstream& file)
{
	if (eriFileName.empty() || !file) return;

	std::ifstream erifile(eriFileName);
	if (!erifile) return;

	file << "\nDifferences, if there are any:\n";


	int countDifferences = 0;


	std::string line;
	while (std::getline(erifile, line))
	{
		std::istringstream lineStream(line);

		int i, j, k, l;
		double val;

		lineStream >> i >> j >> k >> l >> val;

		// indexed from 1 in the files
		--i;
		--j;
		--k;
		--l;

		double calcValue = repo.getElectronElectron(i, j, k, l);

		file.precision(6);
		if (abs(val - calcValue) > 1E-5)
		{
			file << "Differences, integral(" << i + 1 << "," << j + 1 << "," << k + 1 << "," << l + 1 << ")=" << calcValue << " Expected: " << val << std::endl;
			++countDifferences;
		}
		//else file << "NO Differences, integral(" << i + 1 << "," << j + 1 << "," << k + 1 << "," << l + 1 << ")=" << calcValue << " Expected: " << val << std::endl;
	}

	if (0 == countDifferences) file << "No differences!" << std::endl;
}



// TODO: try to load the provided file at https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2303 and do comparisons in the code
// even the geometry of the molecule could be loaded from the file

void Test::TestWater(const std::string& fileName, const std::string& sfileName, const std::string& tfileName, const std::string& vfileName, const std::string& erifileName)
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

	OutputMatrices(molecule, file, sfileName, tfileName, vfileName, erifileName, 8.002367061810450);
}


void Test::TestMethane(const std::string& fileName, const std::string& sfileName, const std::string& tfileName, const std::string& vfileName, const std::string& erifileName)
{
	Systems::AtomWithShells H1, H2, H3, H4, C;

	for (auto& atom : basisSTO3G.atoms)
	{
		if (1 == atom.Z)
		{
			H1 = H2 = H3 = H4 = atom;
		}
		else if (6 == atom.Z)
			C = atom;
	}

	Systems::Molecule molecule;

	C.position.X = 0;
	C.position.Y = 0;
	C.position.Z = 0;

	H1.position.X = 1.183771681898;
	H1.position.Y = -1.183771681898;
	H1.position.Z = -1.183771681898;

	H2.position.X = 1.183771681898;
	H2.position.Y = 1.183771681898;
	H2.position.Z = 1.183771681898;

	H3.position.X = -1.183771681898;
	H3.position.Y = 1.183771681898;
	H3.position.Z = -1.183771681898;

	H4.position.X = -1.183771681898;
	H4.position.Y = -1.183771681898;
	H4.position.Z = 1.183771681898;

	molecule.atoms.push_back(C);
	molecule.atoms.push_back(H1);
	molecule.atoms.push_back(H2);
	molecule.atoms.push_back(H3);
	molecule.atoms.push_back(H4);
	molecule.Init();

	std::ofstream file(fileName);

	OutputMatrices(molecule, file, sfileName, tfileName, vfileName, erifileName, 13.497304462036480);
}
