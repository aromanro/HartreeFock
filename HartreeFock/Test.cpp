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
#include "RestrictedCCSD.h"

#include "RestrictedConfigurationIInteractionSingles.h"

#include "Basis.h"
#include "ChemUtils.h"

#include "Constants.h"

#include <sstream>

#include <fstream>
#include <iostream>
#include <iomanip>


// will test values against examples from here:
// https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2303
// uses STO3G, so that's what is loaded

Test::Test(const std::string& basisFile)
{
	basis.Load(basisFile);
}

void Test::OutputMatricesForAtom(const std::string& atomName, const std::string& basisSetName, const std::string& fileName)
{
	Chemistry::Basis basisCustom;
	basisCustom.Load(basisSetName);

	Systems::AtomWithShells theAtom;

	const int Z = Chemistry::ChemUtils::GetZForAtom(atomName);
	for (auto& atom : basisCustom.atoms)
	{
		if (static_cast<const unsigned int>(Z) == atom.Z)
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



void Test::OutputMatrices(Systems::Molecule& molecule, std::ofstream& file, const std::string& sfileName, const std::string& tfileName, const std::string& vfileName, const std::string& erifileName, bool useDIIS, const double expectedNucEnergy)
{
	HartreeFock::HartreeFockAlgorithm* hartreeFock;
	const bool restricted = molecule.alphaElectrons == molecule.betaElectrons;
	if (restricted)
		//hartreeFock = new HartreeFock::RestrictedHartreeFock();
		hartreeFock = new HartreeFock::RestrictedCCSD();
	else
	{
		hartreeFock = new HartreeFock::UnrestrictedHartreeFock();
		((HartreeFock::UnrestrictedHartreeFock*)hartreeFock)->addAsymmetry = false;
	}

	hartreeFock->UseDIIS = useDIIS;
	hartreeFock->alpha = 0.5;
	hartreeFock->initGuess = 0.7;

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

		if (useDIIS) ((HartreeFock::RestrictedHartreeFock*)hartreeFock)->DIISStep(0, FockMatrix);

		Eigen::MatrixXd FockTransformed = hartreeFock->Vt * FockMatrix * hartreeFock->V; // orthogonalize

		file << "\nInitial transformed Fock Matrix:\n";

		OutputMatrix(FockTransformed, file);

		const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(FockTransformed);
		const Eigen::MatrixXd& Cprime0 = es.eigenvectors();

		file << "\nInitial MO Coefficients:\n";

		OutputMatrix(Cprime0, file);

		Eigen::MatrixXd C = hartreeFock->V * Cprime0; // transform back the eigenvectors into the original non-orthogonalized AO basis

		Eigen::MatrixXd DensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				for (unsigned int vec = 0; vec < ((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied.size(); ++vec) // only eigenstates that are occupied 
					if (((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied[vec]) DensityMatrix(i, j) += C(i, vec) * C(j, vec); // 2 is for the number of electrons in the eigenstate, it's the restricted Hartree-Fock


		file << "\nInitial Density Matrix (no occupation multiplication):\n";

		OutputMatrix(DensityMatrix, file);

		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->eigenvals = es.eigenvalues();

		DensityMatrix *= 2;
		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->CalculateEnergy(((HartreeFock::RestrictedHartreeFock*)hartreeFock)->eigenvals, DensityMatrix/*, FockMatrix*/);

		file.precision(14);
		file << "\nInitial Hartree-Fock electronic energy: " << hartreeFock->GetTotalEnergy() - hartreeFock->nuclearRepulsionEnergy << " Hartrees" << std::endl;

		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->DensityMatrix = DensityMatrix;


		TestIterationsAndPostHF(molecule, hartreeFock, h, FockMatrix, FockTransformed, C, DensityMatrix, file, useDIIS);
	}

	delete hartreeFock;
}


void Test::Iterate(double& rmsD, double& deltaE, Systems::Molecule& molecule, HartreeFock::HartreeFockAlgorithm* hartreeFock, const Eigen::MatrixXd& h, Eigen::MatrixXd& FockMatrix, Eigen::MatrixXd& FockTransformed, Eigen::MatrixXd& C, Eigen::MatrixXd& DensityMatrix, std::ofstream& file, bool useDIIS)
{
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

		bool usedDiis = false;
		if (useDIIS) usedDiis = ((HartreeFock::RestrictedHartreeFock*)hartreeFock)->DIISStep(step, FockMatrix);

		FockTransformed = hartreeFock->Vt * FockMatrix * hartreeFock->V; // orthogonalize
		const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esl(FockTransformed);
		const Eigen::MatrixXd& Cprime = esl.eigenvectors();
		C = hartreeFock->V * Cprime; // transform back the eigenvectors into the original non-orthogonalized AO basis
		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->C = C;

		DensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				for (unsigned int vec = 0; vec < ((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied.size(); ++vec) // only eigenstates that are occupied 
					if (((HartreeFock::RestrictedHartreeFock*)hartreeFock)->occupied[vec]) DensityMatrix(i, j) += 2. * C(i, vec) * C(j, vec); // 2 is for the number of electrons in the eigenstate, it's the restricted Hartree-Fock


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
		file << step << "\t" << totalEnergy - hartreeFock->nuclearRepulsionEnergy << "\t" << totalEnergy << "\t" << deltaE << "\t" << rmsD;
		if (useDIIS) file << "\t" << hartreeFock->lastErrorEst;
		file << std::endl;


		if (rmsD < rmsDConvergence && (abs(deltaE) < (useDIIS ? energyConvergenceDIIS : energyConvergence)) && (!useDIIS || hartreeFock->lastErrorEst < diisConvergence)) break;

		//if (step == 36) break;
		// ***************************************************************************************************
		// go to the next density matrix

		((HartreeFock::RestrictedHartreeFock*)hartreeFock)->DensityMatrix = DensityMatrix;
	}

	file << std::endl;
}


void Test::TestIterationsAndPostHF(Systems::Molecule& molecule, HartreeFock::HartreeFockAlgorithm* hartreeFock, const Eigen::MatrixXd& h, Eigen::MatrixXd& FockMatrix, Eigen::MatrixXd& FockTransformed, Eigen::MatrixXd& C, Eigen::MatrixXd& DensityMatrix, std::ofstream& file, bool useDIIS)
{
	double rmsD = 0;
	double deltaE = 0;

	Iterate(rmsD, deltaE, molecule, hartreeFock, h, FockMatrix, FockTransformed, C, DensityMatrix, file, useDIIS);


	for (int atom = 0; atom < molecule.atoms.size(); ++atom)
		file << "Atom " << atom << " charge: " << hartreeFock->CalculateAtomicCharge(atom) << std::endl;

	file << std::endl;

	const double Emp2 = hartreeFock->CalculateMp2Energy();

	file << "Emp2: " << Emp2 << std::endl;
	file << "Total: " << hartreeFock->GetTotalEnergy() + Emp2 << std::endl;

	((HartreeFock::RestrictedCCSD*)hartreeFock)->InitCC();
	const double CCMP2 = ((HartreeFock::RestrictedCCSD*)hartreeFock)->MP2EnergyFromt4();
	file << "Emp2 from CC: " << CCMP2 << std::endl << std::endl;

	//file << "Emp2 from CC with general formula: " << ((HartreeFock::RestrictedCCSD*)hartreeFock)->CorrelationEnergy() << std::endl; 

	for (int i = 0; i < 100; ++i)
	{
		const double oldCCEnergy = ((HartreeFock::RestrictedCCSD*)hartreeFock)->GetCCEnergy();
		rmsD = ((HartreeFock::RestrictedCCSD*)hartreeFock)->StepCC(i);
		const double newCCEnergy = ((HartreeFock::RestrictedCCSD*)hartreeFock)->GetCCEnergy();

		file.precision(12);
		file << "Iter: " << i + 1 << "\tEcc = " << newCCEnergy << std::endl;


		deltaE = oldCCEnergy - newCCEnergy;

		if (rmsD < rmsDConvergence && (abs(deltaE) < (useDIIS ? CCenergyConvergenceDIIS : energyConvergence)) && (!useDIIS || hartreeFock->lastErrorEst < CCdiisConvergence)) break;
	}

	file << "Total CC: " << hartreeFock->GetTotalEnergy() + ((HartreeFock::RestrictedCCSD*)hartreeFock)->GetCCEnergy() << std::endl;

	const double Te = ((HartreeFock::RestrictedCCSD*)hartreeFock)->TEnergy();

	file << "E(T): " << Te << std::endl;

	file << "Total ECCSD(T): " << hartreeFock->GetTotalEnergy() + ((HartreeFock::RestrictedCCSD*)hartreeFock)->GetCCEnergy() + Te << std::endl;

	Vector3D moment(hartreeFock->GetMoment()); // multiply with 2.541746473 for Debye

	file << std::endl;
	file << "Mu-x: " << moment.X << " au, " << moment.X * 2.541746473 << " D" << std::endl;
	file << "Mu-y: " << moment.Y << " au, " << moment.Y * 2.541746473 << " D" << std::endl;
	file << "Mu-z: " << moment.Z << " au, " << moment.Z * 2.541746473 << " D" << std::endl;

	const double total = moment.Length();
	file << "Total dipole moment: " << total << " au, " << total * 2.541746473 << " D" << std::endl;


	// test CIS

	HartreeFock::RestrictedConfigurationIInteractionSingles restrictedCIS((HartreeFock::RestrictedHartreeFock*)hartreeFock);

	restrictedCIS.Init();

	Eigen::MatrixXd CISH = restrictedCIS.getSpinOrbitalCISMatrix();

	file.precision(5);
	file << "\nCIS Hamiltonian Matrix: \n";

	OutputMatrix(CISH, file);


	const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> escis(CISH, Eigen::DecompositionOptions::EigenvaluesOnly);

	Eigen::VectorXd eigenvals = escis.eigenvalues();

	file << "\nCIS Eigenvalues: \n";

	for (int i = 0; i < eigenvals.rows(); ++i)
		file << eigenvals(i) << std::endl;


	file.precision(5);
	file << "\nSpin-adapted CIS singlet Hamiltonian Matrix: \n";

	CISH = restrictedCIS.getSpinAdaptedCISSinglet();

	OutputMatrix(CISH, file);


	file << "\nCIS Singlet Eigenvalues: \n";

	const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> escisSinglet(CISH, Eigen::DecompositionOptions::EigenvaluesOnly);

	eigenvals = escisSinglet.eigenvalues();

	for (int i = 0; i < eigenvals.rows(); ++i)
		file << eigenvals(i) << std::endl;



	file.precision(5);
	file << "\nSpin-adapted CIS triplet Hamiltonian Matrix: \n";

	CISH = restrictedCIS.getSpinAdaptedCISTriplet();

	OutputMatrix(CISH, file);


	file << "\nCIS Triplet Eigenvalues: \n";

	const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> escisTriplet(CISH, Eigen::DecompositionOptions::EigenvaluesOnly);

	eigenvals = escisTriplet.eigenvalues();

	for (int i = 0; i < eigenvals.rows(); ++i)
		file << eigenvals(i) << std::endl;

	const Eigen::MatrixXd TDHFm = restrictedCIS.getTDHFMatrix();

	file << "\nTDHF / RPA Hamiltonian: \n";
	OutputMatrix(TDHFm, file);

	const Eigen::EigenSolver<Eigen::MatrixXd> TDHFSolver(TDHFm, false);

	const Eigen::VectorXcd& eigenv = TDHFSolver.eigenvalues();
	std::vector<double> vals;
	for (int i = 0; i < eigenv.rows(); ++i)
		vals.push_back(eigenv(i).real());

	std::sort(vals.begin(), vals.end());

	file << "\nRPA Excitation Energies: \n";

	for (int i = 0; i < vals.size(); ++i)
		file << vals[i] << std::endl;


	file << "\nRPA Excitation Energies (method 2): \n";

	const Eigen::MatrixXd TDHFmr = restrictedCIS.getReducedTDHFMatrix();
	const Eigen::EigenSolver<Eigen::MatrixXd> TDHFReducedSolver(TDHFmr, false);

	const Eigen::VectorXcd& eigenv2 = TDHFReducedSolver.eigenvalues();

	vals.clear();
	for (int i = 0; i < eigenv2.rows(); ++i)
		vals.push_back(sqrt(eigenv2(i).real()));
	std::sort(vals.begin(), vals.end());

	for (int i = 0; i < vals.size(); ++i)
		file << vals[i] << std::endl;
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



void Test::InitWaterMolecule(Systems::Molecule& molecule)
{
	Systems::AtomWithShells H1, H2, O;

	for (auto& atom : basis.atoms)
	{
		if (1 == atom.Z)
		{
			H1 = H2 = atom;
		}
		else if (8 == atom.Z)
			O = atom;
	}

	O.position.X = 0;
	O.position.Y = -0.143225816552;
	O.position.Z = 0;

	H1.position.X = 1.638036840407;
	H1.position.Y = 1.136548822547;
	H1.position.Z = 0;

	H2.position.X = -1.638036840407;
	H2.position.Y = 1.136548822547;
	H2.position.Z = 0;

	// translate the molecule with O in zero
	//const Vector3D<double> translate = O.position;
	//O.position = Vector3D<double>(0, 0, 0);
	//H1.position -= translate;
	//H2.position -= translate;

	molecule.atoms.push_back(O);
	molecule.atoms.push_back(H1);
	molecule.atoms.push_back(H2);
	molecule.Init();
}



// TODO: try to load the provided file at https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2303 and do comparisons in the code
// even the geometry of the molecule could be loaded from the file

void Test::TestWater(const std::string& fileName, const std::string& sfileName, const std::string& tfileName, const std::string& vfileName, const std::string& erifileName, bool useDIIS)
{
	Systems::Molecule molecule;

	InitWaterMolecule(molecule);

	std::ofstream file(fileName);

	OutputMatrices(molecule, file, sfileName, tfileName, vfileName, erifileName, useDIIS, 8.002367061810450);
}

void Test::TestWaterDipoleMoment(const std::string& fileName)
{
	Systems::Molecule molecule;

	InitWaterMolecule(molecule);

	std::ofstream file(fileName);

	HartreeFock::HartreeFockAlgorithm* hartreeFock;
	const bool restricted = molecule.alphaElectrons == molecule.betaElectrons;
	if (restricted)
		hartreeFock = new HartreeFock::RestrictedHartreeFock();
		//hartreeFock = new HartreeFock::RestrictedCCSD();
	else
	{
		hartreeFock = new HartreeFock::UnrestrictedHartreeFock();
		((HartreeFock::UnrestrictedHartreeFock*)hartreeFock)->addAsymmetry = false;
	}

	hartreeFock->UseDIIS = false;
	hartreeFock->alpha = 0.5;
	hartreeFock->initGuess = 0.7;

	hartreeFock->Init(&molecule);
	double E = hartreeFock->Calculate();
	Vector3D scfMoment(hartreeFock->GetMoment());

	file << "Scf Mu-x: " << scfMoment.X << " au, " << scfMoment.X * 2.541746473 << " D" << std::endl;
	file << "Scf Mu-y: " << scfMoment.Y << " au, " << scfMoment.Y * 2.541746473 << " D" << std::endl;
	file << "Scf Mu-z: " << scfMoment.Z << " au, " << scfMoment.Z * 2.541746473 << " D" << std::endl;

	double total = scfMoment.Length();
	file << "Scf Total dipole moment: " << total << " au, " << total * 2.541746473 << " D" << std::endl << std::endl;

	// the above method gives basically the same value as the following one (with numerical derivatives), if MP2 energy is not included
	// that is, for STO3G basis and H2O molecule, the above gives 0.603521 au, the following gives 0.60352

	const double deltaE = 0.001;

	const double deltaE2 = 2 * deltaE;

	molecule.ElectricField.X = -deltaE;
	hartreeFock->Init(&molecule);
	double Ex1 = hartreeFock->Calculate();
	Ex1 += hartreeFock->CalculateMp2Energy();

	if (!hartreeFock->converged) file << "Not converged Ex1" << std::endl;


	molecule.ElectricField.X = deltaE;
	hartreeFock->Init(&molecule);
	double Ex2 = hartreeFock->Calculate();
	Ex2 += hartreeFock->CalculateMp2Energy();

	if (!hartreeFock->converged) file << "Not converged Ex2" << std::endl;


	molecule.ElectricField.X = 0;


	molecule.ElectricField.Y = -deltaE;
	hartreeFock->Init(&molecule);
	double Ey1 = hartreeFock->Calculate();
	Ey1 += hartreeFock->CalculateMp2Energy();
	
	if (!hartreeFock->converged) file << "Not converged Ey1" << std::endl;


	molecule.ElectricField.Y = deltaE;
	hartreeFock->Init(&molecule);
	double Ey2 = hartreeFock->Calculate();
	Ey2 += hartreeFock->CalculateMp2Energy();

	if (!hartreeFock->converged) file << "Not converged Ey2" << std::endl;

	molecule.ElectricField.Y = 0;

	molecule.ElectricField.Z = -deltaE;
	hartreeFock->Init(&molecule);
	double Ez1 = hartreeFock->Calculate();
	Ez1 += hartreeFock->CalculateMp2Energy();

	if (!hartreeFock->converged) file << "Not converged Ez1" << std::endl;

	molecule.ElectricField.Z = deltaE;
	hartreeFock->Init(&molecule);
	double Ez2 = hartreeFock->Calculate();
	Ez2 += hartreeFock->CalculateMp2Energy();

	if (!hartreeFock->converged) file << "Not converged Ez2" << std::endl;

	Vector3D mu(hartreeFock->GetNuclearMoment());

	delete hartreeFock;

	mu.X += (Ex1 - Ex2) / deltaE2;
	mu.Y += (Ey1 - Ey2) / deltaE2;
	mu.Z += (Ez1 - Ez2) / deltaE2;

	file << "Mu-x: " << mu.X << " au, " << mu.X * 2.541746473 << " D" << std::endl;
	file << "Mu-y: " << mu.Y << " au, " << mu.Y * 2.541746473 << " D" << std::endl;
	file << "Mu-z: " << mu.Z << " au, " << mu.Z * 2.541746473 << " D" << std::endl;

	total = mu.Length();
	file << "Total dipole moment: " << total << " au, " << total * 2.541746473 << " D" << std::endl;

	file << std::endl;
}


void Test::TestMethane(const std::string& fileName, const std::string& sfileName, const std::string& tfileName, const std::string& vfileName, const std::string& erifileName, bool useDIIS)
{
	Systems::AtomWithShells H1, H2, H3, H4, C;

	for (auto& atom : basis.atoms)
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

	OutputMatrices(molecule, file, sfileName, tfileName, vfileName, erifileName, useDIIS, 13.497304462036480);
}
