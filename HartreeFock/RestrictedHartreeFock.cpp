#include "stdafx.h"
#include "RestrictedHartreeFock.h"

namespace HartreeFock {


	RestrictedHartreeFock::RestrictedHartreeFock(int iterations)
		: HartreeFockAlgorithm(iterations), totalEnergy(std::numeric_limits<double>::infinity()), nrOccupiedLevels(0)
	{
	}


	RestrictedHartreeFock::~RestrictedHartreeFock()
	{
	}


	void RestrictedHartreeFock::Init(Systems::Molecule* molecule)
	{
		HartreeFockAlgorithm::Init(molecule);

		DensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		nrOccupiedLevels = molecule->ElectronsNumber() / 2;

		// this should not happen for a well behaving molecule for the restricted method, but someone might try calculate the Hydrogen atom with it
		// prevent a crash with this
		if (!nrOccupiedLevels) nrOccupiedLevels = 1;

		if (nrOccupiedLevels > static_cast<unsigned int>(numberOfOrbitals)) nrOccupiedLevels = numberOfOrbitals;
	}

	void RestrictedHartreeFock::Step(int iter)
	{
		// *****************************************************************************************************************

		// the Fock matrix
		Eigen::MatrixXd FockMatrix;
		
		InitFockMatrix(iter, FockMatrix);

		// ***************************************************************************************************************************

		// solve the Roothaan-Hall equation

		//diagonalize it - it can be done faster by diagonalizing the overlap matrix outside the loop but for tests this should be enough
		//I leave it here just in case - if all that S, U, s, V seems confusing this should help :)
		//Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(FockMatrix, overlapMatrix.matrix);
		//const Eigen::MatrixXd& C = es.eigenvectors();


		// this hopefully is faster than the one commented above
		
		Eigen::MatrixXd FockTransformed = Vt * FockMatrix * V;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(FockTransformed);
		const Eigen::MatrixXd& Cprime = es.eigenvectors();
		Eigen::MatrixXd C = V * Cprime;

		// normalize it
		//NormalizeC(C, nrOccupiedLevels);
		
		//***************************************************************************************************************

		// calculate the density matrix

		Eigen::MatrixXd newDensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				for (unsigned int vec = 0; vec < nrOccupiedLevels; ++vec) // only eigenstates that are occupied 
					newDensityMatrix(i, j) += 2. * C(i, vec) * C(j, vec); // 2 is for the number of electrons in the eigenstate, it's the restricted Hartree-Fock
															  
		
		//**************************************************************************************************************

		const Eigen::VectorXd& eigenvals = es.eigenvalues();

		CalculateEnergy(eigenvals, newDensityMatrix/*, FockMatrix*/);

		TRACE("Step: %d Energy: %f\n", iter, totalEnergy);
															  
		// ***************************************************************************************************
		// go to the next density matrix

		DensityMatrix = alpha * newDensityMatrix + (1. - alpha) * DensityMatrix;  // use mixing if alpha is set less than 1
	}

	void RestrictedHartreeFock::InitFockMatrix(int iter, Eigen::MatrixXd& FockMatrix) const
	{
		// this could be made faster knowing that the matrix should be symmetric
		// but it would be less expressive so I'll let it as it is
		// maybe I'll improve it later
		// anyway, the slower part is dealing with electron-electron integrals

		if (0 == iter)
		{
			if (initGuess > 0)
			{
				FockMatrix.resize(h.rows(), h.cols());

				for (int i = 0; i < h.rows(); ++i)
					for (int j = 0; j < h.rows(); ++j)
						FockMatrix(i, j) = initGuess * overlapMatrix.matrix(i, j) * (h(i, i) + h(j, j)) / 2.;
			}
			else FockMatrix = h;
		}
		else
		{
			Eigen::MatrixXd G = Eigen::MatrixXd::Zero(h.rows(), h.cols());

			for (int i = 0; i < numberOfOrbitals; ++i)
				for (int j = 0; j < numberOfOrbitals; ++j)
					for (int k = 0; k < numberOfOrbitals; ++k)
						for (int l = 0; l < numberOfOrbitals; ++l)
						{
							double coulomb = integralsRepository.getElectronElectron(i, j, k, l);
							double exchange = integralsRepository.getElectronElectron(i, l, k, j);

							G(i, j) += DensityMatrix(k, l) * (coulomb - 1. / 2. * exchange);
						}

			FockMatrix = h + G;
		}
	}

	void RestrictedHartreeFock::CalculateEnergy(const Eigen::VectorXd& eigenvals, const Eigen::MatrixXd& calcDensityMatrix/*, Eigen::MatrixXd& F*/)
	{
		// one way of calculating the energy

		totalEnergy = 0;

		// simply adding energy levels does not work, because of over counting interaction energy

		for (int i = 1; i < h.rows(); ++i)
			for (int j = 0; j < i; ++j)
				totalEnergy += calcDensityMatrix(i, j) * h(j, i);

		// only the values below the diagonal were added
		// *2 to have the sum of all elements except those on diagonal
		totalEnergy *= 2.; 
		
		// now add the diagonal elements, too
		for (int i = 0; i < h.rows(); ++i) totalEnergy += calcDensityMatrix(i, i) * h(i, i);

		totalEnergy /= 2.;
	
		for (unsigned int level = 0; level < nrOccupiedLevels; ++level)
			totalEnergy += eigenvals(level);

		HOMOEnergy = eigenvals(nrOccupiedLevels - 1);

		

		// the above is equivalent with this commented code, but since we already have eigenvalues for F calculated, we can get rid of the matrix addition to obtain H
/*
		Eigen::MatrixXd H = h + F;

		for (int i = 0; i < H.rows(); ++i)
			for (int j = 0; j < H.cols(); ++j)
				totalEnergy += calcP(i, j) * H(j, i);

		totalEnergy /= 2.;
*/

		totalEnergy += nuclearRepulsionEnergy;
	}

	double RestrictedHartreeFock::GetTotalEnergy() const
	{
		return totalEnergy;
	}
}