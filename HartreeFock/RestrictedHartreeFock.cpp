#include "stdafx.h"
#include "RestrictedHartreeFock.h"

namespace HartreeFock {


	RestrictedHartreeFock::RestrictedHartreeFock(int iterations)
		: HartreeFockAlgorithm(iterations), totalEnergy(std::numeric_limits<double>::infinity()), nrLevels(0)
	{
	}


	RestrictedHartreeFock::~RestrictedHartreeFock()
	{
	}


	void RestrictedHartreeFock::Init(Systems::Molecule* molecule)
	{
		HartreeFockAlgorithm::Init(molecule);

		P = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		nrLevels = molecule->ElectronsNumber() / 2;

		if (nrLevels > (unsigned int)numberOfOrbitals) nrLevels = numberOfOrbitals;
	}

	void RestrictedHartreeFock::Step(int iter)
	{
		// *****************************************************************************************************************

		// the Fock matrix
		Eigen::MatrixXd F;
		
		InitFockMatrix(iter, F);

		// ***************************************************************************************************************************

		// solve the Roothaan-Hall equation

		//diagonalize it - it can be done faster by diagonalizing the overlap matrix outside the loop but for tests this should be enough
		//I leave it here just in case - if all that S, U, s, V seems confusing this should help :)
		//Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(F, overlapMatrix.matrix);
		//const Eigen::MatrixXd& C = es.eigenvectors();


		// this hopefully is faster than the one commented above
		
		Eigen::MatrixXd Fprime = Vt * F * V;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Fprime);
		const Eigen::MatrixXd& Cprime = es.eigenvectors();
		Eigen::MatrixXd C = V * Cprime;
		
		//***************************************************************************************************************

		// calculate the density matrix

		Eigen::MatrixXd newP = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				for (unsigned int vec = 0; vec < nrLevels; ++vec) // only eigenstates that are occupied 
					newP(i, j) += 2. * C(i, vec) * C(j, vec); // 2 is for the number of electrons in the eigenstate, it's the restricted Hartree-Fock
															  
		
		//**************************************************************************************************************

		const Eigen::VectorXd& eigenvals = es.eigenvalues();

		CalculateEnergy(eigenvals, newP/*, F*/);

		TRACE("Step: %d Energy: %f\n", iter, totalEnergy);
															  
		// ***************************************************************************************************
		// go to the next density matrix

		P = alpha * newP + (1. - alpha) * P;  // use mixing is alpha is set less than 1
	}

	void RestrictedHartreeFock::InitFockMatrix(int iter, Eigen::MatrixXd& F) const
	{
		// this could be made twice as fast knowing that the matrix should be symmetric
		// but it would be less expressive so I'll let it as it is

		if (0 == iter)
		{
			if (initGuess > 0)
			{
				F.resize(h.rows(), h.cols());

				for (int i = 0; i < h.rows(); ++i)
					for (int j = 0; j < h.rows(); ++j)
						F(i, j) = initGuess * overlapMatrix.matrix(i, j) * (h(i, i) + h(j, j)) / 2.;
			}
			else F = h;
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

							G(i, j) += P(k, l) * (coulomb - 1. / 2. * exchange);
						}

			F = h + G;
		}
	}

	void RestrictedHartreeFock::CalculateEnergy(const Eigen::VectorXd& eigenvals, const Eigen::MatrixXd& calcP/*, Eigen::MatrixXd& F*/)
	{
		// one way of calculating the energy

		totalEnergy = 0;

		// simply adding energy levels does not work, because of over counting interaction energy
		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				totalEnergy += calcP(i, j) * h(j, i);

		totalEnergy /= 2.;

		for (unsigned int level = 0; level < nrLevels; ++level)
			totalEnergy += eigenvals(level);

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