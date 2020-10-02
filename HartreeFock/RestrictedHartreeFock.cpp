#include "stdafx.h"
#include "RestrictedHartreeFock.h"

namespace HartreeFock {


	RestrictedHartreeFock::RestrictedHartreeFock(int iterations)
		: HartreeFockAlgorithm(iterations), nrOccupiedLevels(0)
	{
	}


	RestrictedHartreeFock::~RestrictedHartreeFock()
	{
	}


	void RestrictedHartreeFock::Init(Systems::Molecule* molecule)
	{
		HartreeFockAlgorithm::Init(molecule);

		DensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		occupied.resize(0); // just in case it was resized before
		nrOccupiedLevels = molecule->ElectronsNumber() / 2;

		// this should not happen for a well behaving molecule for the restricted method, but someone might try calculate the Hydrogen atom with it
		// prevent a crash with this
		if (!nrOccupiedLevels) nrOccupiedLevels = 1;
		if (nrOccupiedLevels > static_cast<unsigned int>(numberOfOrbitals)) nrOccupiedLevels = numberOfOrbitals;

		occupied.resize(nrOccupiedLevels, true);
	}

	double RestrictedHartreeFock::Step(int iter)
	{
		// *****************************************************************************************************************

		// the Fock matrix
		Eigen::MatrixXd FockMatrix;

		InitFockMatrix(iter, FockMatrix);

		// will be used for DIIS

		bool UsedDIIS = false;

		if (UseDIIS && iter && iter < maxDIISiterations)
		{
			// the density matrix should commute with the Fock matrix. The difference is the error.
			const Eigen::MatrixXd errorMatrix = overlapMatrix.matrix * DensityMatrix * FockMatrix - FockMatrix * DensityMatrix * overlapMatrix.matrix;

			// another choice: if reached the limit of kept matrices, replace the ones with the bigger error
			errorMatrices.emplace_back(errorMatrix);
			fockMatrices.push_back(FockMatrix);
			if (errorMatrices.size() > 8)
			{
				errorMatrices.pop_front();
				fockMatrices.pop_front();
			}

			if (errorMatrices.size() > 6)
			{
				// use DIIS
				const size_t nrMatrices = errorMatrices.size();
				Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nrMatrices + 1, nrMatrices + 1);

				auto errorIter1 = errorMatrices.begin();
				for (size_t i = 0; i < nrMatrices; ++i)
				{
					auto errorIter2 = errorMatrices.begin();

					for (size_t j = 0; j < i; ++j)
					{
						B(i, j) = B(j, i) = (*errorIter1).cwiseProduct(*errorIter2).sum();

						++errorIter2;
					}

					B(i, i) = (*errorIter1).cwiseProduct(*errorIter1).sum();

					if (i == nrMatrices - 1) lastErrorEst = sqrt(B(i, i));

					B(nrMatrices, i) = B(i, nrMatrices) = 1;

					++errorIter1;
				}

				// Solve the system of linear equations

				Eigen::VectorXd C = Eigen::VectorXd::Zero(nrMatrices + 1);
				C(nrMatrices) = 1;

				//C = B.fullPivLu().solve(C);
				C = B.colPivHouseholderQr().solve(C);
				//C = B.fullPivHouseholderQr().solve(C);

				// compute the new Fock matrix

				FockMatrix = Eigen::MatrixXd::Zero(FockMatrix.rows(), FockMatrix.cols());

				auto fockIter = fockMatrices.begin();
				for (size_t i = 0; i < nrMatrices; ++i)
				{
					FockMatrix += C(i) * *fockIter;
					++fockIter;
				}

				UsedDIIS = true;
			}
		}
		else lastErrorEst = 0;


		// ***************************************************************************************************************************

		// solve the Roothaan-Hall equation

		//diagonalize it - it can be done faster by diagonalizing the overlap matrix outside the loop but for tests this should be enough
		//I leave it here just in case - if all that S, U, s, V seems confusing this should help :)
		//Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(FockMatrix, overlapMatrix.matrix);
		//const Eigen::MatrixXd& C = es.eigenvectors();


		// this hopefully is faster than the one commented above

		const Eigen::MatrixXd FockTransformed = Vt * FockMatrix * V; // orthogonalize
		const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(FockTransformed);
		const Eigen::MatrixXd& Cprime = es.eigenvectors();
		Ce = V * Cprime; // transform back the eigenvectors into the original non-orthogonalized AO basis

		// normalize it
		//NormalizeC(C, occupied);

		//***************************************************************************************************************

		// calculate the density matrix

		Eigen::MatrixXd newDensityMatrix = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				for (unsigned int vec = 0; vec < occupied.size(); ++vec) // only eigenstates that are occupied 
					if (occupied[vec]) newDensityMatrix(i, j) += 2. * Ce(i, vec) * Ce(j, vec); // 2 is for the number of electrons in the eigenstate, it's the restricted Hartree-Fock


		//**************************************************************************************************************

		eigenvals = es.eigenvalues();

		CalculateEnergy(eigenvals, newDensityMatrix/*, FockMatrix*/);

		TRACE("Step: %d Energy: %f\n", iter, totalEnergy);

		// calculate rms for differences between new and old density matrices, it can be used to check for convergence, too
		Eigen::MatrixXd rmsDensityMatricesDif = newDensityMatrix - DensityMatrix;
		rmsDensityMatricesDif = rmsDensityMatricesDif.cwiseProduct(rmsDensityMatricesDif);
		const double rmsD = sqrt(rmsDensityMatricesDif.sum());

		// ***************************************************************************************************
		// go to the next density matrix

		if (UsedDIIS)
			DensityMatrix = newDensityMatrix;
		else
			DensityMatrix = alpha * newDensityMatrix + (1. - alpha) * DensityMatrix;  // use mixing if alpha is set less than 1

		return rmsD;
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
						FockMatrix(i, j) = initGuess * overlapMatrix.matrix(i, j) * (h(i, i) + h(j, j)) * 0.5;
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
							const double coulomb = integralsRepository.getElectronElectron(i, j, k, l);
							const double exchange = integralsRepository.getElectronElectron(i, l, k, j);

							G(i, j) += DensityMatrix(k, l) * (coulomb - 0.5 * exchange);
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

		totalEnergy *= 0.5;

		for (unsigned int level = 0; level < occupied.size(); ++level)
			if (occupied[level]) totalEnergy += eigenvals(level);

		HOMOEnergy = eigenvals(nrOccupiedLevels - 1);



		// the above is equivalent with this commented code, but since we already have eigenvalues for F calculated, we can get rid of the matrix addition to obtain H
/*
		Eigen::MatrixXd H = h + F;

		for (int i = 0; i < H.rows(); ++i)
			for (int j = 0; j < H.cols(); ++j)
				totalEnergy += calcDensityMatrix(i, j) * H(j, i);

		totalEnergy /= 2.;
*/

		totalEnergy += nuclearRepulsionEnergy;
	}

	double RestrictedHartreeFock::CalculateMp2Energy()
	{
		mp2Energy = 0;

		// TODO: calculate it

		return mp2Energy;
	}

}