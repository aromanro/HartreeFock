#include "stdafx.h"
#include "UnrestrictedHartreeFock.h"


namespace HartreeFock {


	UnrestrictedHartreeFock::UnrestrictedHartreeFock(int iterations)
		: HartreeFockAlgorithm(iterations), nrOccupiedLevelsPlus(0), nrOccupiedLevelsMinus(0), addAsymmetry(true), asymmetry(0.1)
	{
	}


	UnrestrictedHartreeFock::~UnrestrictedHartreeFock()
	{
	}

	void UnrestrictedHartreeFock::Init(Systems::Molecule* molecule)
	{
		HartreeFockAlgorithm::Init(molecule);


		DensityMatrixPlus = Eigen::MatrixXd::Zero(h.rows(), h.cols());
		DensityMatrixMinus = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		occupiedPlus.resize(0);
		occupiedMinus.resize(0);

		const unsigned int electronsNumber = molecule->ElectronsNumber();
		nrOccupiedLevelsMinus = static_cast<unsigned int>(floor(electronsNumber / 2.));
		nrOccupiedLevelsPlus = electronsNumber - nrOccupiedLevelsMinus;

		if (molecule->alphaElectrons > 0 || molecule->betaElectrons > 0)
		{
			nrOccupiedLevelsPlus = molecule->alphaElectrons;
			nrOccupiedLevelsMinus = molecule->betaElectrons;
		}


		// preventing setting too many electrons
		if (nrOccupiedLevelsMinus + nrOccupiedLevelsPlus > static_cast<unsigned int>(2 * numberOfOrbitals))
		{
			nrOccupiedLevelsMinus = numberOfOrbitals;
			nrOccupiedLevelsPlus = numberOfOrbitals;
		}

		if (nrOccupiedLevelsPlus > static_cast<unsigned int>(numberOfOrbitals))
		{
			const int dif = nrOccupiedLevelsPlus - numberOfOrbitals;
			nrOccupiedLevelsPlus = numberOfOrbitals;
			nrOccupiedLevelsMinus += dif;
		}

		if (nrOccupiedLevelsMinus > static_cast<unsigned int>(numberOfOrbitals))
			nrOccupiedLevelsMinus = numberOfOrbitals;

		occupiedPlus.resize(nrOccupiedLevelsPlus, true);
		occupiedMinus.resize(nrOccupiedLevelsMinus, true);
	}

	double UnrestrictedHartreeFock::Step(int iter)
	{
		// *****************************************************************************************************************

		// the Fock matrices
		Eigen::MatrixXd FockMatrixPlus;
		Eigen::MatrixXd FockMatrixMinus;

		InitFockMatrices(iter, FockMatrixPlus, FockMatrixMinus);

		// will be used for DIIS

		
		bool UsedDIIS = false;
		 
		if (UseDIIS && iter && iter < maxDIISiterations)
		{
			const Eigen::MatrixXd errorMatrixPlus = overlapMatrix.matrix * DensityMatrixPlus * FockMatrixPlus - FockMatrixPlus * DensityMatrixPlus * overlapMatrix.matrix;

			errorMatricesPlus.emplace_back(errorMatrixPlus);
			fockMatricesPlus.push_back(FockMatrixPlus);

			const Eigen::MatrixXd errorMatrixMinus = overlapMatrix.matrix * DensityMatrixMinus * FockMatrixMinus - FockMatrixMinus * DensityMatrixMinus * overlapMatrix.matrix;

			errorMatricesMinus.emplace_back(errorMatrixMinus);
			fockMatricesMinus.push_back(FockMatrixMinus);

			if (errorMatricesPlus.size() > 8)
			{
				errorMatricesPlus.pop_front();
				fockMatricesPlus.pop_front();
				errorMatricesMinus.pop_front();
				fockMatricesMinus.pop_front();
			}

			if (errorMatricesPlus.size() > 3)
			{
				// use DIIS
				const size_t nrMatrices = errorMatricesPlus.size();
				Eigen::MatrixXd Bplus = Eigen::MatrixXd::Zero(nrMatrices + 1, nrMatrices + 1);
				Eigen::MatrixXd Bminus = Eigen::MatrixXd::Zero(nrMatrices + 1, nrMatrices + 1);

				auto errorPlusIter1 = errorMatricesPlus.begin();
				auto errorMinusIter1 = errorMatricesMinus.begin();
				for (size_t i = 0; i < nrMatrices; ++i)
				{
					auto errorPlusIter2 = errorMatricesPlus.begin();
					auto errorMinusIter2 = errorMatricesMinus.begin();

					for (size_t j = 0; j < i; ++j)
					{
						Bplus(i, j) = Bplus(j, i) = (*errorPlusIter1).cwiseProduct(*errorPlusIter2).sum();
						Bminus(i, j) = Bminus(j, i) = (*errorMinusIter1).cwiseProduct(*errorMinusIter2).sum();

						++errorPlusIter2;
						++errorMinusIter2;
					}

					Bplus(i, i) = (*errorPlusIter1).cwiseProduct(*errorPlusIter1).sum();
					Bminus(i, i) = (*errorMinusIter1).cwiseProduct(*errorMinusIter1).sum();

					if (i == nrMatrices - 1) lastErrorEst = sqrt(Bplus(i, i) + Bminus(i, i));

					Bplus(nrMatrices, i) = Bplus(i, nrMatrices) = 1;
					Bminus(nrMatrices, i) = Bminus(i, nrMatrices) = 1;

					++errorPlusIter1;
					++errorMinusIter1;
				}

				// Solve the systems of linear equations

				Eigen::VectorXd CPlus = Eigen::VectorXd::Zero(nrMatrices + 1);
				CPlus(nrMatrices) = 1;
				Eigen::VectorXd CMinus = Eigen::VectorXd::Zero(nrMatrices + 1);
				CMinus(nrMatrices) = 1;

				//CPlus = Bplus.fullPivHouseholderQr().solve(CPlus);
				//CMinus = Bminus.fullPivHouseholderQr().solve(CMinus);
				CPlus = Bplus.colPivHouseholderQr().solve(CPlus);
				CMinus = Bminus.colPivHouseholderQr().solve(CMinus);

				// compute the new Fock matrices

				FockMatrixPlus = Eigen::MatrixXd::Zero(FockMatrixPlus.rows(), FockMatrixPlus.cols());
				FockMatrixMinus = Eigen::MatrixXd::Zero(FockMatrixMinus.rows(), FockMatrixMinus.cols());

				auto fockIterPlus = fockMatricesPlus.begin();
				auto fockIterMinus = fockMatricesMinus.begin();
				for (size_t i = 0; i < nrMatrices; ++i)
				{
					FockMatrixPlus += CPlus(i) * (*fockIterPlus);
					FockMatrixMinus += CMinus(i) * (*fockIterMinus);

					++fockIterPlus;
					++fockIterMinus;
				}

				UsedDIIS = true;
			}
		}
		else lastErrorEst = 0;
		

		// ***************************************************************************************************************************

		// solve the Pople-Nesbet–Berthier equations

		 // orthogonalize
		const Eigen::MatrixXd FockMatrixPlusTransformed = Vt * FockMatrixPlus * V;
		const Eigen::MatrixXd FockMatrixMinusTransformed = Vt * FockMatrixMinus * V;

		if (FockMatrixPlusTransformed.rows() > 1)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esplus(FockMatrixPlusTransformed);
			const Eigen::MatrixXd& Cplusprime = esplus.eigenvectors();
			Cplus = V * Cplusprime; // transform back the eigenvectors into the original non-orthogonalized AO basis
			eigenvalsplus = esplus.eigenvalues();
		}
		else
		{
			Cplus = V * Eigen::MatrixXd::Ones(1, 1);
			eigenvalsplus = FockMatrixPlusTransformed(0, 0) * Eigen::VectorXd::Ones(1);
		}

		if (FockMatrixMinusTransformed.rows() > 1)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esminus(FockMatrixMinusTransformed);
			const Eigen::MatrixXd& Cminusprime = esminus.eigenvectors();
			Cminus = V * Cminusprime; // transform back the eigenvectors into the original non-orthogonalized AO basis
			eigenvalsminus = esminus.eigenvalues();
		}
		else
		{
			Cminus = V * Eigen::MatrixXd::Ones(1, 1);
			eigenvalsminus = FockMatrixMinusTransformed(0, 0) * Eigen::VectorXd::Ones(1);
		}

		// normalize them
		//NormalizeC(Cplus, occupiedPlus);
		//NormalizeC(Cminus, occupiedMinus);

		//***************************************************************************************************************

		// calculate the density matrices

		Eigen::MatrixXd newDensityMatrixPlus = Eigen::MatrixXd::Zero(h.rows(), h.cols());
		Eigen::MatrixXd newDensityMatrixMinus = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
			{
				for (unsigned int vec = 0; vec < occupiedPlus.size(); ++vec) // only eigenstates that are occupied 
					if (occupiedPlus[vec]) newDensityMatrixPlus(i, j) += Cplus(i, vec) * Cplus(j, vec);

				for (unsigned int vec = 0; vec < occupiedMinus.size(); ++vec) // only eigenstates that are occupied 
					if (occupiedMinus[vec]) newDensityMatrixMinus(i, j) += Cminus(i, vec) * Cminus(j, vec);
			}

		//**************************************************************************************************************

		CalculateEnergy(eigenvalsplus, eigenvalsminus, newDensityMatrixPlus, newDensityMatrixMinus/*, Fplus, Fminus*/);

		TRACE("Step: %d Energy: %f\n", iter, totalEnergy);

		// calculate rms for differences between new and old density matrices, it can be used to check for convergence, too
		Eigen::MatrixXd rmsDensityMatricesDif = newDensityMatrixPlus - DensityMatrixPlus;
		rmsDensityMatricesDif = rmsDensityMatricesDif.cwiseProduct(rmsDensityMatricesDif);
		double rmsD = sqrt(rmsDensityMatricesDif.sum());

		rmsDensityMatricesDif = newDensityMatrixMinus - DensityMatrixMinus;
		rmsDensityMatricesDif = rmsDensityMatricesDif.cwiseProduct(rmsDensityMatricesDif);
		rmsD += sqrt(rmsDensityMatricesDif.sum());

		// ***************************************************************************************************
		// go to the next density matrices
		// use mixing if alpha is set less then one
		
		if (UsedDIIS)
		{
			DensityMatrixPlus = newDensityMatrixPlus;
			DensityMatrixMinus = newDensityMatrixMinus;
		}
		else
		{
			DensityMatrixPlus = alpha * newDensityMatrixPlus + (1. - alpha) * DensityMatrixPlus;
			DensityMatrixMinus = alpha * newDensityMatrixMinus + (1. - alpha) * DensityMatrixMinus;
		}

		return rmsD;
	}

	void UnrestrictedHartreeFock::InitFockMatrices(int iter, Eigen::MatrixXd& FockMatrixPlus, Eigen::MatrixXd& FockMatrixMinus) const
	{
		// this could be made faster knowing that the matrix should be symmetric
		// but it would be less expressive so I'll let it as it is
		// maybe I'll improve it later
		// anyway, the slower part is dealing with electron-electron integrals

		if (0 == iter)
		{
			if (initGuess > 0)
			{
				FockMatrixPlus.resize(h.rows(), h.cols());
				FockMatrixMinus.resize(h.rows(), h.cols());

				for (int i = 0; i < h.rows(); ++i)
					for (int j = 0; j < h.rows(); ++j)
						FockMatrixPlus(i, j) = FockMatrixMinus(i, j) = initGuess * overlapMatrix.matrix(i, j) * (h(i, i) + h(j, j)) / 2.;
			}
			else
			{
				FockMatrixPlus = h;
				FockMatrixMinus = h;
			}

			if (addAsymmetry && FockMatrixPlus.cols() > 1)
			{
				FockMatrixPlus(0, 1) += asymmetry;
				FockMatrixPlus(1, 0) = FockMatrixPlus(0, 1);
			}
		}
		else
		{
			Eigen::MatrixXd Gplus = Eigen::MatrixXd::Zero(h.rows(), h.cols());
			Eigen::MatrixXd Gminus = Eigen::MatrixXd::Zero(h.rows(), h.cols());

			for (int i = 0; i < numberOfOrbitals; ++i)
				for (int j = 0; j < numberOfOrbitals; ++j)
					for (int k = 0; k < numberOfOrbitals; ++k)
						for (int l = 0; l < numberOfOrbitals; ++l)
						{
							double coulomb = integralsRepository.getElectronElectron(i, j, k, l);
							double exchange = integralsRepository.getElectronElectron(i, l, k, j);

							Gplus(i, j) += DensityMatrixPlus(k, l) * (coulomb - exchange) + DensityMatrixMinus(k, l) * coulomb; // the beta electrons interact with the alpha ones with coulomb interaction, too
							Gminus(i, j) += DensityMatrixMinus(k, l) * (coulomb - exchange) + DensityMatrixPlus(k, l) * coulomb; // the alpha electrons interact with the beta ones with coulomb interaction, too
						}

			FockMatrixPlus = h + Gplus;
			FockMatrixMinus = h + Gminus;
		}
	}


	void UnrestrictedHartreeFock::CalculateEnergy(const Eigen::VectorXd& eigenvalsplus, const Eigen::VectorXd& eigenvalsminus, const Eigen::MatrixXd& calcDensityMatrixPlus, const Eigen::MatrixXd& calcDensityMatrixMinus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/)
	{
		totalEnergy = 0;

		// simply adding energy levels does not work, because of over counting interaction		

		// ************************************************************************************
		// this is simpler but gives the same results as the one below!!!!!!!!

		for (int i = 1; i < h.rows(); ++i)
			for (int j = 0; j < i; ++j)
				totalEnergy += (calcDensityMatrixPlus(i, j) + calcDensityMatrixMinus(i, j)) * h(j, i);

		// only the values below the diagonal were added
		// *2 to have the sum of all elements except those on diagonal
		totalEnergy *= 2.;

		// now add the diagonal elements, too
		for (int i = 0; i < h.rows(); ++i) totalEnergy += (calcDensityMatrixPlus(i, i) + calcDensityMatrixMinus(i, i)) * h(i, i);


		for (unsigned int level = 0; level < occupiedPlus.size(); ++level)
			if (occupiedPlus[level]) totalEnergy += eigenvalsplus(level);
		for (unsigned int level = 0; level < occupiedMinus.size(); ++level)
			if (occupiedMinus[level]) totalEnergy += eigenvalsminus(level);


		HOMOEnergy = max(nrOccupiedLevelsPlus ? eigenvalsplus(nrOccupiedLevelsPlus - 1) : 0, nrOccupiedLevelsMinus ? eigenvalsminus(nrOccupiedLevelsMinus - 1) : 0);

		// ***************************************************************************************

		/*
		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				totalEnergy += (calcPplus(i, j) + calcPminus(i, j)) * h(j, i);


		for (unsigned int i = 0; i < h.rows(); ++i)
			for (unsigned int j = 0; j < h.cols(); ++j)
				totalEnergy += calcPplus(i, j) * Fplus(i, j) + calcPminus(i, j) * Fminus(i, j);
				*/
				// *******************************************************************************************


		totalEnergy /= 2.;

		totalEnergy += nuclearRepulsionEnergy;
	}

	double UnrestrictedHartreeFock::CalculateMp2Energy()
	{
		mp2Energy = 0;

		// TODO: calculate it

		return mp2Energy;
	}

}