#include "stdafx.h"
#include "UnrestrictedHartreeFock.h"


namespace HartreeFock {


	UnrestrictedHartreeFock::UnrestrictedHartreeFock(int iterations)
		: HartreeFockAlgorithm(iterations), totalEnergy(0), nrOccupiedLevelsPlus(0), nrOccupiedLevelsMinus(0), addAsymmetry(true), asymmetry(0.1)
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
	}

	void UnrestrictedHartreeFock::Step(int iter)
	{		
		// *****************************************************************************************************************

		// the Fock matrices
		Eigen::MatrixXd FockMatrixPlus; 
		Eigen::MatrixXd FockMatrixMinus;

		InitFockMatrices(iter, FockMatrixPlus, FockMatrixMinus);

		// ***************************************************************************************************************************

		// solve the Pople-Nesbet–Berthier equations

		Eigen::MatrixXd FockMatrixPlusTransformed = Vt * FockMatrixPlus * V;
		Eigen::MatrixXd FockMatrixMinusTransformed = Vt * FockMatrixMinus * V;

		Eigen::MatrixXd Cplus;
		Eigen::VectorXd eigenvalsplus;

		if (FockMatrixPlusTransformed.rows() > 1)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esplus(FockMatrixPlusTransformed);
			const Eigen::MatrixXd& Cplusprime = esplus.eigenvectors();
			Cplus = V * Cplusprime;
			eigenvalsplus = esplus.eigenvalues();
		}
		else
		{
			Cplus = V * Eigen::MatrixXd::Ones(1, 1);
			eigenvalsplus = FockMatrixPlusTransformed(0, 0) * Eigen::VectorXd::Ones(1);
		}

		
		Eigen::MatrixXd Cminus;
		Eigen::VectorXd eigenvalsminus;

		if (FockMatrixMinusTransformed.rows() > 1)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esminus(FockMatrixMinusTransformed);
			const Eigen::MatrixXd& Cminusprime = esminus.eigenvectors();
			Cminus = V * Cminusprime;
			eigenvalsminus = esminus.eigenvalues();
		}
		else
		{
			Cminus = V * Eigen::MatrixXd::Ones(1, 1);
			eigenvalsminus = FockMatrixMinusTransformed(0, 0) * Eigen::VectorXd::Ones(1);
		}

		// normalize them
		//NormalizeC(Cplus, nrOccupiedLevelsPlus);
		//NormalizeC(Cminus, nrOccupiedLevelsMinus);

		//***************************************************************************************************************

		// calculate the density matrices

		Eigen::MatrixXd newDensityMatrixPlus = Eigen::MatrixXd::Zero(h.rows(), h.cols());
		Eigen::MatrixXd newDensityMatrixMinus = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
			{
				for (unsigned int vec = 0; vec < nrOccupiedLevelsPlus; ++vec) // only eigenstates that are occupied 
					newDensityMatrixPlus(i, j) += Cplus(i, vec) * Cplus(j, vec);

				for (unsigned int vec = 0; vec < nrOccupiedLevelsMinus; ++vec) // only eigenstates that are occupied 
					newDensityMatrixMinus(i, j) += Cminus(i, vec) * Cminus(j, vec);
			}

		//**************************************************************************************************************

		CalculateEnergy(eigenvalsplus, eigenvalsminus, newDensityMatrixPlus, newDensityMatrixMinus/*, Fplus, Fminus*/);

		TRACE("Step: %d Energy: %f\n", iter, totalEnergy);

		// ***************************************************************************************************
		// go to the next density matrices
		// use mixing if alpha is set less then one

		DensityMatrixPlus = alpha * newDensityMatrixPlus + (1. - alpha) * DensityMatrixPlus;
		DensityMatrixMinus = alpha * newDensityMatrixMinus + (1. - alpha) * DensityMatrixMinus;
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


		for (unsigned int level = 0; level < nrOccupiedLevelsPlus; ++level)	totalEnergy += eigenvalsplus(level);
		for (unsigned int level = 0; level < nrOccupiedLevelsMinus; ++level) totalEnergy += eigenvalsminus(level);


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

	double UnrestrictedHartreeFock::GetTotalEnergy() const
	{
		return totalEnergy;
	}

}