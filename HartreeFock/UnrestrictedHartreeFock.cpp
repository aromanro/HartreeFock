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

		
		Pplus = Eigen::MatrixXd::Zero(h.rows(), h.cols());
		Pminus = Eigen::MatrixXd::Zero(h.rows(), h.cols());
		
		const unsigned int electronsNumber = molecule->ElectronsNumber();
		nrOccupiedLevelsMinus = static_cast<unsigned int>(floor(electronsNumber / 2.));
		nrOccupiedLevelsPlus = (electronsNumber % 2 ? nrOccupiedLevelsMinus + 1 : nrOccupiedLevelsMinus);

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
		Eigen::MatrixXd Fplus; 
		Eigen::MatrixXd Fminus;

		InitFockMatrices(iter, Fplus, Fminus);

		// ***************************************************************************************************************************

		// solve the Pople-Nesbet–Berthier equations

		Eigen::MatrixXd Fplusprime = Vt * Fplus * V;
		Eigen::MatrixXd Fminusprime = Vt * Fminus * V;

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esplus(Fplusprime);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esminus(Fminusprime);

		const Eigen::MatrixXd& Cplusprime = esplus.eigenvectors();
		const Eigen::MatrixXd& Cminusprime = esminus.eigenvectors();


		Eigen::MatrixXd Cplus = V * Cplusprime;
		Eigen::MatrixXd Cminus = V * Cminusprime;

		// normalize them
		//NormalizeC(Cplus, nrOccupiedLevelsPlus);
		//NormalizeC(Cminus, nrOccupiedLevelsMinus);

		//***************************************************************************************************************

		// calculate the density matrices

		Eigen::MatrixXd newPplus = Eigen::MatrixXd::Zero(h.rows(), h.cols());
		Eigen::MatrixXd newPminus = Eigen::MatrixXd::Zero(h.rows(), h.cols());

		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
			{
				for (unsigned int vec = 0; vec < nrOccupiedLevelsPlus; ++vec) // only eigenstates that are occupied 
					newPplus(i, j) += Cplus(i, vec) * Cplus(j, vec);

				for (unsigned int vec = 0; vec < nrOccupiedLevelsMinus; ++vec) // only eigenstates that are occupied 
					newPminus(i, j) += Cminus(i, vec) * Cminus(j, vec);
			}

		//**************************************************************************************************************

		const Eigen::VectorXd& eigenvalsplus = esplus.eigenvalues();
		const Eigen::VectorXd& eigenvalsminus = esminus.eigenvalues();

		CalculateEnergy(eigenvalsplus, eigenvalsminus, newPplus, newPminus/*, Fplus, Fminus*/);

		TRACE("Step: %d Energy: %f\n", iter, totalEnergy);

		// ***************************************************************************************************
		// go to the next density matrices
		// use mixing if alpha is set less then one

		Pplus = alpha * newPplus + (1. - alpha) * Pplus;
		Pminus = alpha * newPminus + (1. - alpha) * Pminus;
	}

	void UnrestrictedHartreeFock::InitFockMatrices(int iter, Eigen::MatrixXd& Fplus, Eigen::MatrixXd& Fminus) const
	{
		// this could be made faster knowing that the matrix should be symmetric
		// but it would be less expressive so I'll let it as it is
		// maybe I'll improve it later
		// anyway, the slower part is dealing with electron-electron integrals

		if (0 == iter)
		{
			if (initGuess > 0)
			{
				Fplus.resize(h.rows(), h.cols());
				Fminus.resize(h.rows(), h.cols());

				for (int i = 0; i < h.rows(); ++i)
					for (int j = 0; j < h.rows(); ++j)
						Fplus(i, j) = Fminus(i, j) = initGuess * overlapMatrix.matrix(i, j) * (h(i, i) + h(j, j)) / 2.;
			}
			else
			{
				Fplus = h;
				Fminus = h;
			}			

			if (addAsymmetry)
			{
				Fplus(0, 1) += asymmetry;
				Fplus(1, 0) = Fplus(0, 1);
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

							Gplus(i, j) += Pplus(k, l) * (coulomb - exchange) + Pminus(k, l) * coulomb; // the beta electrons interact with the alpha ones with coulomb interaction, too
							Gminus(i, j) += Pminus(k, l) * (coulomb - exchange) + Pplus(k, l) * coulomb; // the alpha electrons interact with the beta ones with coulomb interaction, too
						}

			Fplus = h + Gplus;
			Fminus = h + Gminus;
		}
	}


	void UnrestrictedHartreeFock::CalculateEnergy(const Eigen::VectorXd& eigenvalsplus, const Eigen::VectorXd& eigenvalsminus, const Eigen::MatrixXd& calcPplus, const Eigen::MatrixXd& calcPminus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/)
	{
		totalEnergy = 0;

		// simply adding energy levels does not work, because of over counting interaction		

		// ************************************************************************************
		// this is simpler but gives the same results as the one below!!!!!!!!

		for (int i = 1; i < h.rows(); ++i)
			for (int j = 0; j < i; ++j)
				totalEnergy += (calcPplus(i, j) + calcPminus(i, j)) * h(j, i);

		// only the values below the diagonal were added
		// *2 to have the sum of all elements except those on diagonal
		totalEnergy *= 2.; 
		
		// now add the diagonal elements, too
		for (int i = 0; i < h.rows(); ++i) totalEnergy += (calcPplus(i, i) + calcPminus(i, i)) * h(i, i);


		for (unsigned int level = 0; level < nrOccupiedLevelsPlus; ++level)	totalEnergy += eigenvalsplus(level);
		for (unsigned int level = 0; level < nrOccupiedLevelsMinus; ++level) totalEnergy += eigenvalsminus(level);


		HOMOEnergy = max(eigenvalsplus(nrOccupiedLevelsPlus - 1), eigenvalsminus(nrOccupiedLevelsMinus - 1));

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