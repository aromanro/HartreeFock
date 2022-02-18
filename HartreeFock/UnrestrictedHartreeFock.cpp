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


	bool UnrestrictedHartreeFock::DIISStep(int iter, Eigen::MatrixXd& FockMatrixPlus, Eigen::MatrixXd& FockMatrixMinus)
	{
		bool UsedDIIS = false;

		if (UseDIIS && iter && iter < maxDIISiterations)
		{
			const Eigen::MatrixXd errorMatrixPlus = overlapMatrix.matrix * DensityMatrixPlus * FockMatrixPlus - FockMatrixPlus * DensityMatrixPlus * overlapMatrix.matrix;
			diisPlus.AddValueAndError(FockMatrixPlus, errorMatrixPlus);

			const Eigen::MatrixXd errorMatrixMinus = overlapMatrix.matrix * DensityMatrixMinus * FockMatrixMinus - FockMatrixMinus * DensityMatrixMinus * overlapMatrix.matrix;

			if (diisMinus.AddValueAndError(FockMatrixMinus, errorMatrixMinus))
			{
				// use DIIS

				lastErrorEst = diisMinus.Estimate(FockMatrixMinus) + diisPlus.Estimate(FockMatrixPlus);
				UsedDIIS = true;
			}
		}
		else lastErrorEst = 0;

		return UsedDIIS;
	}


	double UnrestrictedHartreeFock::Step(int iter)
	{
		// *****************************************************************************************************************

		// the Fock matrices
		Eigen::MatrixXd FockMatrixPlus;
		Eigen::MatrixXd FockMatrixMinus;

		InitFockMatrices(iter, FockMatrixPlus, FockMatrixMinus);

		// will be used for DIIS

		const bool UsedDIIS = DIISStep(iter, FockMatrixPlus, FockMatrixMinus);
		
		// ***************************************************************************************************************************

		// solve the Pople-Nesbet–Berthier equations

		 // orthogonalize
		const Eigen::MatrixXd FockMatrixPlusTransformed = Vt * FockMatrixPlus * V;
		const Eigen::MatrixXd FockMatrixMinusTransformed = Vt * FockMatrixMinus * V;

		if (FockMatrixPlusTransformed.rows() > 1)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esplus(FockMatrixPlusTransformed);
			eigenvalsplus = esplus.eigenvalues();
			const Eigen::MatrixXd& Cplusprime = esplus.eigenvectors();

			Cplus = V * Cplusprime; // transform back the eigenvectors into the original non-orthogonalized AO basis
		}
		else
		{
			eigenvalsplus = FockMatrixPlusTransformed(0, 0) * Eigen::VectorXd::Ones(1);
			Cplus = V * Eigen::MatrixXd::Ones(1, 1);
		}

		if (FockMatrixMinusTransformed.rows() > 1)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esminus(FockMatrixMinusTransformed);
			eigenvalsminus = esminus.eigenvalues();
			const Eigen::MatrixXd& Cminusprime = esminus.eigenvectors();

			Cminus = V * Cminusprime; // transform back the eigenvectors into the original non-orthogonalized AO basis
		}
		else
		{
			eigenvalsminus = FockMatrixMinusTransformed(0, 0) * Eigen::VectorXd::Ones(1);
			Cminus = V * Eigen::MatrixXd::Ones(1, 1);
		}

		// normalize them - in some rare cases it seems to help
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
		double rmsD = rmsDensityMatricesDif.sum();

		rmsDensityMatricesDif = newDensityMatrixMinus - DensityMatrixMinus;
		rmsDensityMatricesDif = rmsDensityMatricesDif.cwiseProduct(rmsDensityMatricesDif);
		rmsD += rmsDensityMatricesDif.sum();

		rmsD = sqrt(rmsD);

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
						FockMatrixPlus(i, j) = FockMatrixMinus(i, j) = initGuess * overlapMatrix.matrix(i, j) * (h(i, i) + h(j, j)) * 0.5;
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


	void UnrestrictedHartreeFock::CalculateEnergy(const Eigen::VectorXd& eigenvalsPlus, const Eigen::VectorXd& eigenvalsMinus, const Eigen::MatrixXd& calcDensityMatrixPlus, const Eigen::MatrixXd& calcDensityMatrixMinus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/)
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
			if (occupiedPlus[level]) totalEnergy += eigenvalsPlus(level);
		for (unsigned int level = 0; level < occupiedMinus.size(); ++level)
			if (occupiedMinus[level]) totalEnergy += eigenvalsMinus(level);


		HOMOEnergy = max(nrOccupiedLevelsPlus ? eigenvalsPlus(nrOccupiedLevelsPlus - 1) : 0, nrOccupiedLevelsMinus ? eigenvalsMinus(nrOccupiedLevelsMinus - 1) : 0);

		// ***************************************************************************************

		/*
		for (int i = 0; i < h.rows(); ++i)
			for (int j = 0; j < h.cols(); ++j)
				totalEnergy += (calcDensityMatrixPlus(i, j) + calcDensityMatrixMinus(i, j)) * h(j, i);


		for (unsigned int i = 0; i < h.rows(); ++i)
			for (unsigned int j = 0; j < h.cols(); ++j)
				totalEnergy += calcDensityMatrixPlus(i, j) * Fplus(i, j) + calcDensityMatrixMinus(i, j) * Fminus(i, j);
				*/
				// *******************************************************************************************


		totalEnergy /= 2.;

		totalEnergy += nuclearRepulsionEnergy;
	}

	double UnrestrictedHartreeFock::CalculateMp2Energy()
	{
		GaussianIntegrals::MolecularOrbitalsIntegralsRepository MP2repo(integralsRepository);

		mp2Energy = CalculateMp2EnergyPlus(MP2repo);
		mp2Energy += CalculateMp2EnergyMinus(MP2repo);

		mp2Energy *= 0.25;

		return mp2Energy;
	}

	double UnrestrictedHartreeFock::CalculateMp2EnergyPlus(GaussianIntegrals::MolecularOrbitalsIntegralsRepository& MP2repo) const
	{
		double mp2Energyt = 0;

		for (int i = 0; i < numberOfOrbitals; ++i)
		{
			if (i >= occupiedPlus.size() || !occupiedPlus[i]) continue; // only occupied

			for (int j = 0; j < numberOfOrbitals; ++j)
			{
				if (j >= occupiedPlus.size() || !occupiedPlus[j]) continue; // only occupied

				for (int a = 0; a < numberOfOrbitals; ++a)
				{
					if (a < occupiedPlus.size() && occupiedPlus[a]) continue; // only unoccupied

					for (int b = 0; b < numberOfOrbitals; ++b)
					{
						if (b < occupiedPlus.size() && occupiedPlus[b]) continue;  // only unoccupied

						const double Esumdif = eigenvalsplus(i) + eigenvalsplus(j) - eigenvalsplus(a) - eigenvalsplus(b);

						const double eeiajb = MP2repo.getElectronElectron(i, a, j, b, Cplus);

						const double partE = eeiajb * eeiajb / Esumdif;

						mp2Energyt += partE;
					}
				}
			}
		}

		return mp2Energyt;
	}

	double UnrestrictedHartreeFock::CalculateMp2EnergyMinus(GaussianIntegrals::MolecularOrbitalsIntegralsRepository& MP2repo) const
	{
		double mp2Energyt = 0;

		for (int i = 0; i < numberOfOrbitals; ++i)
		{
			if (i >= occupiedMinus.size() || !occupiedMinus[i]) continue; // only occupied

			for (int j = 0; j < numberOfOrbitals; ++j)
			{
				if (j >= occupiedMinus.size() || !occupiedMinus[j]) continue; // only occupied

				for (int a = 0; a < numberOfOrbitals; ++a)
				{
					if (a < occupiedMinus.size() && occupiedMinus[a]) continue; // only unoccupied

					for (int b = 0; b < numberOfOrbitals; ++b)
					{
						if (b < occupiedMinus.size() && occupiedMinus[b]) continue;  // only unoccupied

						const double Esumdif = eigenvalsminus(i) + eigenvalsminus(j) - eigenvalsminus(a) - eigenvalsminus(b);

						const double eeiajb = MP2repo.getElectronElectron(i, a, j, b, Cminus);

						const double partE = eeiajb * eeiajb / Esumdif;

						mp2Energyt += partE;
					}
				}
			}
		}

		return mp2Energyt;
	}


	double UnrestrictedHartreeFock::CalculateAtomicCharge(int atom) const
	{
		if (!integralsRepository.m_Molecule || integralsRepository.m_Molecule->atoms.size() <= atom) return 0;

		double result = integralsRepository.m_Molecule->atoms[atom].Z;

		int orbLowLimit = 0;
		int orbHighLimit = 0;
		for (int i = 0; i < integralsRepository.m_Molecule->atoms.size(); ++i)
		{
			const int numBasisFunctions = integralsRepository.m_Molecule->atoms[i].CountNumberOfContractedGaussians();
			if (i == atom)
			{
				orbHighLimit = orbLowLimit + numBasisFunctions;
				break;
			}
			orbLowLimit += numBasisFunctions;
		}

		// this is not so efficient, being done for each atom if needed instead of once for all atoms, but it's ok
		const Eigen::MatrixXd DSPlus = DensityMatrixPlus * overlapMatrix.matrix;
		const Eigen::MatrixXd DSMinus = DensityMatrixMinus * overlapMatrix.matrix;

		for (int i = orbLowLimit; i < orbHighLimit; ++i)
			result -= DSPlus(i, i) + DSMinus(i, i);

		return result;
	}


	Vector3D<double> UnrestrictedHartreeFock::GetMoment() const
	{
		Vector3D<double> moment = GetNuclearMoment();

		const unsigned long long int sz = DensityMatrixPlus.cols();

		for (unsigned int i = 0; i < sz; ++i)
			for (unsigned int j = 0; j < sz; ++j)
			{
				const double dval = DensityMatrixPlus(i, j) + DensityMatrixMinus(i, j);
				moment.X -= dval * momentMatrix.matrix(i, j);
				moment.Y -= dval * momentMatrix.matrixY(i, j);
				moment.Z -= dval * momentMatrix.matrixZ(i, j);
			}

		return moment;
	}

}