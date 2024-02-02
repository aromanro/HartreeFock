#include "stdafx.h"
#include "QuantumMatrix.h"

namespace Matrices {

	QuantumMatrix::QuantumMatrix(GaussianIntegrals::IntegralsRepository* repository)
		: integralsRepository(repository), nrBasis(0)
	{
		Setup();
	}

	void QuantumMatrix::SetRepository(GaussianIntegrals::IntegralsRepository* repository)
	{
		integralsRepository = repository;

		Setup();
	}

	void QuantumMatrix::Setup()
	{
		if (integralsRepository)
		{
			nrBasis = integralsRepository->getMolecule()->CountNumberOfContractedGaussians();
			matrix = Eigen::MatrixXd::Zero(nrBasis, nrBasis);
		}
	}

	void QuantumMatrix::clear()
	{
		matrix.resize(0, 0);
	}


	void OverlapMatrix::Calculate()
	{
		const Systems::Molecule* molecule = integralsRepository->getMolecule();

		int i = 0;
		for (const auto& atom1 : molecule->atoms)
			for (const auto& shell1 : atom1.shells)
				for (const auto& orbital1 : shell1.basisFunctions)
				{
					int j = 0;

					for (const auto& atom2 : molecule->atoms)
						for (const auto& shell2 : atom2.shells)
							for (const auto& orbital2 : shell2.basisFunctions)
							{
								if (j >= i)
								{
									matrix(i, j) = integralsRepository->getOverlap(orbital1, orbital2);
									if (i != j) matrix(j, i) = matrix(i, j);
								}
								
								++j;
							}
					++i;
				}
	}

	void MomentMatrix::Calculate()
	{
		const Systems::Molecule* molecule = integralsRepository->getMolecule();
		
		matrixY = Eigen::MatrixXd::Zero(nrBasis, nrBasis);
		matrixZ = Eigen::MatrixXd::Zero(nrBasis, nrBasis);

		int i = 0;
		for (const auto& atom1 : molecule->atoms)
			for (const auto& shell1 : atom1.shells)
				for (const auto& orbital1 : shell1.basisFunctions)
				{
					int j = 0;

					for (const auto& atom2 : molecule->atoms)
						for (const auto& shell2 : atom2.shells)
							for (const auto& orbital2 : shell2.basisFunctions)
							{
								if (j >= i)
								{
									matrix(i, j) = integralsRepository->getMomentX(orbital1, orbital2);
									matrixY(i, j) = integralsRepository->getMomentY(orbital1, orbital2);
									matrixZ(i, j) = integralsRepository->getMomentZ(orbital1, orbital2);
									if (i != j) 
									{
										matrix(j, i) = matrix(i, j);
										matrixY(j, i) = matrixY(i, j);
										matrixZ(j, i) = matrixZ(i, j);
									}
								}

								++j;
							}
					++i;
				}
	}


	void KineticMatrix::Calculate()
	{
		const Systems::Molecule* molecule = integralsRepository->getMolecule();

		int i = 0;
		for (const auto& atom1 : molecule->atoms)
			for (const auto& shell1 : atom1.shells)
				for (const auto& orbital1 : shell1.basisFunctions)
				{
					int j = 0;

					for (const auto& atom2 : molecule->atoms)
						for (const auto& shell2 : atom2.shells)
							for (const auto& orbital2 : shell2.basisFunctions)
							{
								if (j >= i)
								{
									matrix(i, j) = integralsRepository->getKinetic(orbital1, orbital2);
									if (i != j) matrix(j, i) = matrix(i, j);
								}

								++j;
							}
					++i;
				}
	}



	void NuclearMatrix::Calculate()
	{
		const Systems::Molecule* molecule = integralsRepository->getMolecule();

		int i = 0;
		for (const auto& atom1 : molecule->atoms)
			for (const auto& shell1 : atom1.shells)
				for (const auto& orbital1 : shell1.basisFunctions)
				{
					int j = 0;

					for (const auto& atom2 : molecule->atoms)
						for (const auto& shell2 : atom2.shells)
							for (const auto& orbital2 : shell2.basisFunctions)
							{
								if (j >= i)
								{
									for (const auto& atom : molecule->atoms)
									        matrix(i, j) -= atom.Z * integralsRepository->getNuclear(atom, &orbital1, &orbital2);

									if (i != j) matrix(j, i) = matrix(i, j);
								}

								++j;
							}
					++i;
				}
	}
}


