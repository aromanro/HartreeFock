#include "stdafx.h"
#include "QuantumMatrix.h"

namespace Matrices {

	QuantumMatrix::QuantumMatrix(GaussianIntegrals::IntegralsRepository* repository)
		: integralsRepository(repository)
	{
		Setup();
	}


	QuantumMatrix::~QuantumMatrix()
	{
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
		Systems::Molecule* molecule = integralsRepository->getMolecule();

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

	void KineticMatrix::Calculate()
	{
		Systems::Molecule* molecule = integralsRepository->getMolecule();

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
		Systems::Molecule* molecule = integralsRepository->getMolecule();

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


