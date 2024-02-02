#pragma once

#include <Eigen\eigen>

#include "IntegralsRepository.h"

namespace Matrices {

	class QuantumMatrix
	{
	public:
		Eigen::MatrixXd matrix;

		QuantumMatrix(GaussianIntegrals::IntegralsRepository* repository = nullptr);
		virtual ~QuantumMatrix() = default;

		virtual void Calculate() = 0;

		void SetRepository(GaussianIntegrals::IntegralsRepository* repository);

	protected:
		void Setup();

		GaussianIntegrals::IntegralsRepository* integralsRepository;
		unsigned int nrBasis;

	public:
		void clear();
	};


	class OverlapMatrix : public QuantumMatrix
	{
	public:
		OverlapMatrix(GaussianIntegrals::IntegralsRepository* repository = nullptr) : QuantumMatrix(repository) { if (integralsRepository) Calculate(); }

		virtual void Calculate();
	};


	class MomentMatrix : public QuantumMatrix
	{
	public:
		Eigen::MatrixXd matrixY;
		Eigen::MatrixXd matrixZ;

		MomentMatrix(GaussianIntegrals::IntegralsRepository* repository = nullptr) : QuantumMatrix(repository) { if (integralsRepository) Calculate(); }

		virtual void Calculate();
	};


	class KineticMatrix : public QuantumMatrix 
	{
	public:
		KineticMatrix(GaussianIntegrals::IntegralsRepository* repository = nullptr) : QuantumMatrix(repository) { if (integralsRepository) Calculate(); }

		virtual void Calculate();
	};

	class NuclearMatrix : public QuantumMatrix
	{
	public:
		NuclearMatrix(GaussianIntegrals::IntegralsRepository* repository = nullptr) : QuantumMatrix(repository) { if (integralsRepository) Calculate(); }

		virtual void Calculate();
	};

}

