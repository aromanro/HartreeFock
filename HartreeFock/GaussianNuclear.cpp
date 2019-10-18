#include "stdafx.h"
#include "GaussianNuclear.h"

#include "GaussianOrbital.h"

#include "GaussianTwoElectrons.h"
#include "IntegralsRepository.h"

namespace GaussianIntegrals {

	GaussianNuclear::GaussianNuclear()
	{
	}


	GaussianNuclear::~GaussianNuclear()
	{
	}

	void GaussianNuclear::Reset(IntegralsRepository* repository, double alpha1, double alpha2, const Vector3D<double>& nucleus, const Vector3D<double>& center1, const Vector3D<double>& center2, unsigned int maxL1, unsigned int maxL2, bool calculateHorizontal)
	{
		const double alpha = alpha1 + alpha2;
		const Vector3D<double> Rp = (alpha1 * center1 + alpha2 * center2) / alpha;
		const Vector3D<double> difN = nucleus - Rp;		
		const Vector3D<double> dif = center1 - center2; 

		const unsigned int maxL = maxL1 + maxL2;
		const unsigned int size = maxL + 1;

		// auxiliary integrals
		const BoysFunctions& boys = repository->getBoysFunctions(maxL, alpha * (difN * difN));
		
		Orbitals::QuantumNumbers::QuantumNumbers maxQN(0, 0 , maxL);

		const double factor = 2. * M_PI / alpha * exp(-alpha1 * alpha2 / alpha * dif * dif);
		matrixCalc = Eigen::MatrixXd::Zero(maxQN.GetTotalCanonicalIndex() + 1ULL, size);

		for (unsigned int i = 0; i < size; ++i)	
			matrixCalc(0, i) = factor * boys.functions[i];

		VerticalRecursion(alpha, Rp, center1, difN, maxL);
		
		if (calculateHorizontal) HorizontalRecursion(dif, maxL1, maxL2);
		else matrixCalc = matrixCalc.block(0, 0, matrixCalc.rows(), 1).eval();
	}

	double GaussianNuclear::getNuclear(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const
	{
		return matrixCalc(QN1.GetTotalCanonicalIndex(), QN2.GetTotalCanonicalIndex());
	}


	void GaussianNuclear::VerticalRecursion(double alpha, const Vector3D<double> Rp, const Vector3D<double>& center1, const Vector3D<double>& difN, unsigned int maxL)
	{
		double difScalar, difNScalar;
		double N;

		const unsigned int size = maxL + 1;
		const auto difRp = Rp - center1;
		
		for (auto currentQN = Orbitals::QuantumNumbers::QuantumNumbers(1, 0, 0); currentQN < size; ++currentQN) // for each column
		{
			Orbitals::QuantumNumbers::QuantumNumbers prevQN = currentQN;
			Orbitals::QuantumNumbers::QuantumNumbers prevPrevQN = prevQN;

			const bool addPrevPrev = GaussianTwoElectrons::GetPrevAndPrevPrevAndScalarsForVerticalRecursion(currentQN, difRp, difN, prevQN, prevPrevQN, difScalar, difNScalar, N);

			unsigned int curIndex = currentQN.GetTotalCanonicalIndex();
			unsigned int prevIndex = prevQN.GetTotalCanonicalIndex();

			for (unsigned int m = 0; m < size - currentQN; ++m)
			{
				matrixCalc(curIndex, m) = difScalar * matrixCalc(prevIndex, m) + difNScalar * matrixCalc(prevIndex, m + 1ULL);

				if (addPrevPrev)
				{
					unsigned int prevPrevIndex = prevPrevQN.GetTotalCanonicalIndex();
					matrixCalc(curIndex, m) += N / (2. * alpha) * (matrixCalc(prevPrevIndex, m) - matrixCalc(prevPrevIndex, m + 1ULL));
				}
			}			
		}
	}

	void GaussianNuclear::HorizontalRecursion(const Vector3D<double>& dif, unsigned int maxL1, unsigned int maxL2)
	{
		unsigned int maxL = maxL1 + maxL2;

		Orbitals::QuantumNumbers::QuantumNumbers maxQN1(0, 0, maxL1);
		Orbitals::QuantumNumbers::QuantumNumbers maxQN2(0, 0, maxL2);

		const unsigned int limit = maxQN2.GetTotalCanonicalIndex() + 1;

		Eigen::MatrixXd matrixHoriz = Eigen::MatrixXd::Zero(matrixCalc.rows(), limit);
		matrixHoriz.col(0) = matrixCalc.col(0);
		
		for (auto currentQNj = Orbitals::QuantumNumbers::QuantumNumbers(1, 0, 0); currentQNj <= maxL2; GaussianTwoElectrons::IncrementQNandDecrementLimitIfNeeded(currentQNj, maxL))  //this for walks over the columns of the matrix
		{
			for (auto currentQN = Orbitals::QuantumNumbers::QuantumNumbers(0, 0, 0); currentQN < maxL; ++currentQN) // this for walks over the rows of the matrix
			{
				auto nextQN = currentQN;
				auto prevQNj = currentQNj;

				double difScalar = GaussianTwoElectrons::GetNextAndPrevQNAndScalarDiffForHorizontalRecursion(currentQN, currentQNj, dif, nextQN, prevQNj);

				unsigned int curIndex = currentQN.GetTotalCanonicalIndex();
				unsigned int curIndexJ = currentQNj.GetTotalCanonicalIndex();
				
				unsigned int nextIndex = nextQN.GetTotalCanonicalIndex();
				unsigned int prevIndexJ = prevQNj.GetTotalCanonicalIndex();

				matrixHoriz(curIndex, curIndexJ) = matrixHoriz(nextIndex, prevIndexJ) + difScalar * matrixHoriz(curIndex, prevIndexJ);
			}			
		}

		matrixCalc = matrixHoriz.block(0, 0, maxQN1.GetTotalCanonicalIndex() + 1ULL, limit);
	}

	
}

