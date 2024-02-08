#include "stdafx.h"

#include "GaussianOrbital.h"

#include "GaussianOverlap.h"
#include "GaussianMoment.h"

#include "GaussianKinetic.h"


namespace GaussianIntegrals {

	GaussianKinetic::GaussianKinetic(const Orbitals::GaussianOrbital* gaussian1, const Orbitals::GaussianOrbital* gaussian2, const GaussianOverlap* overlap)
		: m_gaussian1(gaussian1), m_gaussian2(gaussian2), m_overlap(overlap), m_moment(nullptr)
	{		
	}

	GaussianKinetic::GaussianKinetic(const Orbitals::GaussianOrbital* gaussian1, const Orbitals::GaussianOrbital* gaussian2, const GaussianMoment* moment)
		: m_gaussian1(gaussian1), m_gaussian2(gaussian2), m_overlap(nullptr), m_moment(moment)
	{
	}

	void GaussianKinetic::Reset(double alpha1, double alpha2, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2)
	{
		matrixX = Eigen::MatrixXd::Zero(maxQN1.l + 1ULL, maxQN2.l + 1ULL);
		matrixY = Eigen::MatrixXd::Zero(maxQN1.m + 1ULL, maxQN2.m + 1ULL);
		matrixZ = Eigen::MatrixXd::Zero(maxQN1.n + 1ULL, maxQN2.n + 1ULL);

		if (m_moment)
		{
			CalculateKinetic(matrixX, m_moment->matrixX, alpha1, alpha2, maxQN1.l, maxQN2.l);
			CalculateKinetic(matrixY, m_moment->matrixY, alpha1, alpha2, maxQN1.m, maxQN2.m);
			CalculateKinetic(matrixZ, m_moment->matrixZ, alpha1, alpha2, maxQN1.n, maxQN2.n);
		}
		else
		{
			CalculateKinetic(matrixX, m_overlap->matrixX, alpha1, alpha2, maxQN1.l, maxQN2.l);
			CalculateKinetic(matrixY, m_overlap->matrixY, alpha1, alpha2, maxQN1.m, maxQN2.m);
			CalculateKinetic(matrixZ, m_overlap->matrixZ, alpha1, alpha2, maxQN1.n, maxQN2.n);
		}
	}

	double GaussianKinetic::getKinetic(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const
	{		
		if (m_moment)
			return m_moment->factor * (
				matrixX(QN1.l, QN2.l) * m_moment->matrixY(QN1.m, QN2.m) * m_moment->matrixZ(QN1.n, QN2.n) +
				m_moment->matrixX(QN1.l, QN2.l) * matrixY(QN1.m, QN2.m) * m_moment->matrixZ(QN1.n, QN2.n) +
				m_moment->matrixX(QN1.l, QN2.l) * m_moment->matrixY(QN1.m, QN2.m) * matrixZ(QN1.n, QN2.n)
				);


		return m_overlap->factor * (
				matrixX(QN1.l, QN2.l) * m_overlap->matrixY(QN1.m, QN2.m) * m_overlap->matrixZ(QN1.n, QN2.n) +
				m_overlap->matrixX(QN1.l, QN2.l) * matrixY(QN1.m, QN2.m) * m_overlap->matrixZ(QN1.n, QN2.n) +
				m_overlap->matrixX(QN1.l, QN2.l) * m_overlap->matrixY(QN1.m, QN2.m) * matrixZ(QN1.n, QN2.n)
			);
	}

	void GaussianKinetic::CalculateKinetic(Eigen::MatrixXd& matrix, const Eigen::MatrixXd& overlap_matrix, double alpha1, double alpha2, unsigned int maxQN1, unsigned int maxQN2)
	{
		const double twoAlphaProd = 2 * alpha1 * alpha2;
		matrix(0, 0) = twoAlphaProd * overlap_matrix(1, 1);

		for (int i = 1; i <= static_cast<int>(maxQN1); ++i)
			matrix(i, 0) = -i * alpha2 * overlap_matrix(i - 1ULL, 1) + twoAlphaProd * overlap_matrix(i + 1ULL, 1);

		for (int i = 1; i <= static_cast<int>(maxQN2); ++i)
			matrix(0, i) = -i * alpha1 * overlap_matrix(1, i - 1ULL) + twoAlphaProd * overlap_matrix(1, i + 1ULL);

		for (int i = 1; i <= static_cast<int>(maxQN1); ++i)
			for (int j = 1; j <= static_cast<int>(maxQN2); ++j)
			{
				const int iMinusOne = i - 1ULL;
				const int jMinusOne = j - 1ULL;
				const int iPlusOne = i + 1ULL;
				const int jPlusOne = j + 1ULL;

				matrix(i, j) = static_cast<double>(i) * j * overlap_matrix(iMinusOne, jMinusOne) * 0.5
					- j * alpha1 * overlap_matrix(iPlusOne, jMinusOne)
					- i * alpha2 * overlap_matrix(iMinusOne, jPlusOne)
					+ twoAlphaProd * overlap_matrix(iPlusOne, jPlusOne);
			}
	}

}