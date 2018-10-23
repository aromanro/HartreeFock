#include "stdafx.h"

#include "GaussianOverlap.h"
#include "GaussianKinetic.h"

#include "GaussianOrbital.h"

namespace GaussianIntegrals {

	GaussianKinetic::GaussianKinetic(const Orbitals::GaussianOrbital* gaussian1, const Orbitals::GaussianOrbital* gaussian2, const GaussianOverlap* overlap)
		: m_gaussian1(gaussian1), m_gaussian2(gaussian2), m_overlap(overlap)
	{		
	}


	GaussianKinetic::~GaussianKinetic()
	{
	}

	void GaussianKinetic::Reset(double alpha1, double alpha2, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2)
	{
		matrixX = Eigen::MatrixXd::Zero(maxQN1.l + 1, maxQN2.l + 1);
		matrixY = Eigen::MatrixXd::Zero(maxQN1.m + 1, maxQN2.m + 1);
		matrixZ = Eigen::MatrixXd::Zero(maxQN1.n + 1, maxQN2.n + 1);

		CalculateKinetic(matrixX, m_overlap->matrixX, alpha1, alpha2, maxQN1.l, maxQN2.l);
		CalculateKinetic(matrixY, m_overlap->matrixY, alpha1, alpha2, maxQN1.m, maxQN2.m);
		CalculateKinetic(matrixZ, m_overlap->matrixZ, alpha1, alpha2, maxQN1.n, maxQN2.n);
	}

	double GaussianKinetic::getKinetic(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const
	{		
		return m_overlap->factor * (
				matrixX(QN1.l, QN2.l) * m_overlap->matrixY(QN1.m, QN2.m) * m_overlap->matrixZ(QN1.n, QN2.n) +
				m_overlap->matrixX(QN1.l, QN2.l) * matrixY(QN1.m, QN2.m) * m_overlap->matrixZ(QN1.n, QN2.n) +
				m_overlap->matrixX(QN1.l, QN2.l) * m_overlap->matrixY(QN1.m, QN2.m) * matrixZ(QN1.n, QN2.n)
			);
	}

	void GaussianKinetic::CalculateKinetic(Eigen::MatrixXd& matrix, const Eigen::MatrixXd& overlap_matrix, double alpha1, double alpha2, unsigned int maxQN1, unsigned int maxQN2)
	{
		const double alphaProd = alpha1 * alpha2;
		matrix(0, 0) = 2 * alphaProd * overlap_matrix(1, 1);

		for (int i = 1; i <= static_cast<int>(maxQN1); ++i)
			matrix(i, 0) = -i * alpha2 * overlap_matrix(i - 1, 1) + 2. * alphaProd * overlap_matrix(i + 1, 1);

		for (int i = 1; i <= static_cast<int>(maxQN2); ++i)
			matrix(0, i) = -i * alpha1 * overlap_matrix(1, i - 1) + 2. * alphaProd * overlap_matrix(1, i + 1);

		for (int i = 1; i <= static_cast<int>(maxQN1); ++i)
			for (int j = 1; j <= static_cast<int>(maxQN2); ++j)
				matrix(i, j) = i * j * overlap_matrix(i - 1, j - 1) / 2.
					- j * alpha1 * overlap_matrix(i + 1, j - 1) 
					- i * alpha2 * overlap_matrix(i - 1, j + 1) 
					+ 2. * alphaProd * overlap_matrix(i + 1, j + 1);
	}

}