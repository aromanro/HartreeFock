#include "stdafx.h"
#include "GaussianMoment.h"

#include "MathUtils.h"


namespace GaussianIntegrals {


	GaussianMoment::GaussianMoment()
		: factor(0)
	{
	}

	GaussianMoment::GaussianMoment(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Vector3D<double>& center3, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2)
		: factor(0)
	{
		Reset(alpha1, alpha2, center1, center2, center3, maxQN1, maxQN2);
	}


	void GaussianMoment::Reset(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Vector3D<double>& center3, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2)
	{
		matrixX = Eigen::MatrixXd::Zero(3ULL + maxQN1.l + maxQN2.l, maxQN2.l + 2ULL);
		matrixY = Eigen::MatrixXd::Zero(3ULL + maxQN1.m + maxQN2.m, maxQN2.m + 2ULL);
		matrixZ = Eigen::MatrixXd::Zero(3ULL + maxQN1.n + maxQN2.n, maxQN2.n + 2ULL);

		matrixX1 = Eigen::MatrixXd::Zero(2ULL + maxQN1.l + maxQN2.l, maxQN2.l + 2ULL);
		matrixY1 = Eigen::MatrixXd::Zero(2ULL + maxQN1.m + maxQN2.m, maxQN2.m + 2ULL);
		matrixZ1 = Eigen::MatrixXd::Zero(2ULL + maxQN1.n + maxQN2.n, maxQN2.n + 2ULL);

		CalculateMoment(matrixX, matrixX1, alpha1, alpha2, center1.X, center2.X, center3.X, maxQN1.l, maxQN2.l);
		CalculateMoment(matrixY, matrixY1, alpha1, alpha2, center1.Y, center2.Y, center3.Y, maxQN1.m, maxQN2.m);
		CalculateMoment(matrixZ, matrixZ1, alpha1, alpha2, center1.Z, center2.Z, center3.Z, maxQN1.n, maxQN2.n);

		const Vector3D<double> dif = center1 - center2;
		const double oneDivAlpha1PlusAlpha2 = 1. / (alpha1 + alpha2);
		factor = exp(-alpha1 * alpha2 * oneDivAlpha1PlusAlpha2 * dif * dif) * pow(M_PI * oneDivAlpha1PlusAlpha2, 3. / 2.);
	}

	double GaussianMoment::getMoment(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2, bool momentX, bool momentY, bool momentZ) const
	{
		return factor * (momentX ? matrixX1(QN1.l, QN2.l) : matrixX(QN1.l, QN2.l)) * (momentY ? matrixY1(QN1.m, QN2.m) : matrixY(QN1.m, QN2.m)) * (momentZ ? matrixZ1(QN1.n, QN2.n) : matrixZ(QN1.n, QN2.n));
	}

	double GaussianMoment::getOverlap(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const
	{
		return factor * matrixX(QN1.l, QN2.l) * matrixY(QN1.m, QN2.m) * matrixZ(QN1.n, QN2.n); // or getMoment(QN1, QN2, false, false, false) but this is slightly faster
	}

	void GaussianMoment::CalculateMoment(Eigen::MatrixXd& matrix, Eigen::MatrixXd& matrix1, double alpha1, double alpha2, double center1, double center2, double center3, unsigned int maxQN1, unsigned int maxQN2)
	{
		const double alpha = alpha1 + alpha2;
		const double productCenter = (alpha1 * center1 + alpha2 * center2) / alpha;
		const double dif = center1 - center2;
		const double dif1 = center1 - center3;
		const double difCenter = productCenter - center1;
		const double oneDiv2alpha = 1. / (2. * alpha);

		matrix(0, 0) = 1;
		matrix(1, 0) = difCenter;

		// recurrence index
		unsigned int limit = maxQN1 + maxQN2 + 2;

		// vertical recurrence relation - the same as for overlap
		for (unsigned int i = 2; i <= limit; ++i)
			matrix(i, 0) = difCenter * matrix(i - 1, 0) + (i - 1.) * oneDiv2alpha * matrix(i - 2, 0);


		--limit;

		// for first column of matrix1 as well
		for (unsigned int i = 0; i < limit; ++i)
			matrix1(i, 0) = matrix(i + 1ULL, 0) + dif1 * matrix(i, 0);


		// transfer equation - the horizontal recurrence relation
		const unsigned int horzLimit = maxQN2 + 1;
		for (unsigned int j = 1; j <= horzLimit; ++j, --limit)
		{
			const unsigned int jminus1 = j - 1ULL;
			
			// *** this is is exactly as for overlap ********************************

			for (unsigned int i = 0; i <= limit; ++i)
				matrix(i, j) = matrix(i + 1ULL, jminus1) + dif * matrix(i, jminus1);

			// *******************************************************************************

			for (unsigned int i = 0; i < limit; ++i)
				matrix1(i, j) = matrix(i + 1ULL, j) + dif1 * matrix(i, j);
		}

		matrix = matrix.block(0, 0, maxQN1 + 1ULL, maxQN2 + 1ULL).eval();
		matrix1 = matrix1.block(0, 0, maxQN1 + 1ULL, maxQN2 + 1ULL).eval();
	}

}
