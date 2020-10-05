#pragma once


#include <Eigen\eigen>

#include "QuantumNumbers.h"
#include "GaussianIntegral.h"

namespace GaussianIntegrals {

	// The overlap integrals are implemented using guidance from Mathematica Journal
	// Evaluation of Gaussian Molecular Integrals, I. Overlap Integrals
	// Here is a link: http://www.mathematica-journal.com/2012/02/evaluation-of-gaussian-molecular-integrals/

	// overlap integrals are a special case of moment integrals 
	// the recurrence relations are the same but there is one more for transfer of angular momentum from the third center
	// see the GaussianMoment class comments for details

	class GaussianOverlap : public GaussianIntegral {
	public:
		Eigen::MatrixXd matrixX;
		Eigen::MatrixXd matrixY;
		Eigen::MatrixXd matrixZ;

		double factor;

		GaussianOverlap();
		GaussianOverlap(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2);

		void Reset(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2);

		double operator()(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const { return getOverlap(QN1, QN2); }
		double getOverlap(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const;
	protected:
		void CalculateOverlap(Eigen::MatrixXd& matrix, double alpha1, double alpha2, double center1, double center2, unsigned int maxQN1, unsigned int maxQN2);
	};
}