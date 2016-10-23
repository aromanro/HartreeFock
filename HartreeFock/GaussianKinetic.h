#pragma once

#include <Eigen\eigen>

#include "QuantumNumbers.h"
#include "GaussianIntegral.h"

namespace Orbitals {
	class GaussianOrbital;
}

namespace GaussianIntegrals {

	// The kinetic integrals are implemented using guidance from Mathematica Journal
	// Evaluation of Gaussian Molecular Integrals, II Kinetic-Energy Integrals
	// Here is a link: http://www.mathematica-journal.com/2013/01/evaluation-of-gaussian-molecular-integrals-2/


	class GaussianKinetic : public GaussianIntegral
	{
	protected:
		const Orbitals::GaussianOrbital *m_gaussian1;
		const Orbitals::GaussianOrbital *m_gaussian2;
		const GaussianOverlap *m_overlap;

		Eigen::MatrixXd matrixX;
		Eigen::MatrixXd matrixY;
		Eigen::MatrixXd matrixZ;

	public:
		GaussianKinetic(const Orbitals::GaussianOrbital* gaussian1, const Orbitals::GaussianOrbital* gaussian2, const GaussianOverlap* overlap);
		~GaussianKinetic();

		void Reset(double alpha1, double alpha2, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2);

		double operator()(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const { return getKinetic(QN1, QN2); }
		double getKinetic(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const;
	protected:
		void CalculateKinetic(Eigen::MatrixXd& matrix, const Eigen::MatrixXd& overlap_matrix, double alpha1, double alpha2, unsigned int maxQN1, unsigned int maxQN2);
	};

}
