#pragma once

#include <Eigen\eigen>

#include "GaussianIntegral.h"
#include "QuantumNumbers.h"


namespace GaussianIntegrals {

	class IntegralsRepository;

	// for the nuclear integrals I looked first at Mathematica Journal,
	// Evaluation of Gaussian Molecular Integrals, III Nuclear-Electron Attraction Integrals
	// http://www.mathematica-journal.com/2014/12/evaluation-of-gaussian-molecular-integrals-4/
	// there is quite a bit of symbolic computation going on there, so I turned on to 

	// HSERILib: Gaussian integral evaluation
	// link: http://theory.rutgers.edu/~giese/notes/HSERILib.pdf
	// it discusses electron-electron integrals but also other kind of integrals very shortly
	// I took the information from there in order to implement this, but the resemblance
	// with the Mathematica Journal article exists, because both use the same recurrence relations

	class GaussianNuclear : public GaussianIntegral
	{
	public:
		Eigen::MatrixXd matrixCalc;

		void Reset(IntegralsRepository* repository, double alpha1, double alpha2, const Vector3D<double>& nucleus, const Vector3D<double>& center1, const Vector3D<double>& center2, unsigned int maxL1, unsigned int maxL2, bool calculateHorizontal = true);

		double operator()(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const { return getNuclear(QN1, QN2); }
		double getNuclear(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const;	
	protected:
		void VerticalRecursion(double alpha, const Vector3D<double> Rp, const Vector3D<double>& center1, const Vector3D<double>& difN, unsigned int maxL);
		void HorizontalRecursion(const Vector3D<double>& dif, unsigned int maxL1, unsigned int maxL2);
		friend class IntegralsRepository;
	};

}