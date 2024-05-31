#pragma once
#include "GaussianOrbital.h"


#include <vector>

namespace Orbitals {

	// the contained gaussian orbitals all have the same center and quantum numbers, whence the derivation from Orbital
	class ContractedGaussianOrbital :
		public Orbital
	{
	public:
		std::vector<GaussianOrbital> gaussianOrbitals;

		double operator()(const Vector3D<double>& r) const noexcept override;
		Vector3D<double> getGradient(const Vector3D<double>& r) const noexcept override;
		double getLaplacian(const Vector3D<double>& r) const noexcept override;

		void Normalize() noexcept override;

		void AddGaussian(double exponent)
		{
			GaussianOrbital gaussian;

			gaussian.angularMomentum = angularMomentum;
			gaussian.alpha = exponent;

			gaussianOrbitals.emplace_back(std::move(gaussian));
		}
	};

}


