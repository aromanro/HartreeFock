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

		virtual double operator()(const Vector3D<double>& r) const override;

		ContractedGaussianOrbital();
		virtual ~ContractedGaussianOrbital();
		void Normalize();
	};

}