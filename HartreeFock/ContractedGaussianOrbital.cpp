#include "stdafx.h"
#include "ContractedGaussianOrbital.h"

namespace Orbitals {

	ContractedGaussianOrbital::ContractedGaussianOrbital()
	{
	}


	ContractedGaussianOrbital::~ContractedGaussianOrbital()
	{
	}

	double ContractedGaussianOrbital::operator()(const Vector3D<double>& r) const
	{
		double res = 0;

		for (const auto &orbital : gaussianOrbitals)
			res += orbital(r);

		return res;
	}

}

void Orbitals::ContractedGaussianOrbital::Normalize()
{
	for (auto& orb : gaussianOrbitals) orb.Normalize();
}
