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

	Vector3D<double> ContractedGaussianOrbital::getGradient(const Vector3D<double>& r) const
	{
		Vector3D<double> res;

		for (const auto& orbital : gaussianOrbitals)
			res += orbital.getGradient(r);

		return res;
	}

	double ContractedGaussianOrbital::getLaplacian(const Vector3D<double>& r) const
	{
		double res = 0;

		for (const auto& orbital : gaussianOrbitals)
			res += orbital.getLaplacian(r);

		return res;
	}

	void ContractedGaussianOrbital::Normalize()
	{
		for (auto& orb : gaussianOrbitals) orb.Normalize();
	}
}

