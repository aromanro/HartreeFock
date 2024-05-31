#include "stdafx.h"
#include "ContractedGaussianOrbital.h"

namespace Orbitals {

	double ContractedGaussianOrbital::operator()(const Vector3D<double>& r) const noexcept
	{
		double res = 0;

		for (const auto &orbital : gaussianOrbitals)
			res += orbital(r);

		return res;
	}

	Vector3D<double> ContractedGaussianOrbital::getGradient(const Vector3D<double>& r) const noexcept
	{
		Vector3D<double> res;

		for (const auto& orbital : gaussianOrbitals)
			res += orbital.getGradient(r);

		return res;
	}

	double ContractedGaussianOrbital::getLaplacian(const Vector3D<double>& r) const noexcept
	{
		double res = 0;

		for (const auto& orbital : gaussianOrbitals)
			res += orbital.getLaplacian(r);

		return res;
	}

	void ContractedGaussianOrbital::Normalize() noexcept
	{
		for (auto& orb : gaussianOrbitals) orb.Normalize();
	}
}

