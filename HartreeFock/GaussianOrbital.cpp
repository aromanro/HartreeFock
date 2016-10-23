#include "stdafx.h"

#include "MathUtils.h"

//#define _USE_MATH_DEFINES // for C++
#include <cmath>

#include <algorithm>

#include "GaussianOrbital.h"


namespace Orbitals {

	GaussianOrbital::GaussianOrbital()
		: coefficient(1), alpha(1), normalizationFactor(1)
	{
	}


	GaussianOrbital::~GaussianOrbital()
	{
	}


	double GaussianOrbital::getCoefficient() const
	{
		return coefficient;
	}

	double GaussianOrbital::getAlpha() const
	{
		return alpha;
	}

	double GaussianOrbital::getNormalizationFactor() const
	{
		return pow(2. * alpha / M_PI, 3. / 4.) *
			pow (4. * alpha, angularMomentum.AngularMomentum() / 2.) /
			sqrt(GaussianIntegrals::MathUtils::DoubleFactorial(2 * angularMomentum.l - 1) * GaussianIntegrals::MathUtils::DoubleFactorial(2 * angularMomentum.m - 1) * GaussianIntegrals::MathUtils::DoubleFactorial(2 * angularMomentum.n - 1));
	}

	double GaussianOrbital::operator()(const Vector3D<double>& r) const
	{
		Vector3D<double> R = r - center;

		return coefficient * normalizationFactor * pow(R.X, angularMomentum.l) * pow(R.Y, angularMomentum.m) *  pow(R.Z, angularMomentum.n) * exp(-alpha * R * R);
	}

	Vector3D<double> GaussianOrbital::ProductCenter(const GaussianOrbital& other) const
	{
		return (alpha * center + other.alpha * other.center) / (alpha + other.alpha);
	}


	void GaussianOrbital::Normalize()
	{
		normalizationFactor = getNormalizationFactor();
	}
}

