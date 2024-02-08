#include "stdafx.h"

#include "MathUtils.h"

//#define _USE_MATH_DEFINES // for C++
#include <cmath>

#include <algorithm>

#include "GaussianOrbital.h"


namespace Orbitals {

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
		const Vector3D R(r - center);

		return coeffProdNorm * pow(R.X, angularMomentum.l) * pow(R.Y, angularMomentum.m) *  pow(R.Z, angularMomentum.n) * exp(-alpha * R * R);
	}

	Vector3D<double> GaussianOrbital::getGradient(const Vector3D<double>& r) const
	{
		const Vector3D R(r - center);

		const double pRX = pow(R.X, angularMomentum.l);
		const double pRY = pow(R.Y, angularMomentum.m);
		const double pRZ = pow(R.Z, angularMomentum.n);
		const double powProd = pRX * pRY * pRZ;
		const double expaR2 = exp(-alpha * R * R);

		const Vector3D eD(-alpha * 2. * R);

		double valX = eD.X;
		if (angularMomentum.l > 0)
			valX += static_cast<double>(angularMomentum.l) / R.X;

		double valY = eD.Y;
		if (angularMomentum.m > 0)
			valY += static_cast<double>(angularMomentum.m) / R.Y;

		double valZ = eD.Z;
		if (angularMomentum.n > 0)
			valZ += static_cast<double>(angularMomentum.n) / R.Z;

		return coeffProdNorm * expaR2 * powProd * Vector3D(valX, valY, valZ);

		// Numerical, for tests:
		/*
		double h = 0.0001;

		double val = operator()(r);

		double dx = (operator()(Vector3D<double>(r.X + h, r.Y, r.Z)) - val) / h;
		double dy = (operator()(Vector3D<double>(r.X, r.Y + h, r.Z)) - val) / h;
		double dz = (operator()(Vector3D<double>(r.X, r.Y, r.Z + h)) - val) / h;

		return Vector3D<double>(dx, dy, dz);
		*/
	}

	double GaussianOrbital::getLaplacian(const Vector3D<double>& r) const
	{
		const Vector3D R(r - center);

		const double pRX = pow(R.X, angularMomentum.l);
		const double pRY = pow(R.Y, angularMomentum.m);
		const double pRZ = pow(R.Z, angularMomentum.n);
		const double powProd = pRX * pRY * pRZ;
		const double expaR2 = exp(-alpha * R * R);

		const double twoalpha = 2. * alpha;
		const Vector3D eD(-twoalpha * R);

		double valX = eD.X;
		if (angularMomentum.l > 0)
			valX += static_cast<double>(angularMomentum.l) / R.X;
		valX *= -R.X * twoalpha;

		double valY = eD.Y;
		if (angularMomentum.m > 0)
			valY += static_cast<double>(angularMomentum.m) / R.Y;
		valY *= -R.Y * twoalpha;

		double valZ = eD.Z;
		if (angularMomentum.n > 0)
			valZ += static_cast<double>(angularMomentum.n) / R.Z;
		valZ *= -R.Z * twoalpha;

		double valX2 = -twoalpha * (angularMomentum.l + 1.);
		if (angularMomentum.l > 1)
			valX2 += static_cast<double>(angularMomentum.l)* (angularMomentum.l - 1.) / (R.X * R.X);

		double valY2 = -twoalpha * (angularMomentum.m + 1.);
		if (angularMomentum.m > 1)
			valY2 += static_cast<double>(angularMomentum.m)* (angularMomentum.m - 1.) / (R.Y * R.Y);

		double valZ2 = -twoalpha * (angularMomentum.n + 1.);
		if (angularMomentum.n > 1)
			valZ2 += static_cast<double>(angularMomentum.n)* (angularMomentum.n - 1.) / (R.Z * R.Z);

		return (valX + valX2 + valY + valY2 + valZ + valZ2) * coeffProdNorm * powProd * expaR2;

		/*
		double h = 0.0001;
		double h2 = h * h;

		double val = operator()(r);

		double d2x = operator()(Vector3D<double>(r.X + h, r.Y, r.Z));
		double d2y = operator()(Vector3D<double>(r.X, r.Y + h, r.Z));
		double d2z = operator()(Vector3D<double>(r.X, r.Y, r.Z + h));

		d2x += operator()(Vector3D<double>(r.X - h, r.Y, r.Z)) - 2. * val;
		d2y += operator()(Vector3D<double>(r.X, r.Y - h, r.Z)) - 2. * val;
		d2z += operator()(Vector3D<double>(r.X, r.Y, r.Z - h)) - 2. * val;

		d2x /= h2;
		d2y /= h2;
		d2z /= h2;

		return d2x + d2y + d2z;
		*/
	}

	Vector3D<double> GaussianOrbital::ProductCenter(const GaussianOrbital& other) const
	{
		return (alpha * center + other.alpha * other.center) / (alpha + other.alpha);
	}


	void GaussianOrbital::Normalize()
	{
		normalizationFactor = getNormalizationFactor();
		coeffProdNorm = coefficient * normalizationFactor;
	}
}

