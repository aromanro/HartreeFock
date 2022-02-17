#pragma once

#include "Orbital.h"

namespace Orbitals {

	class GaussianOrbital : public Orbital
	{
	public:
		double coefficient;
		double alpha;

		double normalizationFactor;

		GaussianOrbital();
		virtual ~GaussianOrbital();

		virtual double getCoefficient() const;
		virtual double getAlpha() const;

		double operator()(const Vector3D<double>& r) const override;
		Vector3D<double> getGradient(const Vector3D<double>& r) const override;
		double getLaplacian(const Vector3D<double>& r) const override;

		Vector3D<double> ProductCenter(const GaussianOrbital& other) const;

	protected:
		double getNormalizationFactor() const;

		double coeffProdNorm;

	public:
		void Normalize();
	};
}


