#pragma once

#include "Orbital.h"

namespace Orbitals {

	class GaussianOrbital : public Orbital
	{
	public:
		double coefficient = 1;
		double alpha = 1;

		double normalizationFactor = 1;

		virtual double getCoefficient() const noexcept;
		virtual double getAlpha() const noexcept;

		double operator()(const Vector3D<double>& r) const noexcept override;
		Vector3D<double> getGradient(const Vector3D<double>& r) const noexcept override;
		double getLaplacian(const Vector3D<double>& r) const noexcept override;

		Vector3D<double> ProductCenter(const GaussianOrbital& other) const noexcept;

		void Normalize() noexcept override;

	protected:
		double getNormalizationFactor() const noexcept;

		double coeffProdNorm = 0;
	};
}


