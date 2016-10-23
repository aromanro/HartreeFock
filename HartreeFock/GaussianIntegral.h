#pragma once

#include "Vector3D.h"
#include "GaussianOrbital.h"

namespace GaussianIntegrals {

	class GaussianIntegral
	{
	public:
		GaussianIntegral();
		virtual ~GaussianIntegral();

		static Vector3D<double> ProductCenter(const Orbitals::GaussianOrbital& orbital1, const Orbitals::GaussianOrbital& orbital2)
		{
			return orbital1.ProductCenter(orbital2);
		}
	};


}