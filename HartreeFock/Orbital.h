#pragma once

#include "Vector3D.h"

#include "QuantumNumbers.h"

namespace Orbitals {

	class OrbitalBase
	{
	public:
		virtual ~OrbitalBase() = default;

		virtual double operator()(const Vector3D<double>& r) const  noexcept = 0;

		virtual Vector3D<double> getGradient(const Vector3D<double>& r) const  noexcept = 0;
		virtual double getLaplacian(const Vector3D<double>& r) const noexcept = 0;
		virtual void Normalize() noexcept = 0;
	};

	class Orbital : public OrbitalBase
	{
	public:
		Vector3D<double> center;
		QuantumNumbers::QuantumNumbers angularMomentum;
		
		unsigned int ID = 0;
		unsigned int centerID = 0;
		unsigned int shellID = 0;

		const Vector3D<double>& getCenter() const noexcept
		{
			return center;
		}

		char AtomicOrbital() const noexcept { return angularMomentum.AtomicOrbital(); }
	};

}

