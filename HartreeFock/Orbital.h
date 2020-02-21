#pragma once

#include "Vector3D.h"

#include "QuantumNumbers.h"

namespace Orbitals {

	class OrbitalBase
	{
	public:
		OrbitalBase() {}
		virtual ~OrbitalBase() {}

		virtual double operator()(const Vector3D<double>& r) const = 0;

		virtual Vector3D<double> getGradient(const Vector3D<double>& r) const = 0;
		virtual double getLaplacian(const Vector3D<double>& r) const = 0;
	};

	class Orbital : public OrbitalBase
	{
	public:
		Vector3D<double> center;
		QuantumNumbers::QuantumNumbers angularMomentum;
		
		unsigned int ID;
		unsigned int centerID;
		unsigned int shellID;

		Orbital();
		virtual ~Orbital();

		virtual Vector3D<double> getCenter() const;
		char AtomicOrbital() const { return angularMomentum.AtomicOrbital(); }
	};

}

