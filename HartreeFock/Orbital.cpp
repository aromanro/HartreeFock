#include "stdafx.h"
#include "Orbital.h"

namespace Orbitals {

	Orbital::Orbital()
		: ID(0), centerID(0), shellID(0), angularMomentum(0, 0, 0)
	{
	}


	Orbital::~Orbital()
	{
	}


	Vector3D<double> Orbital::getCenter() const
	{
		return center;
	}

}

