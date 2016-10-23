#include "stdafx.h"
#include "PrimitiveShell.h"


namespace Orbitals {

	PrimitiveGaussianShell::PrimitiveGaussianShell()
	{
	}


	PrimitiveGaussianShell::~PrimitiveGaussianShell()
	{
	}

	Vector3D<double> PrimitiveGaussianShell::getCenter() const {
		assert(basisFunctions.size());

		return basisFunctions.front().center;
	}

	double PrimitiveGaussianShell::getAlpha() const {
		assert(basisFunctions.size());

		return basisFunctions.front().alpha;
	}

	double PrimitiveGaussianShell::operator()(const Vector3D<double>& r) const
	{
		double res = 0;

		for (const auto& orbital : basisFunctions)
			res += orbital(r);

		return res;
	}

	ContractedGaussianShell::ContractedGaussianShell()
	{
	}

	ContractedGaussianShell::~ContractedGaussianShell()
	{
	}

	Vector3D<double> ContractedGaussianShell::getCenter() const {
		assert(basisFunctions.size());

		return basisFunctions.front().center;
	}

	std::string ContractedGaussianShell::GetShellString() const
	{
		assert(basisFunctions.size());

		std::string res;

		for (auto &orbital : basisFunctions)
		{
			char c = (char)toupper(orbital.AtomicOrbital());

			if (std::string::npos == res.find(c))
			{
				res += c;
			}
		}

		return res;
	}

	unsigned int ContractedGaussianShell::AdjustOrbitalsCount(char orbital, unsigned int res)
	{
		switch (tolower(orbital))
		{
		default:
		case 's':
			break;
		case 'p':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(1);
			break;
		case 'd':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(2);
			break;
		case 'f':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(3);
			break;
		case 'g':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(4);
			break;
		case 'h':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(5);
			break;
		}

		return res;
	}

	unsigned int ContractedGaussianShell::CountOrbitals(char orbitalChar) const
	{
		unsigned int res = 0;

		for (auto const &orbital : basisFunctions)
		{
			for (auto const &gaussian : orbital.gaussianOrbitals)
				if (orbitalChar == gaussian.AtomicOrbital())
					++res;
		}

		return AdjustOrbitalsCount(orbitalChar, res);
	}

	unsigned int ContractedGaussianShell::CountContractedOrbitals(char orbitalChar) const
	{
		unsigned int res = 0;

		for (auto const &orbital : basisFunctions)
		{
			if (orbitalChar == orbital.gaussianOrbitals.front().AtomicOrbital())
				++res;
		}

		return AdjustOrbitalsCount(orbitalChar, res);
	}

	unsigned int ContractedGaussianShell::CountNumberOfContractedGaussians() const
	{		
		return (unsigned int)basisFunctions.size();
	}

	unsigned int ContractedGaussianShell::CountNumberOfGaussians() const
	{
		unsigned int res = 0;

		for (auto const &orbital : basisFunctions)
			res += (unsigned int)orbital.gaussianOrbitals.size();

		return res;
	}

	double ContractedGaussianShell::operator()(const Vector3D<double>& r) const
	{
		double res = 0;

		for (const auto& orbital : basisFunctions)
			res += orbital(r);

		return res;
	}

}

void Orbitals::ContractedGaussianShell::AddOrbital(char type)
{
	unsigned int L = 0;

	switch (tolower(type))
	{
	case 's':
		L = 0;
		break;
	case 'p':
		L = 1;
		break;
	case 'd':
		L = 2;
		break;
	case 'f':
		L = 3;
		break;
	case 'g':
		L = 4;
		break;
	case 'h':
		L = 5;
		break;
	default:
		AfxMessageBox(L"Unknown orbital type");
		break;
	}

	AddOrbitalsInCanonicalOrder(L);
}



void Orbitals::ContractedGaussianShell::AddGaussians(double exponent)
{
	for (auto &orbital : basisFunctions)
	{
		GaussianOrbital gaussian;

		gaussian.angularMomentum = orbital.angularMomentum;
		gaussian.alpha = exponent;
		
		orbital.gaussianOrbitals.push_back(std::move(gaussian));
	}
}


void Orbitals::ContractedGaussianShell::SetCenters(const Vector3D<double>& center)
{
	for (auto &orbital : basisFunctions)
	{
		orbital.center = center;

		for (auto& gaussian : orbital.gaussianOrbitals)
			gaussian.center = center;
	}
}


void Orbitals::ContractedGaussianShell::AddOrbitalsInCanonicalOrder(unsigned int L)
{
	Orbitals::ContractedGaussianOrbital orbital;
	
	orbital.angularMomentum.l = L;

	while (orbital.angularMomentum == L)
	{
		basisFunctions.push_back(orbital);

		++orbital.angularMomentum;
	}
}


void Orbitals::ContractedGaussianShell::Normalize()
{
	for (auto& orb : basisFunctions) orb.Normalize();
}


void Orbitals::PrimitiveGaussianShell::Normalize()
{
	for (auto& orb : basisFunctions) orb.Normalize();
}