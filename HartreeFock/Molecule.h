#pragma once

#include <vector>
#include "Atom.h"

#include "Basis.h"

namespace Systems {


	class Molecule
	{
	public:
		std::vector<AtomWithShells> atoms;

		unsigned int alphaElectrons = 0;
		unsigned int betaElectrons = 0;

		Vector3D<double> ElectricField;

		unsigned int CountNumberOfContractedGaussians() const;
		unsigned int CountNumberOfGaussians() const;
		void SetIDs();
		double NuclearRepulsionEnergy() const;
		double NuclearElectricFieldEnergy() const;
		unsigned int ElectronsNumber() const;
		unsigned int GetMaxAngularMomentum() const;
		void SetCenterForShells();
		void Normalize();
		void Init();

		bool LoadXYZ(const std::string& fileName, const Chemistry::Basis& basis);
	};

}

