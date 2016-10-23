#pragma once

#include <vector>
#include "Atom.h"

namespace Systems {


	class Molecule
	{
	public:
		std::vector<AtomWithShells> atoms;

		unsigned int alphaElectrons;
		unsigned int betaElectrons;


		Molecule();
		~Molecule();

		unsigned int CountNumberOfContractedGaussians() const;
		unsigned int CountNumberOfGaussians() const;
		void SetIDs();
		double NuclearRepulsionEnergy() const;
		unsigned int ElectronsNumber();
		unsigned int GetMaxAngularMomentum();
		void SetCenterForShells();
		void Normalize();
		void Init();
	};

}

