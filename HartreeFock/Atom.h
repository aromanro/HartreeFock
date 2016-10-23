#pragma once


#include "Vector3D.h"
#include "PrimitiveShell.h"

namespace Systems {

	class Atom
	{
	public:
		Vector3D<double> position;

		unsigned int Z;
		unsigned int electronsNumber;

		unsigned int ID;

		Atom(unsigned int nrZ = 0, int nrElectrons = -1);
		virtual ~Atom();
	};

	class AtomWithShells : public Atom
	{
	public:
		std::vector<Orbitals::ContractedGaussianShell> shells;

		AtomWithShells(unsigned int nrZ = 0, unsigned int nrElectrons = -1) : Atom(nrZ, nrElectrons) {}

		void AddShell(const std::string& name)
		{
			Orbitals::ContractedGaussianShell shell;

			for (auto c : name)	shell.AddOrbital(c);

			shells.push_back(std::move(shell));
		}

		void SetCenterForShells()
		{
			for (auto &shell : shells) shell.SetCenters(position);
		}

		void SetPosition(const Vector3D<double>& pos)
		{
			position = pos;
			SetCenterForShells();
		}

		unsigned int CountNumberOfContractedGaussians() const;
		unsigned int CountNumberOfGaussians() const;

		void GetMaxQN(double alpha, Orbitals::QuantumNumbers::QuantumNumbers& maxQN) const
		{
			maxQN.l = maxQN.m = maxQN.n = 0;

			for (const auto &shell : shells)
				for (const auto &orbital : shell.basisFunctions)
					for (const auto &gaussian : orbital.gaussianOrbitals)
						if (alpha == gaussian.alpha)
						{
							maxQN.l = max(maxQN.l, orbital.angularMomentum.l);
							maxQN.m = max(maxQN.m, orbital.angularMomentum.m);
							maxQN.n = max(maxQN.m, orbital.angularMomentum.n);
						}
		}

		unsigned int GetMaxAngularMomentum() const
		{
			unsigned int L = 0;

			for (const auto &shell : shells)
				for (const auto &orbital : shell.basisFunctions)
					L = max(L, orbital.angularMomentum);

			return L;
		}

		unsigned int GetMaxAngularMomentum(double alpha) const
		{
			unsigned int L = 0;

			for (const auto &shell : shells)
				for (const auto &orbital : shell.basisFunctions)
					for (const auto &gaussian : orbital.gaussianOrbitals)
						if (alpha == gaussian.alpha)
					         L = max(L, orbital.angularMomentum);

			return L;
		}

		void Normalize()
		{
			for (auto &shell : shells) shell.Normalize();
		}
	};

}