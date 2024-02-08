#include "stdafx.h"
#include "Molecule.h"
#include "ChemUtils.h"

#include <algorithm>
#include <numeric>
#include <functional>


#include <sstream>

#include <fstream>
#include <iostream>
#include <iomanip>

namespace Systems {


	unsigned int Molecule::CountNumberOfContractedGaussians() const
	{
		unsigned int res = 0;

		for (const auto& atom : atoms)
			res += atom.CountNumberOfContractedGaussians();

		return res;
	}

	unsigned int Molecule::CountNumberOfGaussians() const
	{
		unsigned int res = 0;

		for (const auto& atom : atoms)
			res += atom.CountNumberOfGaussians();

		return res;
	}

}

void Systems::Molecule::SetIDs()
{
	unsigned int ID = 0;
	unsigned int contractedID = 0;
	unsigned int shellID = 0;
	unsigned int atomID = 0;

	for (auto& atom : atoms)
	{
		atom.ID = atomID++;
		
		for (auto& shell : atom.shells)
		{
			shell.ID = shellID++;
			shell.centerID = atom.ID;

			for (auto& orbital : shell.basisFunctions)
			{
				orbital.ID = contractedID++;
				orbital.centerID = atom.ID;
				orbital.shellID = shell.ID;

				for (auto& gaussian : orbital.gaussianOrbitals)
				{
					gaussian.ID = ID++;
					gaussian.centerID = atom.ID;
					gaussian.shellID = shell.ID;
				}
			}
		}
	}
}


double Systems::Molecule::NuclearRepulsionEnergy() const
{
	double energy = 0;

	for (unsigned int atom1 = 0; atom1 < atoms.size(); ++atom1)
		for (unsigned int atom2 = atom1 + 1; atom2 < atoms.size(); ++atom2)
			energy += static_cast<double>(atoms[atom1].Z) * atoms[atom2].Z / (atoms[atom1].position - atoms[atom2].position).Length();

	return energy;
}

double Systems::Molecule::NuclearElectricFieldEnergy() const
{
	double energy = 0;

	for (unsigned int atom = 0; atom < atoms.size(); ++atom)
		energy -= static_cast<double>(atoms[atom].Z) * atoms[atom].position * ElectricField;

	return energy;
}

unsigned int Systems::Molecule::ElectronsNumber() const
{
	if (alphaElectrons > 0 || betaElectrons > 0) return alphaElectrons + betaElectrons;

	return static_cast<unsigned int>(std::accumulate(atoms.begin(), atoms.end(), 0, [](unsigned int result, const AtomWithShells& atom) { return result + atom.electronsNumber; }));
}


unsigned int Systems::Molecule::GetMaxAngularMomentum() const
{
	unsigned int L = 0;

	for (const auto &atom : atoms)
			L = max(L, atom.GetMaxAngularMomentum());

	return L;
}


void Systems::Molecule::SetCenterForShells()
{
	for (auto &atom : atoms) atom.SetCenterForShells();
}


void Systems::Molecule::Normalize()
{
	for (auto &atom : atoms) atom.Normalize();
}


void Systems::Molecule::Init()
{
	SetCenterForShells();
	SetIDs();

	const unsigned int elNum = ElectronsNumber();

	// not set yet, set them
	if (0 == alphaElectrons && 0 == betaElectrons)
	{
		alphaElectrons = static_cast<unsigned int>(elNum / 2);
		betaElectrons = elNum - alphaElectrons;
	}
}


bool Systems::Molecule::LoadXYZ(const std::string& fileName, const Chemistry::Basis& basis)
{
	try
	{
		std::ifstream mfile(fileName);
		if (!mfile) return false;

		std::string line;

		if (!std::getline(mfile, line)) return false;

		unsigned int nrAtoms = std::stoi(line);

		// read comment
		if (!std::getline(mfile, line)) return false;
		for (unsigned int i = 0; i < nrAtoms; ++i)
		{
			if (!std::getline(mfile, line)) return false; // atom i
			std::istringstream lineStream(line);

			std::string atomName;
			double x, y, z;

			lineStream >> atomName >> x >> y >> z;

			unsigned int Z = Chemistry::ChemUtils::GetZForAtom(atomName);
			if (0 == Z) return false;

			bool found = false;
			for (const auto& atom : basis.atoms)
			{
				if (atom.Z == Z)
				{
					found = true;
					atoms.push_back(atom);

					atoms.back().position.X = x;
					atoms.back().position.Y = y;
					atoms.back().position.Z = z;

					break;
				}
			}

			if (!found) return false;
		}

		Init();

		return true;
	}
	catch (...)
	{
	}

	return false;
}
