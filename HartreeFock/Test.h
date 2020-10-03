#pragma once

#include "ChemUtils.h"
#include "Basis.h"

class Test
{
public:
	Test();

	static void OutputMatricesForAtom(const std::string& atomName, const std::string& basisSetName, const std::string& fileName);
	static void OutputMatrices(Systems::Molecule& molecule, std::ofstream& file);

	void TestWater(const std::string& fileName);

protected:
	Chemistry::Basis basisSTO3G;
};

