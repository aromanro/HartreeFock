#pragma once
class Test
{
public:
	static void OutputMatricesForAtom(const std::string& atomName, const std::string& basisSetName, const std::string& fileName);
	static void OutputMatrices(Systems::Molecule& molecule, std::ofstream& file);
};

