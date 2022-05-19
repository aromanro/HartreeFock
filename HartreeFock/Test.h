#pragma once

#include "ChemUtils.h"
#include "Basis.h"

namespace HartreeFock {
	class HartreeFockAlgorithm;
}


class Test
{
public:
	Test(const std::string& basisFile = "sto3g.txt");

	static void OutputMatricesForAtom(const std::string& atomName, const std::string& basisSetName, const std::string& fileName);
	static void OutputMatrices(Systems::Molecule& molecule, std::ofstream& file, const std::string& sfileName = "", const std::string& tfileName = "", const std::string& vfileName = "", const std::string& erifileName = "", bool useDIIS = false, const double expectedNucEnergy = 0);

	static void CheckDifferences(const Eigen::MatrixXd& matrix, const std::string& matrixFileName, std::ofstream& file);
	static void CheckDifferences(const GaussianIntegrals::IntegralsRepository& repo, const std::string& eriFileName, std::ofstream& file);

	void TestWater(const std::string& fileName, const std::string& sfileName = "", const std::string& tfileName = "", const std::string& vfileName = "", const std::string& erifileName = "", bool useDIIS = false);
	void TestMethane(const std::string& fileName, const std::string& sfileName = "", const std::string& tfileName = "", const std::string& vfileName = "", const std::string& erifileName = "", bool useDIIS = false);

	// specifically for a single test now, maybe I'll generalize it later
	void Test::TestWaterDipoleMoment(const std::string& fileName);

protected:
	static void OutputMatrix(const Eigen::MatrixXd& matrix, std::ofstream& file);

	static void TestIterationsAndPostHF(Systems::Molecule& molecule, HartreeFock::HartreeFockAlgorithm* hartreeFock, const Eigen::MatrixXd& h, Eigen::MatrixXd& FockMatrix, Eigen::MatrixXd& FockTransformed, Eigen::MatrixXd& C, Eigen::MatrixXd& DensityMatrix, std::ofstream& file, bool useDIIS);

	static void Iterate(double& rmsD, double& deltaE, Systems::Molecule& molecule, HartreeFock::HartreeFockAlgorithm* hartreeFock, const Eigen::MatrixXd& h, Eigen::MatrixXd& FockMatrix, Eigen::MatrixXd& FockTransformed, Eigen::MatrixXd& C, Eigen::MatrixXd& DensityMatrix, std::ofstream& file, bool useDIIS);

	// it seems that I might need it for more tests
	// for now, I intend to test the dipole moment computation (the newer method, see pages 376-378 in Szabo & Ostlund book)
	void InitWaterMolecule(Systems::Molecule& molecule);

	Chemistry::Basis basis;
};

