#include "stdafx.h"
#include "HartreeFockAlgorithm.h"


namespace HartreeFock {

	HartreeFockAlgorithm::HartreeFockAlgorithm(int iterations)
		: alpha(0.75), nuclearRepulsionEnergy(0), numberOfOrbitals(0), initGuess(0.75), maxIterations(iterations), terminate(false), converged(false), inited(false)
	{
	}


	HartreeFockAlgorithm::~HartreeFockAlgorithm()
	{
	}

}

void HartreeFock::HartreeFockAlgorithm::Init(Systems::Molecule* molecule)
{
	converged = false;
	integralsRepository.Reset(molecule);

	overlapMatrix.SetRepository(&integralsRepository);
	kineticMatrix.SetRepository(&integralsRepository);
	nuclearMatrix.SetRepository(&integralsRepository);

	overlapMatrix.Calculate();
	kineticMatrix.Calculate();
	nuclearMatrix.Calculate();

	integralsRepository.ClearMatricesMaps();

	h = kineticMatrix.matrix + nuclearMatrix.matrix;

	nuclearRepulsionEnergy = molecule->NuclearRepulsionEnergy();

	numberOfOrbitals = molecule->CountNumberOfContractedGaussians();

	integralsRepository.CalculateElectronElectronIntegrals();
	integralsRepository.ClearAllMaps();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(overlapMatrix.matrix);

	U = solver.eigenvectors();
	s = solver.eigenvalues().cwiseInverse().cwiseSqrt().asDiagonal();

	V = U * s;

	Vt = V.adjoint();

	inited = true;
}

double HartreeFock::HartreeFockAlgorithm::Calculate()
{
	double curEnergy = 0;
	double prevEnergy = std::numeric_limits<double>::infinity();

	if (!inited) return prevEnergy;

	// some big number before bail out
	for (int iter = 0; iter < maxIterations; ++iter)
	{
		Step(iter);

		curEnergy = GetTotalEnergy();

		if (abs(prevEnergy - curEnergy) <= 1E-13) {
			converged = true;
			break;
		}

		if (terminate) break;

		prevEnergy = curEnergy;
	}

	return curEnergy;
}




double HartreeFock::HartreeFockAlgorithm::DiffDensityMatrices(const Eigen::MatrixXd& oldP, const Eigen::MatrixXd& newP)
{
	double res = 0.0;

	assert(oldP.cols() == newP.cols());
	assert(oldP.rows() == newP.rows());

	for (int i = 0; i < oldP.rows(); ++i)
		for (int j = 0; j < oldP.cols(); ++j)
		{
			double val = oldP(i, j) - newP(i, j);
			
			res += val * val;
		}

	return sqrt(res);
}
