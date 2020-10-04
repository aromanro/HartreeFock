#include "stdafx.h"
#include "HartreeFockAlgorithm.h"

#include "Constants.h"


namespace HartreeFock {

	HartreeFockAlgorithm::HartreeFockAlgorithm(int iterations)
		: totalEnergy(std::numeric_limits<double>::infinity()), mp2Energy(0), nuclearRepulsionEnergy(0), numberOfOrbitals(0),  maxIterations(iterations), inited(false), alpha(0.75), initGuess(0.75), terminate(false), converged(false),
		HOMOEnergy(0), lastErrorEst(0), UseDIIS(true), maxDIISiterations(1000), normalIterAfterDIIS(500)
	{
	}


	HartreeFockAlgorithm::~HartreeFockAlgorithm()
	{
	}


	void HartreeFockAlgorithm::Init(Systems::Molecule* molecule)
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

		V = U * s * U.adjoint(); // This is S^-1/2

		Vt = V.adjoint();

		inited = true;
	}



	double HartreeFockAlgorithm::Calculate()
	{
		double curEnergy = 0;
		double prevEnergy = std::numeric_limits<double>::infinity();

		if (!inited) return prevEnergy;

		int iter = 0;
		for (; iter < maxIterations; ++iter)
		{
			const double rmsD = Step(iter);

			curEnergy = GetTotalEnergy();

			if (abs(prevEnergy - curEnergy) <= energyConvergence && rmsD < rmsDConvergence && lastErrorEst < diisConvergence) {
				converged = true;
				break;
			}

			if (terminate) return curEnergy;

			prevEnergy = curEnergy;
		}

		// did it converge with DIIS?
		if (UseDIIS && iter < maxDIISiterations && converged && normalIterAfterDIIS)
		{
			UseDIIS = false;

			// do a loop without checking the convergence, sometimes DIIS gets stuck in a bad position close to the minimum
			for (int i = 0; i < normalIterAfterDIIS; ++i)
			{
				curEnergy = Step(iter + i);
				if (terminate)
				{
					UseDIIS = true; // restore it back
					return curEnergy;
				}
			}

			iter += normalIterAfterDIIS;

			// now continue with normal iteration with convergence checking
			for (; iter < maxIterations; ++iter)
			{
				const double rmsD = Step(iter);

				curEnergy = GetTotalEnergy();

				if (abs(prevEnergy - curEnergy) <= energyConvergence && rmsDConvergence < rmsDConvergence) {
					converged = true;
					UseDIIS = true; // restore it back
					return curEnergy;
				}

				if (terminate)
				{
					UseDIIS = true; // restore it back
					return curEnergy;
				}

				prevEnergy = curEnergy;
			}

			UseDIIS = true; // restore it back
		}

		return curEnergy;
	}

	double HartreeFockAlgorithm::DiffDensityMatrices(const Eigen::MatrixXd& oldP, const Eigen::MatrixXd& newP)
	{
		double res = 0.0;

		assert(oldP.cols() == newP.cols());
		assert(oldP.rows() == newP.rows());

		for (int i = 0; i < oldP.rows(); ++i)
			for (int j = 0; j < oldP.cols(); ++j)
			{
				const double val = oldP(i, j) - newP(i, j);
			
				res += val * val;
			}

		return sqrt(res);
	}


	void HartreeFockAlgorithm::NormalizeC(Eigen::MatrixXd& C, const std::vector<bool>& occupied)
	{
		assert(occupied.size() <= C.cols());
		assert(C.rows() == overlapMatrix.matrix.rows());

		for(int vec = 0; vec < occupied.size(); ++vec) 
		{
			if (!occupied[vec]) continue;

			double factor = 0.0;

			for(int i = 0; i < overlapMatrix.matrix.rows(); ++i)
				for(int j = 0; j < overlapMatrix.matrix.cols(); ++j)
					factor += C(i, vec) * overlapMatrix.matrix(i, j) * C(j, vec);

			C.col(vec) /= sqrt(factor);
		}
	}

}








