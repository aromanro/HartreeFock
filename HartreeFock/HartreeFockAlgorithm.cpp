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
		momentMatrix.SetRepository(&integralsRepository);
		kineticMatrix.SetRepository(&integralsRepository);
		nuclearMatrix.SetRepository(&integralsRepository);

		overlapMatrix.Calculate();
		momentMatrix.Calculate();
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

		//V = U * s;
		V = U * s * U.adjoint(); // This is S^-1/2 in the original basis

		Vt = V.adjoint(); // V is actually symmetric, so this is unnecessary, Vt = V would do (not in the case V = U * s, though)

		// now, a trick that might get rid of some numerical errors from the solver, for the symmetric case (comment it if using the other one)
		V = Vt = (V + Vt) * 0.5; // now it's really symmetrical

		// what will happen is that we'll work in a 'molecular orbitals' basis, not in the original 'atomic orbitals' one
		// the new basis is orthonormal, that is, the overlap is I, the identity matrix

		// the transformed S (overlap matrix) in the new 'molecular orbitals' basis is
		// S^-1/2 * S * S^1/2 = I
		// so the generalized eigenvalue problem becomes the regular eigenvalue problem
		// see the RestrictedHartreeFock::Step for how it's done (with comments)
		// the Test class allows you to verify this

		// U is obviously unitary by construction 
		// it's the matrix with eigenvectors as columns, multiplying U * U^t just gets I
		// the matrix elements being the scalar product of the vectors, which makes it is obviously zero everywhere (different vectors are orthogonal) except the diagonal, which gets 1, that is, the norm, the scalar product of a vector with itself
		// due of the orthonormality
		
		// if you would do the transform with U^t * S * U, you would simply diagonalize S, obviously having the eigenvalues on the diagonal
		// we want more than that, we want to get I, so we multiply left and right with s (S^-1/2 in the diagonal basis)
		// we end up with a multiplication of three matrices, left and right having 1/sqrt(eigenvalue) on the diagonal and the middle one with eigenvalues on the diagonal
		// obviously that gets 1 on the diagonal and 0 elsewhere

		// I = s * U^t * S * U * s
		// we could simply call V = U * s and work like that 
		
		// a not so old version of the program worked with that, but in order to have intermediate results comparable with some results given here:
		// https://github.com/CrawfordGroup/ProgrammingProjects
		// I switched to the current variant, which also has some other nice property, see below:
		
		// since the identity is diagonal 1 in any basis, we can transform it back to original basis
		// I = U * I * U^t = U * s * U^t * S * U * s * U^t
		// and now we have I = Vt * S * V
		// where V = U * s * U^t and Vt = (U * s * U^t)^t = (U^t)^t * s^t * U^t = U * s * U^t = V
		// so now V should be actually symmetric (V = Vt), check it out with the Test class

		// for the case V = U * s though, things are not so nice (but can still work fine):
		// I = s * U^t * S * U * s = Vt * S * V
		// V = U * s and Vt = (U * s)^t = s * U^t which is not the same as V 

		// Note that obviously the 'molecular orbitals' basis would be different if using V = U * s
		// but still orthonormal and things should go fine in that one, too
		
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


	Vector3D<double> HartreeFockAlgorithm::GetNuclearMoment() const
	{
		Vector3D<double> moment;

		for (const Systems::Atom& atom : integralsRepository.getMolecule()->atoms)
			moment += static_cast<double>(atom.Z) * atom.position;

		return moment;
	}


}








