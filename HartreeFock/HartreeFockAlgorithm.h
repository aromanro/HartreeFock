#pragma once

#include "MathUtils.h"

#include "GaussianOrbital.h"
#include "Molecule.h"
#include "IntegralsRepository.h"
#include "QuantumMatrix.h"

#include "BoysFunction.h"


#include <atomic>

namespace HartreeFock {

	class HartreeFockAlgorithm
	{
	protected:
		
		Matrices::OverlapMatrix overlapMatrix;
		Matrices::KineticMatrix kineticMatrix;
		Matrices::NuclearMatrix nuclearMatrix;
		
		Eigen::MatrixXd h;

		double nuclearRepulsionEnergy;

		Eigen::MatrixXd U;
		Eigen::MatrixXd s;
		Eigen::MatrixXd V;
		Eigen::MatrixXd Vt;

		int numberOfOrbitals;

		int maxIterations;

		bool inited;

	public:
		GaussianIntegrals::IntegralsRepository integralsRepository;


		double alpha;

		double initGuess;

		std::atomic_bool terminate;

		bool converged;

		double HOMOEnergy;

		HartreeFockAlgorithm(int iterations = 3000);
		virtual ~HartreeFockAlgorithm();
		
		virtual void Init(Systems::Molecule* molecule);

		double Calculate();
		
		virtual void Step(int iter) = 0;
		virtual double GetTotalEnergy() const = 0;
	protected:
		static double DiffDensityMatrices(const Eigen::MatrixXd& oldP, const Eigen::MatrixXd& newP);
		void NormalizeC(Eigen::MatrixXd& C, int nrOccupiedLevels);
	};


}