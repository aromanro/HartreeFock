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
		friend class Test;
	protected:
		double totalEnergy;
		double mp2Energy;

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

		double lastErrorEst;

		bool UseDIIS;
		int maxDIISiterations;

		int normalIterAfterDIIS;

		HartreeFockAlgorithm(int iterations = 3000);
		virtual ~HartreeFockAlgorithm();
		
		virtual void Init(Systems::Molecule* molecule);

		double Calculate();
		
		virtual double Step(int iter) = 0;
		
		double GetTotalEnergy() const
		{
			return totalEnergy;
		};

		virtual double GetMP2Energy() const
		{
			return mp2Energy;
		}

		virtual double CalculateMp2Energy() = 0;

	protected:
		static double DiffDensityMatrices(const Eigen::MatrixXd& oldP, const Eigen::MatrixXd& newP);
		void NormalizeC(Eigen::MatrixXd& C, const std::vector<bool>& occupied);
	};


}