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
		Matrices::MomentMatrix momentMatrix;
		Matrices::KineticMatrix kineticMatrix;
		Matrices::NuclearMatrix nuclearMatrix;

	public:
		Eigen::MatrixXd h;

	protected:
		double nuclearRepulsionEnergy;
		// for external electric field
		// it's not added up to total energy (for now, at least) because the derivative of energy with respect to the electric field is used to obtain the dipole moment
		// there is a way to analytically compute the nuclear dipole moment instead of relying on numerical derivatives, so that's used instead
		// if this energy is added up (and the analytically computed nuclear dipole moment is not), basically the result with numerical derivative is the same, so I prefer the analytic result
		// see the Test class for details
		double nuclearElectricFieldEnergy; 

		Eigen::MatrixXd U;
		Eigen::MatrixXd s;
		Eigen::MatrixXd V;
		Eigen::MatrixXd Vt;

		int maxIterations;

		bool inited;

	public:
		int numberOfOrbitals;

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

		double GetNuclearElectricFieldEnergy() const
		{
			return nuclearElectricFieldEnergy;			
		}

		virtual double GetMP2Energy() const
		{
			return mp2Energy;
		}

		virtual double CalculateMp2Energy() = 0;

		virtual double CalculateAtomicCharge(int atom) const = 0;

		Vector3D<double> GetNuclearMoment() const;

		virtual Vector3D<double> GetMoment() const = 0;

	protected:
		static double DiffDensityMatrices(const Eigen::MatrixXd& oldP, const Eigen::MatrixXd& newP);
		void NormalizeC(Eigen::MatrixXd& C, const std::vector<bool>& occupied);

		void FirstIterations(int& iter, double& curEnergy, double& prevEnergy);
		void SelfConsistentIterations(int& iter, double& curEnergy, double& prevEnergy);
		void NormalIterations(int& iter, double& curEnergy, double& prevEnergy);
	};


}