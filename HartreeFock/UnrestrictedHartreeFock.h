#pragma once
#include "HartreeFockAlgorithm.h"

#include "DIIS.h"
#include <functional>

namespace HartreeFock {

	class UnrestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
		friend class Test;
	public:
		Eigen::MatrixXd DensityMatrixPlus;
		Eigen::MatrixXd DensityMatrixMinus;

		unsigned int nrOccupiedLevelsPlus;
		unsigned int nrOccupiedLevelsMinus;

		// for now in the program it will be filled up with true up to 'nrOccupiedLevels'
		// could be used to compute excited levels, just enlarge it after init and set to true the occupied levels and false the ones that are not occupied
		// adjust the nrOccupied value above then to get the proper homo energy
		std::vector<bool> occupiedPlus;
		std::vector<bool> occupiedMinus;

		double asymmetry;
		bool addAsymmetry;

		// results that might be needed in the end, after the last step
		Eigen::VectorXd eigenvalsplus;
		Eigen::VectorXd eigenvalsminus;
		Eigen::MatrixXd Cplus;
		Eigen::MatrixXd Cminus;

		UnrestrictedHartreeFock(int iterations = 3000);
		virtual ~UnrestrictedHartreeFock();

		void Init(Systems::Molecule* molecule) override;

		bool DIISStep(int iter, Eigen::MatrixXd& FockMatrixPlus, Eigen::MatrixXd& FockMatrixMinus);
		double Step(int iter) override;

		double CalculateMp2Energy() override;
        // Unique doubles: i<j and a<b for equal spins, all alpha/beta pairs
        // otherwise. Emits the coefficient of a†_a,s a†_b,t a_j,t a_i,s.
        double VisitMp2Amplitudes(
            const std::function<void(bool,bool,int,int,int,int,double)>& emit,
            double minimumGap = 1e-8,
            const std::function<double(bool,int,int,bool,int,int)>& moIntegrals = {});
		double CalculateAtomicCharge(int atom) const override;
		Vector3D<double> GetMoment() const override;

	private:
		DIIS<Eigen::MatrixXd> diisPlus;
		DIIS<Eigen::MatrixXd> diisMinus;

		void CalculateEnergy(const Eigen::VectorXd& eigenvalsplus, const Eigen::VectorXd& eigenvalsminus, const Eigen::MatrixXd& calcDensityMatrixPlus, const Eigen::MatrixXd& calcDensityMatrixMinus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/);
		void InitFockMatrices(int iter, Eigen::MatrixXd& FockMatrixPlus, Eigen::MatrixXd& FockMatrixMinus) const;

	};

}
