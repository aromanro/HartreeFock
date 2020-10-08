#pragma once
#include "HartreeFockAlgorithm.h"

namespace HartreeFock {

	class UnrestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
		friend class Test;
	protected:
		std::list<Eigen::MatrixXd> errorMatricesPlus;
		std::list<Eigen::MatrixXd> errorMatricesMinus;

		std::list<Eigen::MatrixXd> fockMatricesPlus;
		std::list<Eigen::MatrixXd> fockMatricesMinus;

		void CalculateEnergy(const Eigen::VectorXd& eigenvalsplus, const Eigen::VectorXd& eigenvalsminus, const Eigen::MatrixXd& calcDensityMatrixPlus, const Eigen::MatrixXd& calcDensityMatrixMinus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/);
		void InitFockMatrices(int iter, Eigen::MatrixXd& FockMatrixPlus, Eigen::MatrixXd& FockMatrixMinus) const;
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
		Eigen::MatrixXd Cplusprime;
		Eigen::MatrixXd Cminusprime;

		UnrestrictedHartreeFock(int iterations = 3000);
		virtual ~UnrestrictedHartreeFock();

		virtual void Init(Systems::Molecule* molecule) override;

		bool DIISStep(int iter, Eigen::MatrixXd& FockMatrixPlus, Eigen::MatrixXd& FockMatrixMinus);
		virtual double Step(int iter) override;

		virtual double CalculateMp2Energy() override;
		virtual double CalculateAtomicCharge(int atom) const override;
		virtual Vector3D<double> GetMoment() const override;
	};

}