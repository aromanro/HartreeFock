#pragma once
#include "HartreeFockAlgorithm.h"

namespace HartreeFock {

	class UnrestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
		friend class Test;
	protected:
		double totalEnergy;

		std::list<Eigen::MatrixXd> errorMatricesPlus;
		std::list<Eigen::MatrixXd> errorMatricesMinus;

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

		UnrestrictedHartreeFock(int iterations = 3000);
		virtual ~UnrestrictedHartreeFock();

		virtual void Init(Systems::Molecule* molecule);

		virtual double Step(int iter);
		virtual double GetTotalEnergy() const;
	};

}