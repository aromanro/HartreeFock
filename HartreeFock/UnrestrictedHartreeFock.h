#pragma once
#include "HartreeFockAlgorithm.h"

namespace HartreeFock {

	class UnrestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
	protected:
		double totalEnergy;


		void CalculateEnergy(const Eigen::VectorXd& eigenvalsplus, const Eigen::VectorXd& eigenvalsminus, const Eigen::MatrixXd& calcDensityMatrixPlus, const Eigen::MatrixXd& calcDensityMatrixMinus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/);
		void InitFockMatrices(int iter, Eigen::MatrixXd& FockMatrixPlus, Eigen::MatrixXd& FockMatrixMinus) const;
	public:
		Eigen::MatrixXd DensityMatrixPlus;
		Eigen::MatrixXd DensityMatrixMinus;

		unsigned int nrOccupiedLevelsPlus;
		unsigned int nrOccupiedLevelsMinus;

		double asymmetry;
		bool addAsymmetry;

		UnrestrictedHartreeFock(int iterations = 3000);
		virtual ~UnrestrictedHartreeFock();

		virtual void Init(Systems::Molecule* molecule);

		virtual void Step(int iter);
		virtual double GetTotalEnergy() const;
	};

}