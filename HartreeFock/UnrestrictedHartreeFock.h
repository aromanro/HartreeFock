#pragma once
#include "HartreeFockAlgorithm.h"

namespace HartreeFock {

	class UnrestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
	protected:
		double totalEnergy;


		void CalculateEnergy(const Eigen::VectorXd& eigenvalsplus, const Eigen::VectorXd& eigenvalsminus, const Eigen::MatrixXd& calcPplus, const Eigen::MatrixXd& calcPminus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/);
		void InitFockMatrices(int iter, Eigen::MatrixXd& Fplus, Eigen::MatrixXd& Fminus) const;
	public:
		Eigen::MatrixXd Pplus;
		Eigen::MatrixXd Pminus;

		unsigned int nrOccupiedLevelsPlus;
		unsigned int nrOccupiedLevelsMinus;

		UnrestrictedHartreeFock(int iterations = 3000);
		virtual ~UnrestrictedHartreeFock();

		virtual void Init(Systems::Molecule* molecule);

		virtual void Step(int iter);
		virtual double GetTotalEnergy() const;
	};

}