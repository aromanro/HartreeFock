#pragma once
#include "HartreeFockAlgorithm.h"

namespace HartreeFock {

	class UnrestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
	protected:
		double totalEnergy;

		Eigen::MatrixXd Pplus;
		Eigen::MatrixXd Pminus;

		void CalculateEnergy(const Eigen::VectorXd& eigenvalsplus, const Eigen::VectorXd& eigenvalsminus, const Eigen::MatrixXd& calcPplus, const Eigen::MatrixXd& calcPminus/*, const Eigen::MatrixXd& Fplus, const Eigen::MatrixXd& Fminus*/);
		void InitFockMatrices(int iter, Eigen::MatrixXd& Fplus, Eigen::MatrixXd& Fminus) const;
	public:
		unsigned int nrLevelsPlus;
		unsigned int nrLevelsMinus;


		UnrestrictedHartreeFock(int iterations = 3000);
		virtual ~UnrestrictedHartreeFock();

		virtual void Init(Systems::Molecule* molecule);

		virtual void Step(int iter);
		virtual double GetTotalEnergy() const;
	};

}