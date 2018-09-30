#pragma once
#include "HartreeFockAlgorithm.h"

namespace HartreeFock {


	class RestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
	protected:
		double totalEnergy;

		unsigned int nrOccupiedLevels;

		void CalculateEnergy(const Eigen::VectorXd& eigenvals, const Eigen::MatrixXd& calcDensityMatrix/*, Eigen::MatrixXd& F*/);
		void InitFockMatrix(int iter, Eigen::MatrixXd& FockMatrix) const;
	public:
		Eigen::MatrixXd DensityMatrix;

		RestrictedHartreeFock(int iterations = 3000);
		virtual ~RestrictedHartreeFock();
		
		virtual void Init(Systems::Molecule* molecule);

		virtual void Step(int iter);
		virtual double GetTotalEnergy() const;
	};

}