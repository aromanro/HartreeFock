#pragma once
#include "HartreeFockAlgorithm.h"

namespace HartreeFock {


	class RestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
	protected:
		double totalEnergy;

		void CalculateEnergy(const Eigen::VectorXd& eigenvals, const Eigen::MatrixXd& calcDensityMatrix/*, Eigen::MatrixXd& F*/);
		void InitFockMatrix(int iter, Eigen::MatrixXd& FockMatrix) const;
	public:
		Eigen::MatrixXd DensityMatrix;

		unsigned int nrOccupiedLevels;

		// for now in the program it will be filled up with true up to 'nrOccupiedLevels'
		// could be used to compute excited levels, just enlarge it after init and set to true the occupied levels and false the ones that are not occupied
		// adjust the nrOccupied value above then to get the proper homo energy
		std::vector<bool> occupied;

		RestrictedHartreeFock(int iterations = 3000);
		virtual ~RestrictedHartreeFock();
		
		virtual void Init(Systems::Molecule* molecule);

		virtual double Step(int iter);
		virtual double GetTotalEnergy() const;
	};

}