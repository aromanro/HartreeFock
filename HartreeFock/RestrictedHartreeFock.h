#pragma once
#include "HartreeFockAlgorithm.h"

#include <list>

namespace HartreeFock {


	class RestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
		friend class Test;
	protected:
		double totalEnergy;

		std::list<Eigen::MatrixXd> errorMatrices;
		std::list<Eigen::MatrixXd> fockMatrices;

		void CalculateEnergy(const Eigen::VectorXd& eigenvals, const Eigen::MatrixXd& calcDensityMatrix/*, Eigen::MatrixXd& F*/);
		void InitFockMatrix(int iter, Eigen::MatrixXd& FockMatrix) const;
	public:
		Eigen::MatrixXd DensityMatrix;

		unsigned int nrOccupiedLevels;

		// for now in the program it will be filled up with true up to 'nrOccupiedLevels'
		// could be used to compute excited levels, just enlarge it after init and set to true the occupied levels and false the ones that are not occupied
		// adjust the nrOccupied value above then to get the proper homo energy
		std::vector<bool> occupied;

		// results that might be needed in the end, after the last step
		Eigen::VectorXd eigenvals;
		Eigen::MatrixXd Ce; // eigenvectors transformed back into the non-orthogonal AO basis


		RestrictedHartreeFock(int iterations = 3000);
		virtual ~RestrictedHartreeFock();
		
		virtual void Init(Systems::Molecule* molecule);

		virtual double Step(int iter);
		virtual double GetTotalEnergy() const;
	};

}