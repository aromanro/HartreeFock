#pragma once
#include "HartreeFockAlgorithm.h"
#include "DIIS.h"

#include <list>

namespace HartreeFock {


	class RestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
		friend class Test;
	protected:
		DIIS<Eigen::MatrixXd> diis;

		void CalculateEnergy(const Eigen::VectorXd& eigenvals, const Eigen::MatrixXd& calcDensityMatrix/*, Eigen::MatrixXd& F*/);
		void InitFockMatrix(int iter, Eigen::MatrixXd& FockMatrix) const;
	public:
		Eigen::MatrixXd DensityMatrix;

		Eigen::MatrixXd LastMOFockMatrix;

		unsigned int nrOccupiedLevels;

		// for now in the program it will be filled up with true up to 'nrOccupiedLevels'
		// could be used to compute excited levels, just enlarge it after init and set to true the occupied levels and false the ones that are not occupied
		// adjust the nrOccupied value above then to get the proper homo energy
		std::vector<bool> occupied;

		// results that might be needed in the end, after the last step
		Eigen::VectorXd eigenvals;
		Eigen::MatrixXd C; // eigenvectors in AO basis


		RestrictedHartreeFock(int iterations = 3000);
		virtual ~RestrictedHartreeFock();
		
		virtual void Init(Systems::Molecule* molecule) override;

		bool DIISStep(int iter, Eigen::MatrixXd& FockMatrix);
		virtual double Step(int iter) override;

		virtual double CalculateMp2Energy() override;
		virtual double CalculateAtomicCharge(int atom) const override;
		virtual Vector3D<double> GetMoment() const override;
	};

}