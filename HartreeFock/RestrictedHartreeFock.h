#pragma once
#include "HartreeFockAlgorithm.h"
#include "DIIS.h"

#include <list>
#include <functional>

namespace HartreeFock {


	class RestrictedHartreeFock :
		public HartreeFockAlgorithm
	{
		friend class Test;
	public:
		Eigen::MatrixXd DensityMatrix;

		//Eigen::MatrixXd LastMOFockMatrix;

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
		
		void Init(Systems::Molecule* molecule) override;

		bool DIISStep(int iter, Eigen::MatrixXd& FockMatrix);
		double Step(int iter) override;

        double CalculateMp2Energy() override;
        // Same RHF MP2 contraction as CalculateMp2Energy, optionally exposing
        // t(i,j,a,b)=(ia|jb)/(eps_i+eps_j-eps_a-eps_b) to correlated clients.
        // A supplied integral accessor can reuse an already transformed MO array.
        double VisitMp2Amplitudes(const std::function<void(int,int,int,int,double)>& emit,
            double minimumGap = 1e-8,
            const std::function<double(int,int,int,int)>& moIntegrals = {});
		double CalculateAtomicCharge(int atom) const override;
		Vector3D<double> GetMoment() const override;

	private:
		DIIS<Eigen::MatrixXd> diis;

		void CalculateEnergy(const Eigen::VectorXd& eigenvals, const Eigen::MatrixXd& calcDensityMatrix/*, Eigen::MatrixXd& F*/);
		void InitFockMatrix(int iter, Eigen::MatrixXd& FockMatrix) const;
	};

}