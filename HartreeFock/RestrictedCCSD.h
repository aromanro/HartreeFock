#pragma once
#include "RestrictedHartreeFock.h"

#include "SpinOrbitalsElectronElectronIntegralsRepository.h"

#include "DIIS.h"

// implementation guided by this tutorial: https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2305

namespace HartreeFock {

	class RestrictedCCSD :
		public RestrictedHartreeFock
	{
	public:
		RestrictedCCSD(int iterations = 3000);
		virtual ~RestrictedCCSD();

		void Init(Systems::Molecule* molecule) override;

		void InitCC()
		{
			m_spinOrbitalBasisIntegrals->Compute(integralsRepository, C);
			InitialGuessClusterAmplitudes();
			f = getSpinOrbitalFockMatrix();
		}

		// just for checking against the MP2 energy
		double MP2EnergyFromt4() const
		{
			double result = 0;

			int indi = 0;
			for (int i = 0; i < numberOfSpinOrbitals; ++i)
			{
				const int orbi = i / 2;
				if (orbi >= occupied.size() || !occupied[orbi]) continue; // only occupied

				int indj = 0;
				for (int j = 0; j < numberOfSpinOrbitals; ++j)
				{
					const int orbj = j / 2;
					if (orbj >= occupied.size() || !occupied[orbj]) continue; // only occupied

					int inda = 0;
					for (int a = 0; a < numberOfSpinOrbitals; ++a)
					{
						const int orba = a / 2;
						if (orba < occupied.size() && occupied[orba]) continue; // only unoccupied

						int indb = 0;
						for (int b = 0; b < numberOfSpinOrbitals; ++b)
						{
							const int orbb = b / 2;
							if (orbb < occupied.size() && occupied[orbb]) continue;  // only unoccupied

							result += (*m_spinOrbitalBasisIntegrals)(i, j, a, b) * t4(indi, indj, inda, indb);

							++indb;
						}
						++inda;
					}
					++indj;
				}
				++indi;
			}

			return 0.25 * result;
		}
		

		// Step #5: Check for Convergence and Iterate

		// current CC correlation energy (CCSD):

		double CorrelationEnergy() const;

		// triples correction
		double TEnergy() const;

		double StepCC(int iter);

		bool DIISStep(int iter, Eigen::MatrixXd& newt2, Eigen::Tensor<double, 4>& newt4);

		inline double GetCCEnergy() const
		{
			return CCEnergy;
		}

	private:
		static inline int delta(int i, int j)
		{
			return i == j ? 1 : 0;
		}

		static inline int oneminusdelta(int i, int j)
		{
			return i == j ? 0 : 1;
		}

		// Step #1: Preparing the Spin-Orbital Basis Integrals

		inline Eigen::MatrixXd getSpinOrbitalFockMatrix() const
		{
			Eigen::MatrixXd spinOrbitalFockMatrix(numberOfSpinOrbitals, numberOfSpinOrbitals);

			Eigen::MatrixXd FockMatrixMO = C.transpose() * h * C;

			for (int p = 0; p < numberOfSpinOrbitals; ++p)
			{
				const int hp = p / 2;

				for (int q = 0; q < numberOfSpinOrbitals; ++q)
				{
					spinOrbitalFockMatrix(p, q) = (p % 2 == q % 2 ? 1 : 0) * FockMatrixMO(hp, q / 2);
					for (int m = 0; m < numberOfSpinOrbitals; ++m)
					{
						const int hm = m / 2;
						if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

						spinOrbitalFockMatrix(p, q) += (*m_spinOrbitalBasisIntegrals)(p, m, q, m);
					}
				}
			}

			return spinOrbitalFockMatrix;
		}


		// Step #2: Build the Initial-Guess Cluster Amplitudes
		void InitialGuessClusterAmplitudes()
		{
			const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

			t2 = Eigen::MatrixXd::Zero(numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

			t4.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

			int indi = 0;
			for (int i = 0; i < numberOfSpinOrbitals; ++i)
			{
				const int hi = i / 2;
				if (hi >= occupied.size() || !occupied[hi]) continue; // only occupied

				int indj = 0;
				for (int j = 0; j < numberOfSpinOrbitals; ++j)
				{
					const int hj = j / 2;
					if (hj >= occupied.size() || !occupied[hj]) continue; // only occupied

					int inda = 0;
					for (int a = 0; a < numberOfSpinOrbitals; ++a)
					{
						const int ha = a / 2;
						if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied

						int indb = 0;
						for (int b = 0; b < numberOfSpinOrbitals; ++b)
						{
							const int hb = b / 2;
							if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

							t4(indi, indj, inda, indb) = (*m_spinOrbitalBasisIntegrals)(i, j, a, b) / (eigenvals(hi) + eigenvals(hj) - eigenvals(ha) - eigenvals(hb));

							++indb;
						}
						++inda;
					}
					++indj;
				}
				++indi;
			}
		}


		// Step #3: Calculate the CC Intermediates

		// for the formulae, see J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, J. Chem. Phys. volume 94, pp. 4334-4345 (1991)
		// "A direct product decomposition approach for symmetry exploitation in many-body methods. I. Energy calculations"

		// formulae 3 - 13

		// formulae 9 and 10 first, taus are needed in the following calculations
		void CalculateTaus();

		// 3 - 5
		void CalculateFae();
		void CalculateFaeSums(double& sum1, double& sum2, double& sum3, int inda, int a, int e);

		void CalculateFmi();
		void CalculateFmiSums(double& sum1, double& sum2, double& sum3, int indi, int m, int i);

		void CalculateFme();

		// 6 - 8
		void CalculateWmnij();
		void CalculateWmnijInner(int indm, int indn, int m, int n);
		void CalculateWmnijSums(double& sum1, double& sum2, int indi, int indj, int m, int n, int i, int j);
		void CalculateWabef();
		void CalculateWabefInner(int inda, int indb, int a, int b);
		void CalculateWabefSums(double& sum1, double& sum2, int inda, int indb, int inde, int indf, int a, int b, int e, int fi);
		void CalculateWmbej();
		void CalculateWmbejInner(int indb, int indm, int b, int m);
		void CalculateWmbejSums(double& sum1, double& sum2, double& sum3, int indb, int indj, int b, int e, int m, int j);


		void CalculateIntermediates();

		// end of computing intermediates


		// Step #4: Compute the Updated Cluster Amplitudes

		// eqn 1 and 2

		Eigen::MatrixXd ComputeNewt2() const;
		void ComputeNewt2Sum1(double& sum1, int indi, int inda) const;
		void ComputeNewt2Sum26(double& sum2, double& sum3, double& sum4, double& sum5, double& sum6, int indi, int inda, int i, int a) const;
		void ComputeNewt2Sum4(double& sum4, int indm, int m, int i, int a) const;
		void ComputeNewt2Sum56(double& sum5, double& sum6, int inda, int inde, int indi, int indm, int a, int e, int i, int m) const;

		Eigen::Tensor<double, 4> ComputeNewt4() const;
		void ComputeNewt4Sum146(double& sum1, double& sum4, double& sum6, int inda, int indb, int indi, int indj, int a, int b, int i, int j) const;
		void ComputeNewt4Sum2357(double& sum2, double& sum3, double& sum5, double& sum7, int inda, int indb, int indi, int indj, int a, int b, int i, int j) const;

		// compute denominator arrays

		// formulas 12 and 13

		inline double D(int i, int a) const
		{
			return f(i, i) - f(a, a);
		}

		inline double D(int i, int j, int a, int b) const
		{
			return  f(i, i) + f(j, j) - f(a, a) - f(b, b);
		}

		// needed for the triples correction

		inline double D(int i, int j, int k, int a, int b, int c) const
		{
			return  f(i, i) + f(j, j) + f(k, k) - f(a, a) - f(b, b) - f(c, c);
		}

		double TEnergyInner(double& sum, int inda, int indi, int indj, int indk, int a, int i, int j, int k) const;

		int numberOfSpinOrbitals;
		int numberOfOccupiedSpinOrbitals;


		double CCEnergy;

		Eigen::MatrixXd f;


		Eigen::MatrixXd t2;
		//Eigen::MatrixXd nonExtrapolatedt2;

		// they should be better implemented than my own implementation
		Eigen::Tensor<double, 4> t4;
		//Eigen::Tensor<double, 4> nonExtrapolatedt4;

		GaussianIntegrals::SpinOrbitalsElectronElectronIntegralsRepository* m_spinOrbitalBasisIntegrals;

		// intermediates
		Eigen::MatrixXd Fae; // unoccupied, unoccupied
		Eigen::MatrixXd Fmi; // occupied, occupied
		Eigen::MatrixXd Fme; // occupied, unoccupied

		Eigen::Tensor<double, 4> Wmnij; // occupied, occupied, occupied, occupied
		Eigen::Tensor<double, 4> Wabef; // unoccupied, unoccupied, unoccupied, unoccupied
		Eigen::Tensor<double, 4> Wmbej; // occupied, unoccupied, unoccupied, occupied

		// effective doubles (two particle excitation operators)
		Eigen::Tensor<double, 4> tau; // occupied, occupied, unoccupied, unoccupied
		Eigen::Tensor<double, 4> taut; // occupied, occupied, unoccupied, unoccupied

		//DIIS<Eigen::MatrixXd, 8, 3, false> diisT2;
		// needs specialization for tensors
		//DIIS<Eigen::Tensor<double, 4>, 8, 3, false> diisT4;

		DIIS<Eigen::VectorXd, 8, 2, false> diist;
		Eigen::VectorXd errorVector;
		Eigen::VectorXd valueVector;
	};

}

