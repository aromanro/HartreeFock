#pragma once

#include "RestrictedHartreeFock.h"
#include "SpinOrbitalsElectronElectronIntegralsRepository.h"

namespace HartreeFock {

	class RestrictedConfigurationIInteractionSingles
	{
	public:
		RestrictedConfigurationIInteractionSingles(RestrictedHartreeFock* hf = nullptr);
		virtual ~RestrictedConfigurationIInteractionSingles();

		bool Init();

		inline Eigen::MatrixXd getSpinOrbitalCISMatrix() const
		{
			const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;
			const int matSize = numberOfOccupiedSpinOrbitals * numberOfUnoccupiedSpinOrbitals;

			Eigen::MatrixXd H;
			H.resize(matSize, matSize);

			// TODO: can be optimized, as the matrix is symmetric
			int indi = 0;
			for (int i = 0; i < numberOfSpinOrbitals; ++i)
			{
				const int hi = i / 2;
				if (hi >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[hi]) continue; // only occupied

				const int indibase = indi * numberOfUnoccupiedSpinOrbitals;

				int indj = 0;
				for (int j = 0; j < numberOfSpinOrbitals; ++j)
				{
					const int hj = j / 2;
					if (hj >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[hj]) continue; // only occupied

					const int indjbase = indj * numberOfUnoccupiedSpinOrbitals;

					int inda = 0;
					for (int a = 0; a < numberOfSpinOrbitals; ++a)
					{
						const int ha = a / 2;
						if (ha < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[ha]) continue; // only unoccupied

						const int ind1 = indibase + inda;

						int indb = 0;
						for (int b = 0; b < numberOfSpinOrbitals; ++b)
						{
							const int hb = b / 2;
							if (hb < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[hb]) continue; // only unoccupied

							const int ind2 = indjbase + indb;
							H(ind1, ind2) = (*m_spinOrbitalBasisIntegrals)(a, j, i, b);

							if (delta(i, j)) H(ind1, ind2) += f(a, b);
							if (delta(a, b)) H(ind1, ind2) -= f(i, j);

							++indb;
						}
						++inda;
					}
					++indj;
				}
				++indi;
			}

			return H;
		}



		inline Eigen::MatrixXd getSpinAdaptedCISSinglet()
		{
			const int numberOfUnoccupiedOrbitals = m_HartreeFock->numberOfOrbitals - m_HartreeFock->nrOccupiedLevels;
			const int matSize = m_HartreeFock->nrOccupiedLevels * numberOfUnoccupiedOrbitals;

			Eigen::MatrixXd H;
			H.resize(matSize, matSize);

			int indi = 0;
			for (int i = 0; i < m_HartreeFock->numberOfOrbitals; ++i)
			{
				if (i >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[i]) continue; // only occupied

				const int ti = 2 * i;
				const int indibase = indi * numberOfUnoccupiedOrbitals;

				int indj = 0;
				for (int j = 0; j < m_HartreeFock->numberOfOrbitals; ++j)
				{
					if (j >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[j]) continue; // only occupied

					const int tj = 2 * j;
					const int indjbase = indj * numberOfUnoccupiedOrbitals;

					int inda = 0;
					for (int a = 0; a < m_HartreeFock->numberOfOrbitals; ++a)
					{
						if (a < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[a]) continue; // only unoccupied

						const int ta = 2 * a;
						const int ind1 = indibase + inda;

						int indb = 0;
						for (int b = 0; b < m_HartreeFock->numberOfOrbitals; ++b)
						{
							if (b < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[b]) continue; // only unoccupied

							const int tb = 2 * b;
							const int ind2 = indjbase + indb;
							
							//H(ind1, ind2) = (*m_spinOrbitalBasisIntegrals)(ta, tj, ti, tb) + (*m_spinOrbitalBasisIntegrals)(ta, tj + 1, ti, tb + 1);

							// as above or using directly the molecular spatial orbitals
							H(ind1, ind2) = 2. * m_molecularEEintegrals->getElectronElectron(a, i, j, b, m_HartreeFock->C) - m_molecularEEintegrals->getElectronElectron(a, b, j, i, m_HartreeFock->C);

							if (delta(i, j)) H(ind1, ind2) += f(ta, tb);
							if (delta(a, b)) H(ind1, ind2) -= f(ti, tj);

							++indb;
						}
						++inda;
					}
					++indj;
				}
				++indi;
			}


			return H;
		}



		inline Eigen::MatrixXd getSpinAdaptedCISTriplet()
		{
			const int numberOfUnoccupiedOrbitals = m_HartreeFock->numberOfOrbitals - m_HartreeFock->nrOccupiedLevels;
			const int matSize = m_HartreeFock->nrOccupiedLevels * numberOfUnoccupiedOrbitals;

			Eigen::MatrixXd H;
			H.resize(matSize, matSize);

			int indi = 0;
			for (int i = 0; i < m_HartreeFock->numberOfOrbitals; ++i)
			{
				if (i >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[i]) continue; // only occupied

				const int ti = 2 * i;
				const int indibase = indi * numberOfUnoccupiedOrbitals;

				int indj = 0;
				for (int j = 0; j < m_HartreeFock->numberOfOrbitals; ++j)
				{
					if (j >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[j]) continue; // only occupied

					const int tj = 2 * j;
					const int indjbase = indj * numberOfUnoccupiedOrbitals;

					int inda = 0;
					for (int a = 0; a < m_HartreeFock->numberOfOrbitals; ++a)
					{
						if (a < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[a]) continue; // only unoccupied

						const int ta = 2 * a;
						const int ind1 = indibase + inda;

						int indb = 0;
						for (int b = 0; b < m_HartreeFock->numberOfOrbitals; ++b)
						{
							if (b < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[b]) continue; // only unoccupied

							const int tb = 2 * b;
							const int ind2 = indjbase + indb;

							//H(ind1, ind2) = (*m_spinOrbitalBasisIntegrals)(ta, tj, ti, tb) - (*m_spinOrbitalBasisIntegrals)(ta, tj + 1, ti, tb + 1);

							// as above or using directly the molecular spatial orbitals
							H(ind1, ind2) = -m_molecularEEintegrals->getElectronElectron(a, b, j, i, m_HartreeFock->C);

							if (delta(i, j)) H(ind1, ind2) += f(ta, tb);
							if (delta(a, b)) H(ind1, ind2) -= f(ti, tj);

							++indb;
						}
						++inda;
					}
					++indj;
				}
				++indi;
			}

			return H;
		}



		inline Eigen::MatrixXd getTDHFMatrix() const
		{
			Eigen::MatrixXd A = getSpinOrbitalCISMatrix();
			Eigen::MatrixXd B = getB();

			Eigen::MatrixXd TDHFm;

			const unsigned int matSize = static_cast<unsigned int>(2 * A.cols());
			TDHFm.resize(matSize, matSize);

			TDHFm.block(0, 0, A.rows(), A.cols()) = A;
			TDHFm.block(A.rows(), A.cols(), A.rows(), A.cols()) = -A;
			TDHFm.block(0, A.cols(), A.rows(), A.cols()) = B;
			TDHFm.block(A.rows(), 0, A.rows(), A.cols()) = -B;

			return TDHFm;
		}



	private:
		
		// there is some common functionality with coupled cluster - move it in a common base class or member?

		static inline int delta(int i, int j)
		{
			return i == j ? 1 : 0;
		}

		inline Eigen::MatrixXd getSpinOrbitalFockMatrix() const
		{
			Eigen::MatrixXd spinOrbitalFockMatrix(numberOfSpinOrbitals, numberOfSpinOrbitals);

			for (int p = 0; p < numberOfSpinOrbitals; ++p)
			{
				const int hp = p / 2;

				for (int q = 0; q < numberOfSpinOrbitals; ++q)
				{
					spinOrbitalFockMatrix(p, q) = (p % 2 == q % 2 ? 1 : 0) * FockMatrixMO(hp, q / 2);
					for (int m = 0; m < numberOfSpinOrbitals; ++m)
					{
						const int hm = m / 2;
						if (hm >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[hm]) continue; // only occupied

						spinOrbitalFockMatrix(p, q) += (*m_spinOrbitalBasisIntegrals)(p, m, q, m);
					}
				}
			}

			return spinOrbitalFockMatrix;
		}


		inline Eigen::MatrixXd getB() const
		{
			const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;
			const int matSize = numberOfOccupiedSpinOrbitals * numberOfUnoccupiedSpinOrbitals;

			Eigen::MatrixXd B;
			B.resize(matSize, matSize);

			// TODO: can be optimized, as the matrix is symmetric
			int indi = 0;
			for (int i = 0; i < numberOfSpinOrbitals; ++i)
			{
				const int hi = i / 2;
				if (hi >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[hi]) continue; // only occupied

				const int indibase = indi * numberOfUnoccupiedSpinOrbitals;

				int indj = 0;
				for (int j = 0; j < numberOfSpinOrbitals; ++j)
				{
					const int hj = j / 2;
					if (hj >= m_HartreeFock->occupied.size() || !m_HartreeFock->occupied[hj]) continue; // only occupied

					const int indjbase = indj * numberOfUnoccupiedSpinOrbitals;

					int inda = 0;
					for (int a = 0; a < numberOfSpinOrbitals; ++a)
					{
						const int ha = a / 2;
						if (ha < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[ha]) continue; // only unoccupied

						const int ind1 = indibase + inda;

						int indb = 0;
						for (int b = 0; b < numberOfSpinOrbitals; ++b)
						{
							const int hb = b / 2;
							if (hb < m_HartreeFock->occupied.size() && m_HartreeFock->occupied[hb]) continue; // only unoccupied

							const int ind2 = indjbase + indb;
							B(ind1, ind2) = (*m_spinOrbitalBasisIntegrals)(a, b, i, j);

							++indb;
						}
						++inda;
					}
					++indj;
				}
				++indi;
			}

			return B;
		}

		RestrictedHartreeFock* m_HartreeFock;

		int numberOfSpinOrbitals;
		int numberOfOccupiedSpinOrbitals;


		GaussianIntegrals::MolecularOrbitalsIntegralsRepository* m_molecularEEintegrals;
		GaussianIntegrals::SpinOrbitalsElectronElectronIntegralsRepository* m_spinOrbitalBasisIntegrals;

		Eigen::MatrixXd FockMatrixMO;
		Eigen::MatrixXd f;
	};


}


