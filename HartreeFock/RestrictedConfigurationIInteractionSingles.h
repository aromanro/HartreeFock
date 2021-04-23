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

	private:
		
		// there is some common functionality with coupled cluster - move it in a common base class or member?

		static inline int delta(int i, int j)
		{
			return i == j ? 1 : 0;
		}

		inline Eigen::MatrixXd getSpinOrbitalFockMatrix()
		{
			Eigen::MatrixXd spinOrbitalFockMatrix(numberOfSpinOrbitals, numberOfSpinOrbitals);

			Eigen::MatrixXd FockMatrixMO = m_HartreeFock->C.transpose() * m_HartreeFock->h * m_HartreeFock->C;

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

		RestrictedHartreeFock* m_HartreeFock;

		int numberOfSpinOrbitals;
		int numberOfOccupiedSpinOrbitals;
		
		GaussianIntegrals::SpinOrbitalsElectronElectronIntegralsRepository* m_spinOrbitalBasisIntegrals;

		Eigen::MatrixXd f;
	};


}


