#include "stdafx.h"
#include "RestrictedCCSD.h"

namespace HartreeFock {

	double RestrictedCCSD::TEnergy() const
	{
		static const double prefactor = 1. / 36.;

		double sum = 0;

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

				int indk = 0;
				for (int k = 0; k < numberOfSpinOrbitals; ++k)
				{
					const int hk = k / 2;
					if (hk >= occupied.size() || !occupied[hk]) continue; // only occupied

					int inda = 0;
					for (int a = 0; a < numberOfSpinOrbitals; ++a)
					{
						const int ha = a / 2;
						if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied

						TEnergyInner(sum, inda, indi, indj, indk, a, i, j, k);

						++inda;
					}
					++indk;
				}
				++indj;
			}
			++indi;
		}

		return prefactor * sum;
	}


	double RestrictedCCSD::TEnergyInner(double& sum, int inda, int indi, int indj, int indk, int a, int i, int j, int k) const
	{
		int indb = 0;
		for (int b = 0; b < numberOfSpinOrbitals; ++b)
		{
			const int hb = b / 2;
			if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

			int indc = 0;
			for (int c = 0; c < numberOfSpinOrbitals; ++c)
			{
				const int hc = c / 2;
				if (hc < occupied.size() && occupied[hc]) continue; // only unoccupied

				const double Dijkabc = D(i, j, k, a, b, c);

				// connected and disconnected triples

				const double td = (t2(indi, inda) * (*m_spinOrbitalBasisIntegrals)(j, k, b, c) - t2(indi, indb) * (*m_spinOrbitalBasisIntegrals)(j, k, a, c) - t2(indi, indc) * (*m_spinOrbitalBasisIntegrals)(j, k, b, a)
					- (t2(indj, inda) * (*m_spinOrbitalBasisIntegrals)(i, k, b, c) - t2(indj, indb) * (*m_spinOrbitalBasisIntegrals)(i, k, a, c) - t2(indj, indc) * (*m_spinOrbitalBasisIntegrals)(i, k, b, a))
					- (t2(indk, inda) * (*m_spinOrbitalBasisIntegrals)(j, i, b, c) - t2(indk, indb) * (*m_spinOrbitalBasisIntegrals)(j, i, a, c) - t2(indk, indc) * (*m_spinOrbitalBasisIntegrals)(j, i, b, a))
					);
				// / Dijkabc;

				double sum1 = 0;

				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int orbe = e / 2;
					if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied

					sum1 += t4(indj, indk, inda, inde) * (*m_spinOrbitalBasisIntegrals)(e, i, b, c) - t4(indj, indk, indb, inde) * (*m_spinOrbitalBasisIntegrals)(e, i, a, c) - t4(indj, indk, indc, inde) * (*m_spinOrbitalBasisIntegrals)(e, i, b, a)
						- (t4(indi, indk, inda, inde) * (*m_spinOrbitalBasisIntegrals)(e, j, b, c) - t4(indi, indk, indb, inde) * (*m_spinOrbitalBasisIntegrals)(e, j, a, c) - t4(indi, indk, indc, inde) * (*m_spinOrbitalBasisIntegrals)(e, j, b, a))
						- (t4(indj, indi, inda, inde) * (*m_spinOrbitalBasisIntegrals)(e, k, b, c) - t4(indj, indi, indb, inde) * (*m_spinOrbitalBasisIntegrals)(e, k, a, c) - t4(indj, indi, indc, inde) * (*m_spinOrbitalBasisIntegrals)(e, k, b, a));

					++inde;
				}

				double sum2 = 0;

				int indm = 0;
				for (int m = 0; m < numberOfSpinOrbitals; ++m)
				{
					const int orbm = m / 2;
					if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

					sum2 += t4(indi, indm, indb, indc) * (*m_spinOrbitalBasisIntegrals)(m, a, j, k) - t4(indi, indm, inda, indc) * (*m_spinOrbitalBasisIntegrals)(m, b, j, k) - t4(indi, indm, indb, inda) * (*m_spinOrbitalBasisIntegrals)(m, c, j, k)
						- (t4(indj, indm, indb, indc) * (*m_spinOrbitalBasisIntegrals)(m, a, i, k) - t4(indj, indm, inda, indc) * (*m_spinOrbitalBasisIntegrals)(m, b, i, k) - t4(indj, indm, indb, inda) * (*m_spinOrbitalBasisIntegrals)(m, c, i, k))
						- (t4(indk, indm, indb, indc) * (*m_spinOrbitalBasisIntegrals)(m, a, j, i) - t4(indk, indm, inda, indc) * (*m_spinOrbitalBasisIntegrals)(m, b, j, i) - t4(indk, indm, indb, inda) * (*m_spinOrbitalBasisIntegrals)(m, c, j, i));

					++indm;
				}

				const double tc = (sum1 - sum2);
				// / Dijkabc;

				//sum += tc * Dijkabc * (tc + td);
				sum += tc * (tc + td) / Dijkabc; // simplified formula, one division instead of two

				++indc;
			}
			++indb;
		}
	}

}