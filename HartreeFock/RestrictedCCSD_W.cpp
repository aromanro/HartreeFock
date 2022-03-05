#include "stdafx.h"
#include "RestrictedCCSD.h"


namespace HartreeFock {
	// formula 6
	void RestrictedCCSD::CalculateWmnij()
	{
		Wmnij.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals);

		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int hm = m / 2;
			if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

			int indn = 0;
			for (int n = 0; n < numberOfSpinOrbitals; ++n)
			{
				const int hn = n / 2;
				if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied

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

						double sum1 = 0;
						double sum2 = 0;

						int inde = 0;
						for (int e = 0; e < numberOfSpinOrbitals; ++e)
						{
							const int he = e / 2;
							if (he < occupied.size() && occupied[he]) continue; // only unoccupied

							sum1 += t2(indj, inde) * (*m_spinOrbitalBasisIntegrals)(m, n, i, e) - t2(indi, inde) * (*m_spinOrbitalBasisIntegrals)(m, n, j, e);


							int indf = 0;
							for (int f = 0; f < numberOfSpinOrbitals; ++f)
							{
								const int hf = f / 2;
								if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

								sum2 += tau(indi, indj, inde, indf) * (*m_spinOrbitalBasisIntegrals)(m, n, e, f);

								++indf;
							}

							++inde;
						}

						Wmnij(indm, indn, indi, indj) = (*m_spinOrbitalBasisIntegrals)(m, n, i, j) + sum1 + 0.25 * sum2;

						++indj;
					}
					++indi;
				}
				++indn;
			}
			++indm;
		}
	}

	// formula 7
	void RestrictedCCSD::CalculateWabef()
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;
		Wabef.resize(numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

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

				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int he = e / 2;
					if (he < occupied.size() && occupied[he]) continue; // only unoccupied

					int indf = 0;
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int hf = f / 2;
						if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

						double sum1 = 0;
						double sum2 = 0;

						int indm = 0;
						for (int m = 0; m < numberOfSpinOrbitals; ++m)
						{
							const int hm = m / 2;
							if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

							sum1 += t2(indm, indb) * (*m_spinOrbitalBasisIntegrals)(a, m, e, f) - t2(indm, inda) * (*m_spinOrbitalBasisIntegrals)(b, m, e, f);

							int indn = 0;
							for (int n = 0; n < numberOfSpinOrbitals; ++n)
							{
								const int orbn = n / 2;
								if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied

								sum2 += tau(indm, indn, inda, indb) * (*m_spinOrbitalBasisIntegrals)(m, n, e, f);

								++indn;
							}

							++indm;
						}

						Wabef(inda, indb, inde, indf) = (*m_spinOrbitalBasisIntegrals)(a, b, e, f) - sum1 + 0.25 * sum2;

						++indf;
					}
					++inde;
				}
				++indb;
			}
			++inda;
		}
	}

	// formula 8
	void RestrictedCCSD::CalculateWmbej()
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;
		Wmbej.resize(numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals);

		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int hm = m / 2;
			if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

			int indb = 0;
			for (int b = 0; b < numberOfSpinOrbitals; ++b)
			{
				const int hb = b / 2;
				if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int he = e / 2;
					if (he < occupied.size() && occupied[he]) continue; // only unoccupied

					int indj = 0;
					for (int j = 0; j < numberOfSpinOrbitals; ++j)
					{
						const int hj = j / 2;
						if (hj >= occupied.size() || !occupied[hj]) continue; // only occupied

						double sum1 = 0;
						double sum2 = 0;
						double sum3 = 0;

						// sum 1 and sum 3
						int indf = 0;
						for (int f = 0; f < numberOfSpinOrbitals; ++f)
						{
							const int hf = f / 2;
							if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

							sum1 += t2(indj, indf) * (*m_spinOrbitalBasisIntegrals)(m, b, e, f);

							int indn = 0;
							for (int n = 0; n < numberOfSpinOrbitals; ++n)
							{
								const int hn = n / 2;
								if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied

								sum3 += (0.5 * t4(indj, indn, indf, indb) + t2(indj, indf) * t2(indn, indb)) * (*m_spinOrbitalBasisIntegrals)(m, n, e, f);

								++indn;
							}

							++indf;
						}


						// sum2
						int indn = 0;
						for (int n = 0; n < numberOfSpinOrbitals; ++n)
						{
							const int hn = n / 2;
							if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied

							sum2 += t2(indn, indb) * (*m_spinOrbitalBasisIntegrals)(m, n, e, j);

							++indn;
						}

						Wmbej(indm, indb, inde, indj) = (*m_spinOrbitalBasisIntegrals)(m, b, e, j) + sum1 - sum2 - sum3;

						++indj;
					}
					++inde;
				}
				++indb;
			}
			++indm;
		}
	}

}
