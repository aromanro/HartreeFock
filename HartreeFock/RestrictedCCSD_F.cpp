#include "stdafx.h"
#include "RestrictedCCSD.h"


namespace HartreeFock {
	// formula 3
	void RestrictedCCSD::CalculateFae()
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;
		Fae.resize(numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

		int inda = 0;
		for (int a = 0; a < numberOfSpinOrbitals; ++a)
		{
			const int ha = a / 2;
			if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied

			int inde = 0;
			for (int e = 0; e < numberOfSpinOrbitals; ++e)
			{
				const int he = e / 2;
				if (he < occupied.size() && occupied[he]) continue; // only unoccupied

				double sum1 = 0;
				double sum2 = 0;
				double sum3 = 0;

				CalculateFaeSums(sum1, sum2, sum3, inda, a, e);

				Fae(inda, inde) = oneminusdelta(inda, inde) * f(a, e) - 0.5 * (sum1 + sum3) + sum2;

				++inde;
			}
			++inda;
		}
	}

	void RestrictedCCSD::CalculateFaeSums(double& sum1, double& sum2, double& sum3, int inda, int a, int e)
	{
		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int hm = m / 2;
			if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

			sum1 += f(m, e) * t2(indm, inda);

			int indf = 0;
			for (int fi = 0; fi < numberOfSpinOrbitals; ++fi)
			{
				const int hf = fi / 2;
				if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

				sum2 += t2(indm, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, fi, e);

				int indn = 0;
				for (int n = 0; n < numberOfSpinOrbitals; ++n)
				{
					const int hn = n / 2;
					if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied

					sum3 += taut(indm, indn, inda, indf) * (*m_spinOrbitalBasisIntegrals)(m, n, e, fi);

					++indn;
				}
				++indf;
			}
			++indm;
		}
	}

	// formula 4
	void RestrictedCCSD::CalculateFmi()
	{
		Fmi.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals);

		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int hm = m / 2;
			if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

			int indi = 0;
			for (int i = 0; i < numberOfSpinOrbitals; ++i)
			{
				const int hi = i / 2;
				if (hi >= occupied.size() || !occupied[hi]) continue; // only occupied

				double sum1 = 0;
				double sum2 = 0;
				double sum3 = 0;

				CalculateFmiSums(sum1, sum2, sum3, indi, m, i);

				Fmi(indm, indi) = oneminusdelta(indm, indi) * f(m, i) + 0.5 * (sum1 + sum3) + sum2;

				++indi;
			}
			++indm;
		}
	}

	void RestrictedCCSD::CalculateFmiSums(double& sum1, double& sum2, double& sum3, int indi, int m, int i)
	{
		int inde = 0;
		for (int e = 0; e < numberOfSpinOrbitals; ++e)
		{
			const int he = e / 2;
			if (he < occupied.size() && occupied[he]) continue; // only unoccupied

			sum1 += t2(indi, inde) * f(m, e);

			int indn = 0;
			for (int n = 0; n < numberOfSpinOrbitals; ++n)
			{
				const int hn = n / 2;
				if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied

				sum2 += t2(indn, inde) * (*m_spinOrbitalBasisIntegrals)(m, n, i, e);

				int indf = 0;
				for (int fi = 0; fi < numberOfSpinOrbitals; ++fi)
				{
					const int hf = fi / 2;
					if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

					sum3 += taut(indi, indn, inde, indf) * (*m_spinOrbitalBasisIntegrals)(m, n, e, fi);

					++indf;
				}
				++indn;
			}
			++inde;
		}
	}


	// formula 5
	void RestrictedCCSD::CalculateFme()
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

		Fme.resize(numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int hm = m / 2;
			if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

			int inde = 0;
			for (int e = 0; e < numberOfSpinOrbitals; ++e)
			{
				const int he = e / 2;
				if (he < occupied.size() && occupied[he]) continue; // only unoccupied

				double nfSum = 0;

				int indn = 0;
				for (int n = 0; n < numberOfSpinOrbitals; ++n)
				{
					const int hn = n / 2;
					if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied

					int indf = 0;
					for (int fi = 0; fi < numberOfSpinOrbitals; ++fi)
					{
						const int hf = fi / 2;
						if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

						nfSum += t2(indn, indf) * (*m_spinOrbitalBasisIntegrals)(m, n, e, fi);

						++indf;
					}
					++indn;
				}


				Fme(indm, inde) = f(m, e) + nfSum;

				++inde;
			}
			++indm;
		}
	}

}
