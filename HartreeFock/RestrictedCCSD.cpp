#include "stdafx.h"
#include "RestrictedCCSD.h"


namespace HartreeFock {

	RestrictedCCSD::RestrictedCCSD(int iterations)
		: RestrictedHartreeFock(iterations), numberOfSpinOrbitals(0), numberOfOccupiedSpinOrbitals(0), m_spinOrbitalBasisIntegrals(nullptr)
	{

	}

	RestrictedCCSD::~RestrictedCCSD()
	{
		delete m_spinOrbitalBasisIntegrals;
	}

	void RestrictedCCSD::Init(Systems::Molecule* molecule)
	{
		RestrictedHartreeFock::Init(molecule);

		numberOfSpinOrbitals = 2 * numberOfOrbitals;
		numberOfOccupiedSpinOrbitals = 2 * nrOccupiedLevels;

		if (m_spinOrbitalBasisIntegrals) delete m_spinOrbitalBasisIntegrals;
		m_spinOrbitalBasisIntegrals = new GaussianIntegrals::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(numberOfOrbitals);		
	}

	void RestrictedCCSD::CalculateTaus()
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

		tau.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);
		taut.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

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

						const double term2 = t2(indi, inda) * t2(indj, indb) - t2(indi, indb) * t2(indj, inda);

						tau(indi, indj, inda, indb) = t4(indi, indj, inda, indb) + 0.5 * term2;
						taut(indi, indj, inda, indb) = t4(indi, indj, inda, indb) + term2;

						++indb;
					}

					++inda;
				}

				++indj;
			}

			++indi;
		}

	}

	// formula 3
	void RestrictedCCSD::CalculateFae()
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;
		Fae.resize(numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

		int inda = 0;
		for (int a = 0; a < numberOfSpinOrbitals; ++a)
		{
			const int orba = a / 2;
			if (orba < occupied.size() && occupied[orba]) continue; // only unoccupied

			int inde = 0;
			for (int e = 0; e < numberOfSpinOrbitals; ++e)
			{
				const int orbe = e / 2;
				if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied
	
				double sum1 = 0;
				double sum2 = 0;
				double sum3 = 0;
			
				int indm = 0;
				for (int m = 0; m < numberOfSpinOrbitals; ++m)
				{
					const int orbm = m / 2;
					if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

					sum1 += f(m, e) * t2(indm, inda);

					int indf = 0;
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int orbf = f / 2;
						if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied
					
						sum2 += t2(indm, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, f, e);

						int indn = 0;
						for (int n = 0; n < numberOfSpinOrbitals; ++n)
						{
							const int orbn = n / 2;
							if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied
						
							sum3 += taut(indm, indn, inda, indf) * (*m_spinOrbitalBasisIntegrals)(m, n, e, f);

							++indn;
						}					
						++indf;
					}
					++indm;
				}
								
				Fae(inda, inde) = oneminusdelta(inda, inde) * f(a, e) - 0.5 * (sum1 + sum3) + sum2;
				
				++inde;
			}

			++inda;
		}
	}

	// formula 4
	void RestrictedCCSD::CalculateFmi()
	{
		Fmi.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals);

		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int orbm = m / 2;
			if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied
		
			int indi = 0;
			for (int i = 0; i < numberOfSpinOrbitals; ++i)
			{
				const int hi = i / 2;
				if (hi >= occupied.size() || !occupied[hi]) continue; // only occupied
			
				double sum1 = 0;
				double sum2 = 0;
				double sum3 = 0;

				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int orbe = e / 2;
					if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied

					sum1 += f(m, e) * t2(indi, inde);

					int indn = 0;
					for (int n = 0; n < numberOfSpinOrbitals; ++n)
					{
						const int orbn = n / 2;
						if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied

						sum2 += t2(indn, inde) * (*m_spinOrbitalBasisIntegrals)(m, n, i, e);

						int indf = 0;
						for (int f = 0; f < numberOfSpinOrbitals; ++f)
						{
							const int orbf = f / 2;
							if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied

							sum3 += taut(indi, indn, inde, indf) * (*m_spinOrbitalBasisIntegrals)(m, n, e, f);

							++indf;
						}
						++indn;
					}

					++inde;
				}


				Fmi(indm, indi) = oneminusdelta(indm, indi) * f(m, i) + 0.5 * (sum1 + sum3) + sum2;

				++indi;
			}

			++indm;
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
			const int orbm = m / 2;
			if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

			int inde = 0;
			for (int e = 0; e < numberOfSpinOrbitals; ++e)
			{
				const int orbe = e / 2;
				if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied

				double nfSum = 0;

				int indn = 0;
				for (int n = 0; n < numberOfSpinOrbitals; ++n)
				{
					const int orbn = n / 2;
					if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied

					int indf = 0;
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int orbf = f / 2;
						if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied
					
						nfSum += t2(indn, indf) * (*m_spinOrbitalBasisIntegrals)(m, n, e, f);

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



	// formula 6
	void RestrictedCCSD::CalculateWmnij()
	{
		Wmnij.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals);

		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int orbm = m / 2;
			if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

			int indn = 0;
			for (int n = 0; n < numberOfSpinOrbitals; ++n)
			{
				const int orbn = n / 2;
				if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied

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
							const int orbe = e / 2;
							if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied

							sum1 += t2(indj, inde) * (*m_spinOrbitalBasisIntegrals)(m, n, i, e) - t2(indi, inde) * (*m_spinOrbitalBasisIntegrals)(m, n, j, e);


							int indf = 0;
							for (int f = 0; f < numberOfSpinOrbitals; ++f)
							{
								const int orbf = f / 2;
								if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied
							
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
			const int orba = a / 2;
			if (orba < occupied.size() && occupied[orba]) continue; // only unoccupied
		
			int indb = 0;
			for (int b = 0; b < numberOfSpinOrbitals; ++b)
			{
				const int hb = b / 2;
				if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied
			
				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int orbe = e / 2;
					if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied

					int indf = 0;
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int orbf = f / 2;
						if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied

						double sum1 = 0;
						double sum2 = 0;

						int indm = 0;
						for (int m = 0; m < numberOfSpinOrbitals; ++m)
						{
							const int orbm = m / 2;
							if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

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
			const int orbm = m / 2;
			if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

			int indb = 0;
			for (int b = 0; b < numberOfSpinOrbitals; ++b)
			{
				const int hb = b / 2;
				if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int orbe = e / 2;
					if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied

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
							const int orbf = f / 2;
							if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied
						
							sum1 += t2(indj, indf) * (*m_spinOrbitalBasisIntegrals)(m, b, e, f);

							int indn = 0;
							for (int n = 0; n < numberOfSpinOrbitals; ++n)
							{
								const int orbn = n / 2;
								if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied
							
								sum3 += (0.5 * t4(indj, indn, indf, indb) + t2(indj, indf) * t2(indn, indb)) * (*m_spinOrbitalBasisIntegrals)(m, n, e, f);
							
								++indn;
							}

							++f;
						}
					

						// sum2
						int indn = 0;
						for (int n = 0; n < numberOfSpinOrbitals; ++n)
						{
							const int orbn = n / 2;
							if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied
						
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


	void RestrictedCCSD::CalculateIntermediates()
	{
		// arrays
		CalculateTaus();

		CalculateFae();
		CalculateFmi();
		CalculateFme();

		// 4 indexes tensors

		CalculateWmnij();
		CalculateWabef();
		CalculateWmbej();
	}


	// formula 1
	Eigen::MatrixXd RestrictedCCSD::ComputeNewt2() const
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

		Eigen::MatrixXd newt2(numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

		int indi = 0;
		for (int i = 0; i < numberOfSpinOrbitals; ++i)
		{
			const int hi = i / 2;
			if (hi >= occupied.size() || !occupied[hi]) continue; // only occupied
		
			int inda = 0;
			for (int a = 0; a < numberOfSpinOrbitals; ++a)
			{
				const int orba = a / 2;
				if (orba < occupied.size() && occupied[orba]) continue; // only unoccupied
			
				double sum1 = 0;
				double sum2 = 0;
				double sum3 = 0;
				double sum4 = 0;
				double sum5 = 0;
				double sum6 = 0;
			
				// TODO: implement it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				// sum1
				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int orbe = e / 2;
					if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied
				
					sum1 += t2(indi, inde) * Fme(inda, inde);

					++inde;
				}

				int indm = 0;
				for (int m = 0; m < numberOfSpinOrbitals; ++m)
				{
					const int orbm = m / 2;
					if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

					// sum2
					sum2 += t2(indm, inda) * Fmi(indm, indi);


					// sum3
					int inde = 0;
					for (int e = 0; e < numberOfSpinOrbitals; ++e)
					{
						const int orbe = e / 2;
						if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied
					
						sum3 += t4(indi, indm, inda, inde) * Fme(indm, inde);

						// sum5
						int indf = 0;
						for (int f = 0; f < numberOfSpinOrbitals; ++f)
						{
							const int orbf = f / 2;
							if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied
												
							sum5 += t4(indi, indm, inde, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, e, f);

							++indf;
						}

						// sum6
						int indn = 0;
						for (int n = 0; n < numberOfSpinOrbitals; ++n)
						{
							const int orbn = n / 2;
							if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied
						
							sum6 += t4(indm, indn, inda, inde) * (*m_spinOrbitalBasisIntegrals)(n, m, e, i);

							++indn;
						}

						++inde;
					}


					// sum4
					int indf = 0;
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int orbf = f / 2;
						if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied
											
						sum4 += t2(indm, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, i, f);

						++indf;
					}
					++indm;
				}

				newt2(indi, inda) = (f(i, a) + sum1 - sum2 + sum3 - sum4 - 0.5 * (sum5 + sum6)) / D(i, a);
			
				++inda;
			}

			++indi;
		}


		return newt2;
	}

	
	// formula 2
	Eigen::Tensor<double, 4> RestrictedCCSD::ComputeNewt4() const
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

		Eigen::Tensor<double, 4> newt4(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

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
					const int orba = a / 2;
					if (orba < occupied.size() && occupied[orba]) continue; // only unoccupied

					int indb = 0;
					for (int b = 0; b < numberOfSpinOrbitals; ++b)
					{
						const int hb = b / 2;
						if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied
					
					
						// TODO: implement it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

						// this is going to be painful :)
					

						//newt4(indi, indj, inda, indb) = /* a lot of terms */ / D(i, j, a, b);
					
						++indb;
					}
				
					++inda;
				}

				++indj;
			}
			++indi;
		}


		return newt4;
	}


}
