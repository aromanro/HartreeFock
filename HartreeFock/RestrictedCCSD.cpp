#include "stdafx.h"
#include "RestrictedCCSD.h"


namespace HartreeFock {

	RestrictedCCSD::RestrictedCCSD(int iterations)
		: RestrictedHartreeFock(iterations), CCEnergy(std::numeric_limits<double>::infinity()), numberOfSpinOrbitals(0), numberOfOccupiedSpinOrbitals(0), m_spinOrbitalBasisIntegrals(nullptr)
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
		m_spinOrbitalBasisIntegrals = new GaussianIntegrals::SpinOrbitalsElectronElectronIntegralsRepository(numberOfOrbitals);		
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

						const double term1 = t4(indi, indj, inda, indb);
						const double term2 = t2(indi, inda) * t2(indj, indb) - t2(indi, indb) * t2(indj, inda);

						taut(indi, indj, inda, indb) = term1 + 0.5 * term2;
						tau(indi, indj, inda, indb) = term1 + term2;

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
			
				int indm = 0;
				for (int m = 0; m < numberOfSpinOrbitals; ++m)
				{
					const int hm = m / 2;
					if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

					sum1 += f(m, e) * t2(indm, inda);

					int indf = 0;
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int hf = f / 2;
						if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied
					
						sum2 += t2(indm, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, f, e);

						int indn = 0;
						for (int n = 0; n < numberOfSpinOrbitals; ++n)
						{
							const int hn = n / 2;
							if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied
						
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
						for (int f = 0; f < numberOfSpinOrbitals; ++f)
						{
							const int hf = f / 2;
							if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

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
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int hf = f / 2;
						if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied
					
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
				const int ha = a / 2;
				if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied
			
				double sum1 = 0;
				double sum2 = 0;
				double sum3 = 0;
				double sum4 = 0;
				double sum5 = 0;
				double sum6 = 0;
			
				// sum1
				int inde = 0;
				for (int e = 0; e < numberOfSpinOrbitals; ++e)
				{
					const int he = e / 2;
					if (he < occupied.size() && occupied[he]) continue; // only unoccupied
				
					sum1 += t2(indi, inde) * Fae(inda, inde);

					++inde;
				}

				int indm = 0;
				for (int m = 0; m < numberOfSpinOrbitals; ++m)
				{
					const int hm = m / 2;
					if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

					// sum2
					sum2 += t2(indm, inda) * Fmi(indm, indi);

					// sum3
					inde = 0;
					for (int e = 0; e < numberOfSpinOrbitals; ++e)
					{
						const int he = e / 2;
						if (he < occupied.size() && occupied[he]) continue; // only unoccupied
					
						sum3 += t4(indi, indm, inda, inde) * Fme(indm, inde);

						// sum5
						int indf = 0;
						for (int f = 0; f < numberOfSpinOrbitals; ++f)
						{
							const int hf = f / 2;
							if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied
												
							sum5 += t4(indi, indm, inde, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, e, f);

							++indf;
						}

						// sum6
						int indn = 0;
						for (int n = 0; n < numberOfSpinOrbitals; ++n)
						{
							const int hn = n / 2;
							if (hn >= occupied.size() || !occupied[hn]) continue; // only occupied
						
							sum6 += t4(indm, indn, inda, inde) * (*m_spinOrbitalBasisIntegrals)(n, m, e, i);

							++indn;
						}
						++inde;
					}


					// sum4
					int indf = 0;
					for (int f = 0; f < numberOfSpinOrbitals; ++f)
					{
						const int hf = f / 2;
						if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied
											
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
					const int ha = a / 2;
					if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied

					int indb = 0;
					for (int b = 0; b < numberOfSpinOrbitals; ++b)
					{
						const int hb = b / 2;
						if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied
										
						double sum1 = 0;
						double sum2 = 0;
						double sum3 = 0;
						double sum4 = 0;
						double sum5 = 0;
						double sum6 = 0;
						double sum7 = 0;

						// sum1 is a sum over e: can be used also for sum4 (e and f) and sum6 (also over e)

						int inde = 0;
						for (int e = 0; e < numberOfSpinOrbitals; ++e)
						{
							const int orbe = e / 2;
							if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied
						
						
							// sum1 also has a sum over m inside:

							double psum1 = 0;
							double psum2 = 0;

							int indm = 0;
							for (int m = 0; m < numberOfSpinOrbitals; ++m)
							{
								const int orbm = m / 2;
								if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

								const double F = Fme(indm, inde);

								psum1 += t2(indm, indb) * F;
								psum2 += t2(indm, inda) * F;

								++indm;
							}

							sum1 += t4(indi, indj, inda, inde) * (Fae(indb, inde) - 0.5 * psum1) - t4(indi, indj, indb, inde) * (Fae(inda, inde) - 0.5 * psum2);

							// sum4:
						
							int indf = 0;
							for (int f = 0; f < numberOfSpinOrbitals; ++f)
							{
								const int orbf = f / 2;
								if (orbf < occupied.size() && occupied[orbf]) continue; // only unoccupied
							
								sum4 += tau(indi, indj, inde, indf) * Wabef(inda, indb, inde, indf);
							
								++indf;
							}

							// sum6:

							sum6 += t2(indi, inde) * (*m_spinOrbitalBasisIntegrals)(a, b, e, j) - t2(indj, inde) * (*m_spinOrbitalBasisIntegrals)(a, b, e, i);

							++inde;
						}


						

						// sum2 is a sum over m: can be used also for sum3 (over m and n), sum5 (over m and e), sum7 (over m)

						int indm = 0;
						for (int m = 0; m < numberOfSpinOrbitals; ++m)
						{
							const int orbm = m / 2;
							if (orbm >= occupied.size() || !occupied[orbm]) continue; // only occupied

							// sum2 also has a sum over e inside:

							double psum1 = 0;
							double psum2 = 0;

							inde = 0;
							for (int e = 0; e < numberOfSpinOrbitals; ++e)
							{
								const int orbe = e / 2;
								if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied

								const double F = Fme(indm, inde);

								psum1 += t2(indj, inde) * F;
								psum2 += t2(indi, inde) * F;

								++inde;
							}

							sum2 += t4(indi, indm, inda, indb) * (Fmi(indm, indj) + 0.5 * psum1) - t4(indj, indm, inda, indb) * (Fmi(indm, indi) + 0.5 * psum2);

							// sum3:

							int indn = 0;
							for (int n = 0; n < numberOfSpinOrbitals; ++n)
							{
								const int orbn = n / 2;
								if (orbn >= occupied.size() || !occupied[orbn]) continue; // only occupied
							
								sum3 += tau(indm, indn, inda, indb) * Wmnij(indm, indn, indi, indj);

								++indn;
							}

							// sum5:

							// this sum has permutation for both i, j and a, b

							inde = 0;
							for (int e = 0; e < numberOfSpinOrbitals; ++e)
							{
								const int orbe = e / 2;
								if (orbe < occupied.size() && occupied[orbe]) continue; // only unoccupied
							
								sum5 += t4(indi, indm, inda, inde) * Wmbej(indm, indb, inde, indj) - t2(indi, inde) * t2(indm, inda) * (*m_spinOrbitalBasisIntegrals)(m, b, e, j)
									 - (t4(indi, indm, indb, inde) * Wmbej(indm, inda, inde, indj) - t2(indi, inde) * t2(indm, indb) * (*m_spinOrbitalBasisIntegrals)(m, a, e, j))
									 - (t4(indj, indm, inda, inde) * Wmbej(indm, indb, inde, indi) - t2(indj, inde) * t2(indm, inda) * (*m_spinOrbitalBasisIntegrals)(m, b, e, i))
									 + (t4(indj, indm, indb, inde) * Wmbej(indm, inda, inde, indi) - t2(indj, inde) * t2(indm, indb) * (*m_spinOrbitalBasisIntegrals)(m, a, e, i));

								++inde;
							}

							// sum7:

							sum7 += t2(indm, inda) * (*m_spinOrbitalBasisIntegrals)(m, b, i, j) - t2(indm, indb) * (*m_spinOrbitalBasisIntegrals)(m, a, i, j);

							++indm;
						}


						newt4(indi, indj, inda, indb) = ((*m_spinOrbitalBasisIntegrals)(i, j, a, b) + sum1 - sum2 + 0.5 * (sum3 + sum4) + sum5 + sum6 - sum7) / D(i, j, a, b);
					
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



	double RestrictedCCSD::CorrelationEnergy() const
	{
		const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

		double sum1 = 0;
		double sum2 = 0;
		double sum3 = 0;

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
			
				sum1 += f(i, a) * t2(indi, inda);

				int indj = 0;
				for (int j = 0; j < numberOfSpinOrbitals; ++j)
				{
					const int hj = j / 2;
					if (hj >= occupied.size() || !occupied[hj]) continue; // only occupied
					int indb = 0;
					for (int b = 0; b < numberOfSpinOrbitals; ++b)
					{
						const int hb = b / 2;
						if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

						const double integral = (*m_spinOrbitalBasisIntegrals)(i, j, a, b);

						sum2 += integral * t4(indi, indj, inda, indb);
						sum3 += integral * t2(indi, inda) * t2(indj, indb);

						++indb;
					}
					++indj;
				}
				++inda;
			}
			++indi;
		}

		return sum1 + 0.25 * sum2 + 0.5 * sum3;
	}


	bool RestrictedCCSD::DIISStep(int iter, Eigen::MatrixXd& newt2, Eigen::Tensor<double, 4>& newt4)
	{
		bool UsedDIIS = false;

		if (UseDIIS && iter && iter < maxDIISiterations)
		{
			/*
			{
				const Eigen::MatrixXd errorMatrix = newt2 - t2;

				const double errorNorm = errorMatrix.cwiseProduct(errorMatrix).sum();

				const bool diisCanBeDone = diisT2.AddValueAndError(newt2, errorMatrix);

				if (diisCanBeDone && errorNorm > 1E-14 && iter % 2 == 0)
					lastErrorEst = diisT2.Estimate(newt2);
				else
					lastErrorEst = errorNorm;
			}

			const Eigen::Tensor<double, 4> errorTensor = newt4 - t4;

			const Eigen::Tensor<double, 0> errorNorm = (errorTensor * errorTensor).sum();
			const double errorNormVal = errorNorm();

			const bool diisCanBeDone = diisT4.AddValueAndError(newt4, errorTensor);

			if (diisCanBeDone && errorNormVal > 1E-14 && iter % 2 == 0)
			{
				// use DIIS
				lastErrorEst += diisT4.Estimate(newt4);
				UsedDIIS = true;
			}
			else lastErrorEst += errorNormVal;
			*/

			const long long int vectorSize = newt2.size() + newt4.size();

			errorVector.resize(vectorSize);
			valueVector.resize(vectorSize);

			int index = 0;
			for (int i = 0; i < newt2.rows(); ++i)
				for (int j = 0; j < newt2.cols(); ++j)
				{
					valueVector[index] = newt2(i, j);
					errorVector[index] = newt2(i, j) - t2(i, j);
					++index;
				}

			const std::array<Eigen::Index, 1> one_dim{ { newt4.size() } };
			const Eigen::Tensor<double, 1> t4reshaped = t4.reshape(one_dim);
			Eigen::Tensor<double, 1> newt4reshaped = newt4.reshape(one_dim);

			for (int i = 0; i < newt4.size(); ++i)
			{
				valueVector[index] = newt4reshaped[i];
				errorVector[index] = newt4reshaped[i] - t4reshaped[i];
				++index;
			}

			const double errorNorm = errorVector.cwiseProduct(errorVector).sum();
			const bool diisCanBeDone = diist.AddValueAndError(valueVector, errorVector);

			if (diisCanBeDone && errorNorm > 1E-18 /*&& iter % 2 == 0*/)
			{
				// use DIIS
				lastErrorEst = diist.Estimate(valueVector);
				UsedDIIS = true;

				// now copy back results
				index = 0;
				for (int i = 0; i < newt2.rows(); ++i)
					for (int j = 0; j < newt2.cols(); ++j)
					{
						newt2(i, j) = valueVector[index];
						++index;
					}
				
				for (int i = 0; i < newt4.size(); ++i)
				{
					newt4reshaped[i] = valueVector[index];
					++index;
				}

				newt4 = newt4reshaped.reshape(newt4.dimensions());
			}
			else lastErrorEst = 0;
		}
		else lastErrorEst = 0;

		return UsedDIIS;
	}


	double RestrictedCCSD::StepCC(int iter)
	{
		CalculateIntermediates();

		Eigen::MatrixXd newt2 = ComputeNewt2();
		Eigen::Tensor<double, 4> newt4 = ComputeNewt4();

		/*
		Eigen::MatrixXd savenewt2;
		Eigen::Tensor<double, 4> savenewt4;
		if (UseDIIS)
		{
			savenewt2 = newt2;
			savenewt4 = newt4;
		}
		*/

		DIISStep(iter, newt2, newt4);

		// only for t2, for t4 it needs too many computations
		const Eigen::MatrixXd rmst2dif = newt2 - t2;
		const double rmsD = sqrt(rmst2dif.cwiseProduct(rmst2dif).sum());

		t2 = newt2;
		t4 = newt4;
		/*
		if (UseDIIS)
		{
			nonExtrapolatedt2 = savenewt2;
			nonExtrapolatedt4 = savenewt4;
		}
		*/

		CCEnergy = CorrelationEnergy();

		return rmsD;
	}


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


}



