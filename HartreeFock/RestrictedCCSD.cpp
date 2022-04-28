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
				ComputeNewt2Sum1(sum1, indi, inda);

				double sum2 = 0;
				double sum3 = 0;
				double sum4 = 0;
				double sum5 = 0;
				double sum6 = 0;
			
				ComputeNewt2Sum26(sum2, sum3, sum4, sum5, sum6, indi, inda, i, a);

				newt2(indi, inda) = (f(i, a) + sum1 - sum2 + sum3 - sum4 - 0.5 * (sum5 + sum6)) / D(i, a);
			
				++inda;
			}
			++indi;
		}

		return newt2;
	}


	void RestrictedCCSD::ComputeNewt2Sum1(double& sum1, int indi, int inda) const
	{
		// sum1
		int inde = 0;
		for (int e = 0; e < numberOfSpinOrbitals; ++e)
		{
			const int he = e / 2;
			if (he < occupied.size() && occupied[he]) continue; // only unoccupied

			sum1 += t2(indi, inde) * Fae(inda, inde);

			++inde;
		}
	}
	
	void RestrictedCCSD::ComputeNewt2Sum26(double& sum2, double& sum3, double& sum4, double& sum5, double& sum6, int indi, int inda, int i, int a) const
	{
		int indm = 0;
		for (int m = 0; m < numberOfSpinOrbitals; ++m)
		{
			const int hm = m / 2;
			if (hm >= occupied.size() || !occupied[hm]) continue; // only occupied

			// sum2
			sum2 += t2(indm, inda) * Fmi(indm, indi);

			// sum3
			int inde = 0;
			for (int e = 0; e < numberOfSpinOrbitals; ++e)
			{
				const int he = e / 2;
				if (he < occupied.size() && occupied[he]) continue; // only unoccupied

				sum3 += t4(indi, indm, inda, inde) * Fme(indm, inde);

				// sum5
				int indf = 0;
				for (int fi = 0; fi < numberOfSpinOrbitals; ++fi)
				{
					const int hf = fi / 2;
					if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

					sum5 += t4(indi, indm, inde, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, e, fi);

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
			for (int fi = 0; fi < numberOfSpinOrbitals; ++fi)
			{
				const int hf = fi / 2;
				if (hf < occupied.size() && occupied[hf]) continue; // only unoccupied

				sum4 += t2(indm, indf) * (*m_spinOrbitalBasisIntegrals)(m, a, i, fi);

				++indf;
			}
			++indm;
		}
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
							for (int fi = 0; fi < numberOfSpinOrbitals; ++fi)
							{
								const int orbf = fi / 2;
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

}



