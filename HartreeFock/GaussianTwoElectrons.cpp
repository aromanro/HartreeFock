#include "stdafx.h"
#include "GaussianTwoElectrons.h"

#include "IntegralsRepository.h"

namespace GaussianIntegrals {


	GaussianTwoElectrons::GaussianTwoElectrons()
	{
	}


	GaussianTwoElectrons::~GaussianTwoElectrons()
	{
	}

	double GaussianTwoElectrons::getValue(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2, const Orbitals::QuantumNumbers::QuantumNumbers& QN3, const Orbitals::QuantumNumbers::QuantumNumbers& QN4)
	{
		return tensor4Calc(QN1.GetCanonicalIndex(), QN2.GetCanonicalIndex(), QN3.GetCanonicalIndex(), QN4.GetCanonicalIndex());
	}



	void GaussianTwoElectrons::Reset(IntegralsRepository* repository, double alpha1, double alpha2, double alpha3, double alpha4, 
		const Vector3D<double>& center1, const Vector3D<double>& center2, const Vector3D<double>& center3, const Vector3D<double>& center4, 
		unsigned int maxL1, unsigned int maxL2, unsigned int maxL3, unsigned int maxL4)
	{
		const double alpha12 = alpha1 + alpha2;
		const double alpha34 = alpha3 + alpha4;
		const double alphaProd = alpha12 * alpha34;
		const double alphaSum = alpha12 + alpha34;
		const double alpha = alphaProd / alphaSum;

		const unsigned int maxL12 = maxL1 + maxL2;
		const unsigned int maxL34 = maxL3 + maxL4;


		const unsigned int maxL = maxL12 + maxL34;
		const unsigned int size = maxL + 1;

		const Orbitals::QuantumNumbers::QuantumNumbers maxQN(0, 0, maxL);
		const unsigned int maxIndex = maxQN.GetTotalCanonicalIndex();

		const Vector3D<double> R12 = center1 - center2;
		const Vector3D<double> R34 = center3 - center4;

		const Vector3D<double> Rp = (alpha1 * center1 + alpha2 * center2) / alpha12;
		const Vector3D<double> Rq = (alpha3 * center3 + alpha4 * center4) / alpha34;
		const Vector3D<double> Rpq = Rp - Rq;

		// auxiliary integrals

		const double exponent = -alpha1 * alpha2 / alpha12 * (R12 * R12) - alpha3 * alpha4 / alpha34 * (R34 * R34);
		const double factor = 2. * pow(M_PI, 5. / 2.) / (alphaProd * sqrt(alphaSum)) * exp(exponent);
		const double T = alpha * (Rpq * Rpq);

		matrixCalc = Eigen::MatrixXd::Zero(maxIndex + 1, size);

		const BoysFunctions& boys = repository->getBoysFunctions(maxL, T);

		for (unsigned int i = 0; i < size; ++i)
			matrixCalc(0, i) = factor * boys.functions[i];

	    // at this point the first row of matrixCalc contains the Boys functions from m = 0 up to m = L1 + L2 + L3 + L4

	    // *********************************************************

		VerticalRecursion(alpha, alpha12, Rp - center1, -alpha / alpha12 * Rpq, maxL);

		// after the above call (apart from m != 0 intermediary results), matrixCalc contains on the 0 column the integrals from (s, s | s, s) to (L1 + L2 + L3 + L4, s | s, s)

		// ************************************************************

		const Orbitals::QuantumNumbers::QuantumNumbers maxQN12(0, 0, maxL12);
		const unsigned int maxIndex12 = maxQN12.GetTotalCanonicalIndex();

		const Orbitals::QuantumNumbers::QuantumNumbers maxQN34(0, 0, maxL34);
		const unsigned int maxIndex34 = maxQN34.GetTotalCanonicalIndex();

		const unsigned int indexL1 = Orbitals::QuantumNumbers::QuantumNumbers(maxL1, 0, 0).GetTotalCanonicalIndex();
		const unsigned int indexL3 = Orbitals::QuantumNumbers::QuantumNumbers(maxL3, 0, 0).GetTotalCanonicalIndex();


		// to start the electron transfer, start with the above calculated 0 column, but enlarge the matrixCalc to hold all 'transferred' integrals in columns
		// transfer from (L1 + L2 + L3 + L4, s | s, s) -> (L1 + L2, s | L3 + L4, s)
		// so it needs as many lines as before and maxIndex34 columns to be able to hold all results, including the intermediary ones

		Eigen::MatrixXd electronTransfer = Eigen::MatrixXd::Zero(maxIndex + 1, maxIndex34 + 1);
		electronTransfer.col(0) = matrixCalc.col(0);  // copy the 0 column into the new matrixCalc
		matrixCalc = electronTransfer;

		const Vector3D<double> Delta = -(alpha2 * R12 + alpha4 * R34) / alpha34;

		// ***************************************************************************************************************************

		ElectronTransfer(alpha12, alpha34, Delta, maxL, maxL34);

		// at this point the matrix holds 0 -> L1 + L2 + L3 + L4 rows (the canonical index) and 0 -> L3 + L4 (the canonical index, not this value) columns
		// not all of them are valid values, see the electron transfer calculation for details

		// **************************************************************************************************************************

		// need to hold only L1 -> L1 + L2 range of rows and L3 -> L3 + L4 range of columns
		matrixCalc = matrixCalc.block(indexL1, indexL3, maxIndex12 - indexL1 + 1, maxIndex34 - indexL3 + 1).eval();
	}

	void GaussianTwoElectrons::VerticalRecursion(double alpha, double alpha12, const Vector3D<double>& Rpa, const Vector3D<double>& Rwp, unsigned int maxL)
	{
		double RpaScalar, RwpScalar;
		double N;

		const unsigned int size = maxL + 1;

		for (auto currentQN = Orbitals::QuantumNumbers::QuantumNumbers(1, 0, 0); currentQN < size; ++currentQN) // for each 'column' starting from 1
		{
			Orbitals::QuantumNumbers::QuantumNumbers prevQN = currentQN;
			Orbitals::QuantumNumbers::QuantumNumbers prevPrevQN = prevQN;

			bool addPrevPrev = GetPrevAndPrevPrevAndScalarsForVerticalRecursion(currentQN, Rpa, Rwp, prevQN, prevPrevQN, RpaScalar, RwpScalar, N);

			unsigned int curIndex = currentQN.GetTotalCanonicalIndex();
			unsigned int prevIndex = prevQN.GetTotalCanonicalIndex();

			for (unsigned int m = 0; m < size - currentQN; ++m)
			{
				// ********************************************************************************************************************************
				// The Vertical Recurrence Relation

				matrixCalc(curIndex, m) = RpaScalar * matrixCalc(prevIndex, m) + RwpScalar * matrixCalc(prevIndex, m + 1);

				if (addPrevPrev)
				{
					unsigned int prevPrevIndex = prevPrevQN.GetTotalCanonicalIndex();
					matrixCalc(curIndex, m) += N / (2. * alpha12) * (matrixCalc(prevPrevIndex, m) - alpha / alpha12 * matrixCalc(prevPrevIndex, m + 1));
				}

				// ********************************************************************************************************************************
			}
		}
	}



	void GaussianTwoElectrons::ElectronTransfer(double alpha12, double alpha34, const Vector3D<double>& delta, unsigned int maxL, unsigned int maxL34)
	{
		unsigned int maxIndex;
		double deltaScalar;
		double Nx = 0;
		double Ny = 0; // assignment is just to keep the compiler happy

		for (auto currentQN2 = Orbitals::QuantumNumbers::QuantumNumbers(1, 0, 0); currentQN2 <= maxL34; ++currentQN2) // for each 'column' starting from 1
		{
			unsigned int curL1 = currentQN2;
			unsigned int curIndexQN2 = currentQN2.GetTotalCanonicalIndex();

			auto prevQN2 = currentQN2;
			auto prevPrevQN2 = prevQN2;

			bool addPrevPrev = GetPrevAndPrevPrevAndScalarForElectronTransfer(currentQN2, delta, prevQN2, prevPrevQN2, deltaScalar, Ny, maxIndex);

			unsigned int prevIndexQN2 = prevQN2.GetTotalCanonicalIndex();
			unsigned int prevPrevIndexQN2 = (addPrevPrev ? prevPrevQN2.GetTotalCanonicalIndex() : 0);

			assert(curL1 <= maxL - curL1);
			assert(maxL - curL1 > 0);

			for (auto currentQN1 = Orbitals::QuantumNumbers::QuantumNumbers(curL1 - 1, 0, 0); currentQN1 <= maxL - curL1; ++currentQN1)
			{
				unsigned int curIndexQN1 = currentQN1.GetTotalCanonicalIndex();

				auto nextQN1 = currentQN1;
				auto prevQN1 = currentQN1;

				bool addPrev;

				if (currentQN2.l == maxIndex)
					addPrev = IncNextDecPrevSetN(nextQN1.l, prevQN1.l, Nx);
				else if (currentQN2.m == maxIndex)
					addPrev = IncNextDecPrevSetN(nextQN1.m, prevQN1.m, Nx);
				else
					addPrev = IncNextDecPrevSetN(nextQN1.n, prevQN1.n, Nx);

				// ********************************************************************************************************************************
				// Electron Transfer Relation

				matrixCalc(curIndexQN1, curIndexQN2) = deltaScalar * matrixCalc(curIndexQN1, prevIndexQN2) - alpha12 / alpha34 * matrixCalc(nextQN1.GetTotalCanonicalIndex(), prevIndexQN2);

				if (addPrev)
					matrixCalc(curIndexQN1, curIndexQN2) += Nx / (2. * alpha34) * matrixCalc(prevQN1.GetTotalCanonicalIndex(), prevIndexQN2);

				if (addPrevPrev)
					matrixCalc(curIndexQN1, curIndexQN2) += Ny / (2. * alpha34) * matrixCalc(curIndexQN1, prevPrevIndexQN2);

				// ********************************************************************************************************************************					
			}
		}
	}





	void GaussianTwoElectrons::HorizontalRecursion1(const Vector3D<double>& dif, unsigned int L1, unsigned int L2, unsigned int L3, unsigned int L4)
	{
		// some values needed later
		const unsigned int ind2Limit = Orbitals::QuantumNumbers::QuantumNumbers(0, 0, L2).GetTotalCanonicalIndex() + 1;
		const unsigned int ind3Limit = (unsigned int)matrixCalc.cols();

		const Orbitals::QuantumNumbers::QuantumNumbers QN1Start(L1, 0, 0);
		const unsigned int QN1Base = QN1Start.GetTotalCanonicalIndex();

		unsigned int L12 = L1 + L2;

		// assertions to check the correctness of some bounds
		assert(matrixCalc.rows() == Orbitals::QuantumNumbers::QuantumNumbers(0, 0, L12).GetTotalCanonicalIndex() - Orbitals::QuantumNumbers::QuantumNumbers(L1, 0, 0).GetTotalCanonicalIndex() + 1);
		assert(ind3Limit == Orbitals::QuantumNumbers::QuantumNumbers(0, 0, L3 + L4).GetTotalCanonicalIndex() - Orbitals::QuantumNumbers::QuantumNumbers(L3, 0, 0).GetTotalCanonicalIndex() + 1);

		// copy the results from the vertical recurrence and electron transfer into a work tensor

		// this holds in the beginning in the first index the L1 -> L1 + L2 range, on the second it holds only s, that is, only the 0 index is filled (but space is reserved to be able to hold all L2 values)
		// on the third index it has the L3 -> L3 + L4 range
		// the later range is not touched, the formula is simply applied on all index values and nothing more, see the 'for' for the ind3 index
		// the second is filled with computation results, as a consequence the first range will decrease to L1 only, the '+ L2' is moved into the second index

		// the formula does this: (L1 + L2, s | L3 + L4, s) -> (L1, L2 | L3 + L4, s)

		Tensors::TensorOrder3<double> workTensor((unsigned int)matrixCalc.rows(), ind2Limit, ind3Limit);

		// copy the values from matrixCalc
		for (unsigned int i = 0; i < matrixCalc.rows(); ++i)
			for (unsigned int j = 0; j < matrixCalc.cols(); ++j)
				workTensor(i, 0, j) = matrixCalc(i, j);

		// clean up the matrixCalc since it's not needed anymore
		matrixCalc.resize(0, 0);

		// the real work - the value for 's' is already in there, here 's' is incremented until reaches L2
		for (auto QN2 = Orbitals::QuantumNumbers::QuantumNumbers(1, 0, 0); QN2 <= L2; IncrementQNandDecrementLimitIfNeeded(QN2, L12))    // for each 'column' starting from 1
		{
			for (auto QN1 = QN1Start; QN1 < L12; ++QN1)
			{
				// ***********************************************************************************************************
				auto nextQN1 = QN1;
				auto prevQN2 = QN2;

				const double difScalar = GetNextAndPrevQNAndScalarDiffForHorizontalRecursion(QN1, QN2, dif, nextQN1, prevQN2);

				const unsigned int curIndex1 = QN1.GetTotalCanonicalIndex() - QN1Base;
				const unsigned int nextIndex1 = nextQN1.GetTotalCanonicalIndex() - QN1Base;

				const unsigned int curIndex2 = QN2.GetTotalCanonicalIndex();
				const unsigned int prevIndex2 = prevQN2.GetTotalCanonicalIndex();

				// ***********************************************************************************************************
				// Horizontal Recurrence Relation 1

				for (unsigned int ind3 = 0; ind3 < ind3Limit; ++ind3)
					workTensor(curIndex1, curIndex2, ind3) = workTensor(nextIndex1, prevIndex2, ind3) + difScalar * workTensor(curIndex1, prevIndex2, ind3);

				// ********************************************************************************************************************************				
			}
		}

		Orbitals::QuantumNumbers::QuantumNumbers QN1(L1, 0, 0);
		Orbitals::QuantumNumbers::QuantumNumbers QN2(L2, 0, 0);
		const unsigned int QN2Base = QN2.GetTotalCanonicalIndex();


		// now copy the values from workTensor into tensor3Calc
		// only L2 values are needed, not the whole 0 -> L2 range from the second index
		// also from the first index only the L1 values are needed, the rest of the values up to L1 + L2 are ignored
		tensor3Calc = Tensors::TensorOrder3<double>(QN1.NumOrbitals(), QN2.NumOrbitals(), workTensor.GetDim(2));

		for (unsigned int i = 0; i < tensor3Calc.GetDim(0); ++i)
			for (unsigned int j = 0; j < tensor3Calc.GetDim(1); ++j)
				for (unsigned int k = 0; k < tensor3Calc.GetDim(2); ++k)
					tensor3Calc(i, j, k) = workTensor(i, QN2Base + j, k);
	}


	// this is very similar with the above, it just transforms (L1, L2 | L3 + L4, s) -> (L1, L2 | L3, L4)

	void GaussianTwoElectrons::HorizontalRecursion2(const Vector3D<double>& dif, unsigned int L1, unsigned int L2, unsigned int L3, unsigned int L4)
	{
		// some values that are needed later
		const Orbitals::QuantumNumbers::QuantumNumbers QN1lim(0, 0, L1);
		const Orbitals::QuantumNumbers::QuantumNumbers QN2lim(0, 0, L2);
		const Orbitals::QuantumNumbers::QuantumNumbers QN3lim(0, 0, L3);
		const Orbitals::QuantumNumbers::QuantumNumbers QN4lim(0, 0, L4);

		const unsigned int limit1 = QN1lim.NumOrbitals();
		const unsigned int limit2 = QN2lim.NumOrbitals();
		const unsigned int limit3 = (unsigned int)tensor3Calc.GetDim(2);
		const unsigned int limit4 = QN4lim.GetTotalCanonicalIndex() + 1;

		unsigned int L34 = L3 + L4;

		const Orbitals::QuantumNumbers::QuantumNumbers QN3Start(L3, 0, 0);
		const unsigned int QN3Base = QN3Start.GetTotalCanonicalIndex();

		const Orbitals::QuantumNumbers::QuantumNumbers QN4Start(L4, 0, 0);
		const unsigned int QN4Base = QN4Start.GetTotalCanonicalIndex();

		// assertions to check that some things are correct

		assert(limit1 == tensor3Calc.GetDim(0));
		assert(limit2 == tensor3Calc.GetDim(1));
		assert(limit3 == Orbitals::QuantumNumbers::QuantumNumbers(0, 0, L3 + L4).GetTotalCanonicalIndex() - QN3Base + 1);

		// copy the results from the vertical recurrence and electron transfer and the first horizontal recursion into a work tensor
		// now an order 4 tensor is needed, to hold all values

		Tensors::TensorOrder4<double> workTensor(limit1, limit2, limit3, limit4);

		for (unsigned int i = 0; i < limit1; ++i)
			for (unsigned int j = 0; j < limit2; ++j)
				for (unsigned int k = 0; k < limit3; ++k)
					workTensor(i, j, k, 0) = tensor3Calc(i, j, k); // in the beginning only s values on the last position

		// don't need the previous results anymore 
		tensor3Calc.Clear();


		// here is the work

		for (auto QN4 = Orbitals::QuantumNumbers::QuantumNumbers(1, 0, 0); QN4 <= L4; IncrementQNandDecrementLimitIfNeeded(QN4, L34))  // for each 'column' starting from 1
		{
			for (auto QN3 = QN3Start; QN3 < L34; ++QN3)
			{
				// ***********************************************************************************************************
				auto nextQN3 = QN3;
				auto prevQN4 = QN4;

				const double difScalar = GetNextAndPrevQNAndScalarDiffForHorizontalRecursion(QN3, QN4, dif, nextQN3, prevQN4);

				const unsigned int curIndex3 = QN3.GetTotalCanonicalIndex() - QN3Base;
				const unsigned int nextIndex3 = nextQN3.GetTotalCanonicalIndex() - QN3Base;

				const unsigned int curIndex4 = QN4.GetTotalCanonicalIndex();
				const unsigned int prevIndex4 = prevQN4.GetTotalCanonicalIndex();

				// ***********************************************************************************************************
				// Horizontal Recurrence Relation 2

				for (unsigned int i = 0; i < limit1; ++i)
					for (unsigned int j = 0; j < limit2; ++j)
						workTensor(i, j, curIndex3, curIndex4) = workTensor(i, j, nextIndex3, prevIndex4) + difScalar * workTensor(i, j, curIndex3, prevIndex4);

				// ********************************************************************************************************************************
			}			
		}


		// copy the result from the work tensor
		// the first indices were not touched, the third range decreased from L3 -> L3 + L4 to L3 only
		// the fourth is now restricted to L4 only, although in the workTensor it started from 0 up to L4

		tensor4Calc = Tensors::TensorOrder4<double>(limit1, limit2, QN3lim.NumOrbitals(), QN4lim.NumOrbitals());

		for (unsigned int i = 0; i < tensor4Calc.GetDim(0); ++i)
			for (unsigned int j = 0; j < tensor4Calc.GetDim(1); ++j)
				for (unsigned int k = 0; k < tensor4Calc.GetDim(2); ++k)
					for (unsigned int l = 0; l < tensor4Calc.GetDim(3); ++l)
						tensor4Calc(i, j, k, l) = workTensor(i, j, k, QN4Base + l);
	}


}