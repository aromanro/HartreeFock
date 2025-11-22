#pragma once

#undef min
#undef max
#include <Eigen\eigen>

#include "GaussianIntegral.h"

#include "TensorOrder3.h"
#include "TensorOrder4.h"

namespace GaussianIntegrals {


	// the main source of information for implementing this is:
	// HSERILib: Gaussian integral evaluation
	// link: http://theory.rutgers.edu/~giese/notes/HSERILib.pdf
	// it discusses electron-electron integrals but also other kind of integrals very shortly

	class IntegralsRepository;

	class GaussianTwoElectrons : public GaussianIntegral
	{
	public:
		Eigen::MatrixXd matrixCalc;
		Tensors::TensorOrder3<double> tensor3Calc;
		Tensors::TensorOrder4<double> tensor4Calc;


		double getValue(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2, const Orbitals::QuantumNumbers::QuantumNumbers& QN3, const Orbitals::QuantumNumbers::QuantumNumbers& QN4);

		void Reset(IntegralsRepository* repository, double alpha1, double alpha2, double alpha3, double alpha4, 
			const Vector3D<double>& center1, const Vector3D<double>& center2, const Vector3D<double>& center3, const Vector3D<double>& center4, 			
			unsigned int maxL1, unsigned int maxL2, unsigned int maxL3, unsigned int maxL4);

		void HorizontalRecursion1(const Vector3D<double>& dif, unsigned int L1, unsigned int L2, unsigned int L3, unsigned int L4);

		void HorizontalRecursion2(const Vector3D<double>& dif, unsigned int L1, unsigned int L2, unsigned int L3, unsigned int L4);

		inline static bool GetPrevAndPrevPrevAndScalarsForVerticalRecursion(const Orbitals::QuantumNumbers::QuantumNumbers& currentQN, const Vector3D<double>& Rpa, const Vector3D<double>& Rwp, Orbitals::QuantumNumbers::QuantumNumbers& prevQN, Orbitals::QuantumNumbers::QuantumNumbers& prevPrevQN, double& RpaScalar, double& RwpScalar, double& N)
		{
			prevPrevQN = prevQN = currentQN;

			const unsigned int maxIndex = currentQN.MaxComponentVal();

			N = 0;

			bool addPrevPrev;

			if (currentQN.l == maxIndex)
			{
				addPrevPrev = DecrementPrevAndPrevPrevAndSetN(prevQN.l, prevPrevQN.l, N);

				RpaScalar = Rpa.X;
				RwpScalar = Rwp.X;
			}
			else if (currentQN.m == maxIndex)
			{
				addPrevPrev = DecrementPrevAndPrevPrevAndSetN(prevQN.m, prevPrevQN.m, N);

				RpaScalar = Rpa.Y;
				RwpScalar = Rwp.Y;
			}
			else
			{
				addPrevPrev = DecrementPrevAndPrevPrevAndSetN(prevQN.n, prevPrevQN.n, N);

				RpaScalar = Rpa.Z;
				RwpScalar = Rwp.Z;
			}

			return addPrevPrev;
		}

		inline static double GetNextAndPrevQNAndScalarDiffForHorizontalRecursion(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2, const Vector3D<double>& dif, Orbitals::QuantumNumbers::QuantumNumbers& nextQN1, Orbitals::QuantumNumbers::QuantumNumbers& prevQN2)
		{
			nextQN1 = QN1;
			prevQN2 = QN2;

			const unsigned int maxIndex = QN2.MaxComponentVal();

			double difScalar;

			if (QN2.l == maxIndex)
			{
				++nextQN1.l;
				--prevQN2.l;

				difScalar = dif.X;
			}
			else if (QN2.m == maxIndex)
			{
				++nextQN1.m;
				--prevQN2.m;

				difScalar = dif.Y;
			}
			else
			{
				++nextQN1.n;
				--prevQN2.n;

				difScalar = dif.Z;
			}

			return difScalar;
		}

		inline static void IncrementQNandDecrementLimitIfNeeded(Orbitals::QuantumNumbers::QuantumNumbers& QN, unsigned int& limit)
		{
			const unsigned int oldL = QN;

			if (++QN != oldL) {
				assert(limit > 0);
				--limit;
			}
		}

	private:
		void VerticalRecursion(double alpha, double alpha12, const Vector3D<double>& Rpa, const Vector3D<double>& Rwp, unsigned int maxL);
		void ElectronTransfer(double alpha12, double alpha34, const Vector3D<double>& delta, unsigned int maxL, unsigned int maxL2);

		inline static bool DecrementPrevAndPrevPrevAndSetN(unsigned int& prev, unsigned int& prevPrev, double& N)
		{
			if (--prev > 0) {
				N = prev;
				prevPrev -= 2;

				return true;
			}
			
			return false;			
		}

		inline static bool IncNextDecPrevSetN(unsigned int& next, unsigned int& prev, double& N)
		{
			++next;
			if (prev > 0) {
				N = prev;
				--prev;
				
				return true;
			}

			return false;
		}

		inline static bool GetPrevAndPrevPrevAndScalarForElectronTransfer(const Orbitals::QuantumNumbers::QuantumNumbers& currentQN, const Vector3D<double>& delta, Orbitals::QuantumNumbers::QuantumNumbers &prevQN, Orbitals::QuantumNumbers::QuantumNumbers &prevPrevQN, double &deltaScalar, double& N, unsigned int& maxIndex)
		{
			prevPrevQN = prevQN = currentQN;

			maxIndex = currentQN.MaxComponentVal();

			N = 0;

			bool addPrevPrev;

			if (currentQN.l == maxIndex)
			{
				addPrevPrev = DecrementPrevAndPrevPrevAndSetN(prevQN.l, prevPrevQN.l, N);

				deltaScalar = delta.X;
			}
			else if (currentQN.m == maxIndex)
			{
				addPrevPrev = DecrementPrevAndPrevPrevAndSetN(prevQN.m, prevPrevQN.m, N);

				deltaScalar = delta.Y;
			}
			else
			{
				addPrevPrev = DecrementPrevAndPrevPrevAndSetN(prevQN.n, prevPrevQN.n, N);

				deltaScalar = delta.Z;
			}

			return addPrevPrev;
		}
	};

}
