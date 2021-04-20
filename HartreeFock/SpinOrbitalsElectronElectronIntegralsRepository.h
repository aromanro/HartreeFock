#pragma once


#include "TensorOrder4.h"
#include "IntegralsRepository.h"

namespace GaussianIntegrals {
	
	class SpinOrbitalsElectronElectronIntegralsRepository
	{
	public:
		SpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals);
		SpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals, IntegralsRepository& repository, const Eigen::MatrixXd& C);


		void Compute(IntegralsRepository& repository, const Eigen::MatrixXd& C);

		const double operator()(unsigned int index1, unsigned int index2, unsigned int index3, unsigned int index4) const {
			return m_integralsTensor(index1, index2, index3, index4);
		}
	protected:
		Tensors::TensorOrder4<double> m_integralsTensor;
	};

}

