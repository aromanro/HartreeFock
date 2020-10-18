#pragma once


#include "TensorOrder4.h"
#include "IntegralsRepository.h"

namespace GaussianIntegrals {

	class CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository
	{
	public:
		CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals);
		CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals, IntegralsRepository& repository);


		void Compute(IntegralsRepository& repository);


	protected:
		Tensors::TensorOrder4<double> m_integralsTensor;
	};

}

