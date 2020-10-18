#include "stdafx.h"

#include <assert.h>

#include "CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository.h"


namespace GaussianIntegrals {

	CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals)
		: m_integralsTensor(2 * numOrbitals, 2 * numOrbitals, 2 * numOrbitals, 2 * numOrbitals)
	{
	}
	
	CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals, IntegralsRepository& repository)
		: CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(numOrbitals)
	{
		Compute(repository);
	}


	void CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository::Compute(IntegralsRepository& repository)
	{
		const unsigned int numberSpinOrbitals = m_integralsTensor.GetDim(0);

		for (unsigned int p = 0; p < numberSpinOrbitals; ++p)
			for (unsigned int q = 0; q < numberSpinOrbitals; q++)
				for (unsigned int r = 0; r < numberSpinOrbitals; r++)
					for (unsigned int s = 0; s < numberSpinOrbitals; s++) 
					{
						const double value1 = repository.getElectronElectron(p, r, q, s) * (p % 2 == r % 2) * (q % 2 == s % 2);
						const double value2 = repository.getElectronElectron(p, s, q, r) * (p % 2 == s % 2) * (q % 2 == r % 2);
						m_integralsTensor(p, q, r, s) = value1 - value2;
					}
	}

}
