#include "stdafx.h"

#include <assert.h>

#include "CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository.h"


namespace GaussianIntegrals {

	CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals)
		: m_integralsTensor(2 * numOrbitals, 2 * numOrbitals, 2 * numOrbitals, 2 * numOrbitals)
	{
	}
	
	CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals, IntegralsRepository& repository, const Eigen::MatrixXd& C)
		: CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(numOrbitals)
	{
		Compute(repository, C);
	}


	void CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository::Compute(IntegralsRepository& repository, const Eigen::MatrixXd& C)
	{
		const unsigned int numberSpinOrbitals = m_integralsTensor.GetDim(0);

		MP2MolecularOrbitalsIntegralsRepository molecularEEintegrals(repository);

		for (unsigned int p = 0; p < numberSpinOrbitals; ++p)
			for (unsigned int q = 0; q < numberSpinOrbitals; ++q)
				for (unsigned int r = 0; r < numberSpinOrbitals; ++r)
					for (unsigned int s = 0; s < numberSpinOrbitals; ++s) 
					{
						const unsigned int hp = p / 2;
						const unsigned int hr = r / 2;
						const unsigned int hq = q / 2;
						const unsigned int hs = s / 2;

						const double value1 = molecularEEintegrals.getElectronElectron(hp, hr, hq, hs, C) * (p % 2 == r % 2) * (q % 2 == s % 2);
						const double value2 = molecularEEintegrals.getElectronElectron(hp, hs, hq, hr, C) * (p % 2 == s % 2) * (q % 2 == r % 2);
						m_integralsTensor(p, q, r, s) = value1 - value2;
					}
	}

}
