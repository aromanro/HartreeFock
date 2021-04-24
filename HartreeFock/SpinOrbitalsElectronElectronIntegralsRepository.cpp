#include "stdafx.h"

#include <assert.h>

#include "SpinOrbitalsElectronElectronIntegralsRepository.h"


namespace GaussianIntegrals {

	SpinOrbitalsElectronElectronIntegralsRepository::SpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals)
		: m_integralsTensor(2 * numOrbitals, 2 * numOrbitals, 2 * numOrbitals, 2 * numOrbitals)
	{
	}
	
	SpinOrbitalsElectronElectronIntegralsRepository::SpinOrbitalsElectronElectronIntegralsRepository(int numOrbitals, IntegralsRepository& repository, const Eigen::MatrixXd& C)
		: SpinOrbitalsElectronElectronIntegralsRepository(numOrbitals)
	{
		Compute(repository, C);
	}


	void SpinOrbitalsElectronElectronIntegralsRepository::Compute(IntegralsRepository& repository, const Eigen::MatrixXd& C)
	{
		MolecularOrbitalsIntegralsRepository molecularEEintegrals(repository);

		Compute(molecularEEintegrals, C);
	}


	void SpinOrbitalsElectronElectronIntegralsRepository::Compute(MolecularOrbitalsIntegralsRepository& repository, const Eigen::MatrixXd& C)
	{
		const unsigned int numberSpinOrbitals = m_integralsTensor.GetDim(0);

		for (unsigned int p = 0; p < numberSpinOrbitals; ++p)
		{
			const unsigned int hp = p / 2;
			for (unsigned int q = 0; q < numberSpinOrbitals; ++q)
			{
				const unsigned int hq = q / 2;
				for (unsigned int r = 0; r < numberSpinOrbitals; ++r)
				{
					const unsigned int hr = r / 2;
					for (unsigned int s = 0; s < numberSpinOrbitals; ++s)
					{
						const unsigned int hs = s / 2;

						m_integralsTensor(p, q, r, s) = ((p % 2 == r % 2 && q % 2 == s % 2) ? repository.getElectronElectron(hp, hr, hq, hs, C) : 0) 
							- ((p % 2 == s % 2 && q % 2 == r % 2) ? repository.getElectronElectron(hp, hs, hq, hr, C) : 0);
					}
				}
			}
		}
	}


}
