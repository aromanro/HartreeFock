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
		const size_t numberSpinOrbitals = m_integralsTensor.GetDim(0);

		for (size_t p = 0; p < numberSpinOrbitals; ++p)
		{
			const unsigned int hp = static_cast<unsigned int>(p / 2);
			for (size_t q = 0; q < numberSpinOrbitals; ++q)
			{
				const unsigned int hq = static_cast<unsigned int>(q / 2);
				for (size_t r = 0; r < numberSpinOrbitals; ++r)
				{
					const unsigned int hr = static_cast<unsigned int>(r / 2);
					for (size_t s = 0; s < numberSpinOrbitals; ++s)
					{
						const unsigned int hs = static_cast<unsigned int>(s / 2);

						m_integralsTensor(p, q, r, s) = ((p % 2 == r % 2 && q % 2 == s % 2) ? repository.getElectronElectron(hp, hr, hq, hs, C) : 0) 
							- ((p % 2 == s % 2 && q % 2 == r % 2) ? repository.getElectronElectron(hp, hs, hq, hr, C) : 0);
					}
				}
			}
		}
	}


}
