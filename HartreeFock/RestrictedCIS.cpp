#include "stdafx.h"
#include "RestrictedCIS.h"


namespace HartreeFock {

	RestrictedCIS::RestrictedCIS(RestrictedHartreeFock* hf)
		: m_HartreeFock(hf), numberOfSpinOrbitals(0), numberOfOccupiedSpinOrbitals(0), m_spinOrbitalBasisIntegrals(nullptr)
	{		
	}

	RestrictedCIS::~RestrictedCIS()
	{
		delete m_spinOrbitalBasisIntegrals;
	}


	bool RestrictedCIS::Init()
	{
		if (!m_HartreeFock) return false;

		numberOfSpinOrbitals = 2 * m_HartreeFock->numberOfOrbitals;
		numberOfOccupiedSpinOrbitals = 2 * m_HartreeFock->nrOccupiedLevels;

		if (m_spinOrbitalBasisIntegrals) delete m_spinOrbitalBasisIntegrals;
		m_spinOrbitalBasisIntegrals = new GaussianIntegrals::SpinOrbitalsElectronElectronIntegralsRepository(m_HartreeFock->numberOfOrbitals);

		m_spinOrbitalBasisIntegrals->Compute(m_HartreeFock->integralsRepository, m_HartreeFock->C);
		f = getSpinOrbitalFockMatrix();

		return true;
	}
}

