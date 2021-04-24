#include "stdafx.h"
#include "RestrictedConfigurationIInteractionSingles.h"


namespace HartreeFock {

	RestrictedConfigurationIInteractionSingles::RestrictedConfigurationIInteractionSingles(RestrictedHartreeFock* hf)
		: m_HartreeFock(hf), numberOfSpinOrbitals(0), numberOfOccupiedSpinOrbitals(0), m_molecularEEintegrals(nullptr), m_spinOrbitalBasisIntegrals(nullptr)
	{		
	}

	RestrictedConfigurationIInteractionSingles::~RestrictedConfigurationIInteractionSingles()
	{
		delete m_spinOrbitalBasisIntegrals;
	}


	bool RestrictedConfigurationIInteractionSingles::Init()
	{
		if (!m_HartreeFock) return false;

		numberOfSpinOrbitals = 2 * m_HartreeFock->numberOfOrbitals;
		numberOfOccupiedSpinOrbitals = 2 * m_HartreeFock->nrOccupiedLevels;

		if (m_molecularEEintegrals) delete m_molecularEEintegrals;
		m_molecularEEintegrals = new GaussianIntegrals::MolecularOrbitalsIntegralsRepository(m_HartreeFock->integralsRepository);

		if (m_spinOrbitalBasisIntegrals) delete m_spinOrbitalBasisIntegrals;
		m_spinOrbitalBasisIntegrals = new GaussianIntegrals::SpinOrbitalsElectronElectronIntegralsRepository(m_HartreeFock->numberOfOrbitals);

		m_spinOrbitalBasisIntegrals->Compute(*m_molecularEEintegrals, m_HartreeFock->C); // also computes the molecular e-e integrals

		FockMatrixMO = m_HartreeFock->C.transpose() * m_HartreeFock->h * m_HartreeFock->C;

		f = getSpinOrbitalFockMatrix();

		return true;
	}
}

