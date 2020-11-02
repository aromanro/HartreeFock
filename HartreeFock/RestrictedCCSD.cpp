#include "stdafx.h"
#include "RestrictedCCSD.h"


namespace HartreeFock {

	RestrictedCCSD::RestrictedCCSD(int iterations)
		: RestrictedHartreeFock(iterations), numberOfSpinOrbitals(0), numberOfOccupiedSpinOrbitals(0), m_spinOrbitalBasisIntegrals(nullptr)
	{

	}

	RestrictedCCSD::~RestrictedCCSD()
	{
		delete m_spinOrbitalBasisIntegrals;
	}

	void RestrictedCCSD::Init(Systems::Molecule* molecule)
	{
		RestrictedHartreeFock::Init(molecule);

		numberOfSpinOrbitals = 2 * numberOfOrbitals;
		numberOfOccupiedSpinOrbitals = 2 * nrOccupiedLevels;

		if (m_spinOrbitalBasisIntegrals) delete m_spinOrbitalBasisIntegrals;
		m_spinOrbitalBasisIntegrals = new GaussianIntegrals::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository(numberOfOrbitals);		
	}

}
