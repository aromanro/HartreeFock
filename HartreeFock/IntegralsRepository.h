#pragma once

#include "MathUtils.h"
#include "GaussianOverlap.h"
#include "GaussianKinetic.h"
#include "GaussianNuclear.h"
#include "GaussianTwoElectrons.h"
#include "BoysFunctions.h"

#include <map>
#include <tuple>
#include <valarray>

#include "Molecule.h"
#include "ContractedGaussianOrbital.h"

namespace GaussianIntegrals {

	class IntegralsRepository
	{
	protected:
		Systems::Molecule* m_Molecule;

		std::map< double, BoysFunctions > boysFunctions;


		std::map < std::tuple<unsigned int, unsigned int, double, double>, GaussianOverlap> overlapIntegralsMap;
		std::map < std::tuple<unsigned int, unsigned int, double, double >, GaussianKinetic> kineticIntegralsMap;

		std::map < std::tuple<unsigned int, unsigned int, unsigned int, double, double >, GaussianNuclear > nuclearVerticalIntegralsMap;
		std::map < std::tuple<unsigned int, unsigned int, unsigned int>, GaussianNuclear> nuclearIntegralsContractedMap;
		
		std::map < std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, double, double, double, double>, GaussianTwoElectrons> electronElectronIntegralsVerticalAndTransferMap;
		std::map < std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>, GaussianTwoElectrons> electronElectronIntegralsContractedMap;
		std::valarray<double> electronElectronIntegrals;

	public:
		bool useLotsOfMemory;

		IntegralsRepository(Systems::Molecule *molecule = nullptr);
		~IntegralsRepository();

		void Reset(Systems::Molecule* molecule = nullptr);

		const BoysFunctions& getBoysFunctions(unsigned int L, double T);


		double getOverlap(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2, bool extendForKinetic = true);
		double getOverlap(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2, bool extendForKinetic = true);

		double getKinetic(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2);
		double getKinetic(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2);


		double getNuclear(const Systems::Atom& atom, const Orbitals::ContractedGaussianOrbital* orbital1, const Orbitals::ContractedGaussianOrbital* orbital2);


		double getElectronElectron(const Orbitals::ContractedGaussianOrbital* orbital1, const Orbitals::ContractedGaussianOrbital* orbital2, const Orbitals::ContractedGaussianOrbital* orbital3, const Orbitals::ContractedGaussianOrbital* orbital4);
		const GaussianTwoElectrons& getElectronElectronVerticalAndTransfer(const Orbitals::GaussianOrbital* orbital1, const Orbitals::GaussianOrbital* orbital2, const Orbitals::GaussianOrbital* orbital3, const Orbitals::GaussianOrbital* orbital4, bool& swapped);



		Systems::Molecule* getMolecule() const { return m_Molecule; }


		void ClearMatricesMaps() 
		{
			overlapIntegralsMap.clear();
			kineticIntegralsMap.clear();

			nuclearVerticalIntegralsMap.clear();
			nuclearIntegralsContractedMap.clear();
		}

		void ClearElectronElectronIntermediaries()
		{
			electronElectronIntegralsVerticalAndTransferMap.clear();
		}

		void ClearElectronElectronMaps()
		{
			ClearElectronElectronIntermediaries();

			electronElectronIntegralsContractedMap.clear();
		}

		void ClearAllMaps()
		{
			ClearMatricesMaps();
			ClearElectronElectronMaps();
			boysFunctions.clear();
		}
	protected:
		const GaussianNuclear& getNuclearVertical(const Systems::Atom& atom, const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2);

		inline void CalculateElectronElectronIntegrals23(int i, const Orbitals::ContractedGaussianOrbital& orb1);
		inline void CalculateElectronElectronIntegrals4(int i, int j, int k, int ij, const Orbitals::ContractedGaussianOrbital& orb1, const Orbitals::ContractedGaussianOrbital& orb2, const Orbitals::ContractedGaussianOrbital& orb3);
	
	
		inline static int GetElectronElectronIndex(int ind1, int ind2, int ind3, int ind4)
		{
			if (ind1 < ind2) std::swap(ind1, ind2);
			if (ind3 < ind4) std::swap(ind3, ind4);

			int ind12 = ind1 * (ind1 + 1) / 2 + ind2;
			int ind34 = ind3 * (ind3 + 1) / 2 + ind4;

			if (ind12 < ind34) std::swap(ind12, ind34);

			return ind12 * (ind12 + 1) / 2 + ind34;
		}

		template<class Orb> static void SwapOrbitals(Orb **orb1, Orb **orb2, Orb **orb3, Orb **orb4);
	public:
		void CalculateElectronElectronIntegrals();

		inline double getElectronElectron(int orbital1, int orbital2, int orbital3, int orbital4) const
		{
			return electronElectronIntegrals[GetElectronElectronIndex(orbital1, orbital2, orbital3, orbital4)];
		}

	};

}

