#pragma once

#include "MathUtils.h"
#include "GaussianOverlap.h"
#include "GaussianKinetic.h"
#include "GaussianNuclear.h"
#include "GaussianTwoElectrons.h"
#include "GaussianMoment.h"
#include "BoysFunctions.h"

#include <map>
#include <unordered_map>
#include <tuple>
#include <valarray>

#include "Molecule.h"
#include "ContractedGaussianOrbital.h"

namespace GaussianIntegrals {

	// this one is here just to be able to use an unordered_map with a tuple of integers as a key
	template <typename... Tp> class TupleHash {
	public:
		size_t operator()(const std::tuple<Tp...>& t) const
		{
			size_t res = 1;
			std::apply([&res](auto&& ... args) {
				auto compute = [&res](const auto& x) {
					res = 31 * res + x;
				};

				(compute(args), ...);
				}, t);

			return res;
		}
	};

	typedef std::tuple<unsigned int, unsigned int, unsigned int> ThreeOrbitalIndicesTuple;
	typedef TupleHash<unsigned int, unsigned int, unsigned int> ThreeOrbitalIndicesTupleHash;
	typedef std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> FourOrbitalIndicesTuple;
	typedef TupleHash<unsigned int, unsigned int, unsigned int, unsigned int> FourOrbitalIndicesTupleHash;


	class IntegralsRepository
	{
	public:
		IntegralsRepository(Systems::Molecule *molecule = nullptr);

		void Reset(Systems::Molecule* molecule = nullptr);

		const BoysFunctions& getBoysFunctions(unsigned int L, double T);


		double getOverlap(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2, bool extendForKinetic = true);
		double getOverlap(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2, bool extendForKinetic = true);

		double getMoment(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2, bool momentX, bool momentY, bool momentZ);
		double getMoment(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2, bool momentX, bool momentY, bool momentZ);

		double getMomentX(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2)
		{
			return getMoment(gaussian1, gaussian2, true, false, false);
		}

		double getMomentY(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2)
		{
			return getMoment(gaussian1, gaussian2, false, true, false);
		}

		double getMomentZ(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2)
		{
			return getMoment(gaussian1, gaussian2, false, false, true);
		}

		double getMomentX(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2)
		{
			return getMoment(orbital1, orbital2, true, false, false);
		}

		double getMomentY(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2)
		{
			return getMoment(orbital1, orbital2, false, true, false);
		}

		double getMomentZ(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2)
		{
			return getMoment(orbital1, orbital2, false, false, true);
		}


		double getKinetic(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2);
		double getKinetic(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2);


		double getNuclear(const Systems::Atom& atom, const Orbitals::ContractedGaussianOrbital* orbital1, const Orbitals::ContractedGaussianOrbital* orbital2);


		double getElectronElectron(const Orbitals::ContractedGaussianOrbital* orbital1, const Orbitals::ContractedGaussianOrbital* orbital2, const Orbitals::ContractedGaussianOrbital* orbital3, const Orbitals::ContractedGaussianOrbital* orbital4);
		const GaussianTwoElectrons& getElectronElectronVerticalAndTransfer(const Orbitals::GaussianOrbital* orbital1, const Orbitals::GaussianOrbital* orbital2, const Orbitals::GaussianOrbital* orbital3, const Orbitals::GaussianOrbital* orbital4, bool& swapped);


		Systems::Molecule* getMolecule() const { return m_Molecule; }


		void ClearMatricesMaps() 
		{
			//overlapIntegralsMap.clear();
			momentIntegralsMap.clear();
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

		void CalculateElectronElectronIntegrals();

		inline double getElectronElectron(int orbital1, int orbital2, int orbital3, int orbital4) const
		{
			return electronElectronIntegrals[GetElectronElectronIndex(orbital1, orbital2, orbital3, orbital4)];
		}

		void SetUseLotsOfMemory(bool useLots) { useLotsOfMemory = useLots; }

	private:
		const GaussianNuclear& getNuclearVertical(const Systems::Atom& atom, const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2);

		inline void CalculateElectronElectronIntegrals23(int i, const Orbitals::ContractedGaussianOrbital& orb1);
		inline void CalculateElectronElectronIntegrals4(int i, int j, int k, long long int ij, const Orbitals::ContractedGaussianOrbital& orb1, const Orbitals::ContractedGaussianOrbital& orb2, const Orbitals::ContractedGaussianOrbital& orb3);

		inline static long long int GetTwoIndex(long long int i, long long int j)
		{
			return i < j ? j * (j + 1ULL) / 2 + i : i * (i + 1ULL) / 2 + j;
		}

		inline static long long int GetElectronElectronIndex(unsigned int ind1, unsigned int ind2, unsigned int ind3, unsigned int ind4)
		{
			long long int ind12 = GetTwoIndex(ind1, ind2);
			long long int ind34 = GetTwoIndex(ind3, ind4);

			return GetTwoIndex(ind12, ind34);
		}

		template<class Orb> static void SwapOrbitals(Orb** orb1, Orb** orb2, Orb** orb3, Orb** orb4);

		Systems::Molecule* m_Molecule;

		std::map<double, BoysFunctions> boysFunctions;


		// the momentIntegralsMap replaces this, as it also computes overlap
		//std::map < std::tuple<unsigned int, unsigned int, double, double>, GaussianOverlap> overlapIntegralsMap;
		std::map<std::tuple<unsigned int, unsigned int, double, double>, GaussianMoment> momentIntegralsMap; // also contains the overlap, might replace the overlap map as well

		std::map<std::tuple<unsigned int, unsigned int, double, double>, GaussianKinetic> kineticIntegralsMap;

		std::map<std::tuple<unsigned int, unsigned int, unsigned int, double, double>, GaussianNuclear> nuclearVerticalIntegralsMap;
		std::unordered_map<ThreeOrbitalIndicesTuple, GaussianNuclear, ThreeOrbitalIndicesTupleHash> nuclearIntegralsContractedMap;

		std::map<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, double, double, double, double>, GaussianTwoElectrons> electronElectronIntegralsVerticalAndTransferMap;
		std::unordered_map<FourOrbitalIndicesTuple, GaussianTwoElectrons, FourOrbitalIndicesTupleHash> electronElectronIntegralsContractedMap;
		std::valarray<double> electronElectronIntegrals;

		bool useLotsOfMemory = true;
	};


	// this class transforms the integrals that are needed for MP2 from AO basis into MO basis
	// so the used coefficients are those in MO basis, not the transformed ones into the AO basis which are used to obtain the density matrix
	class MolecularOrbitalsIntegralsRepository
	{
	public:
		MolecularOrbitalsIntegralsRepository(const IntegralsRepository& repository) 
			: m_repo(repository)
		{
		}

		inline double getElectronElectron(unsigned int orbital1, unsigned int orbital2, unsigned int orbital3, unsigned int orbital4, const Eigen::MatrixXd& C)
		{
			const FourOrbitalIndicesTuple indTuple = std::make_tuple(orbital1, orbital2, orbital3, orbital4);
			if (m_fourthLevelIntegrals.find(indTuple) != m_fourthLevelIntegrals.end())
				return m_fourthLevelIntegrals.at(indTuple);

			// don't have it yet, compute it
			double result = 0;

			// The very slow method, gets the same results as the faster one
			/*
			for (int i = 0; i < C.cols(); ++i)
			{
				double res1 = 0;
				for (int j = 0; j < C.cols(); ++j)
				{
					double res2 = 0;
					for (int k = 0; k < C.cols(); ++k)
					{
						double res3 = 0;
						for (int l = 0; l < C.cols(); ++l)
							res3 += C(l, orbital4) * m_repo.getElectronElectron(i, j, k, l);
						res2 += C(k, orbital3) * res3;
					}
					res1 += C(j, orbital2) * res2;
				}
				result += C(i, orbital1) * res1;
			}
			*/
			
			for (unsigned int mu = 0; mu < C.cols(); ++mu)
				result += C(mu, orbital1) * getElectronElectronThirdLevel(mu, orbital2, orbital3, orbital4, C);

			m_fourthLevelIntegrals[indTuple] = result;

			return result;
		}

	private:
		inline double getElectronElectronFirstLevel(unsigned int orbital1, unsigned int orbital2, unsigned int orbital3, unsigned int orbital4, const Eigen::MatrixXd& C)
		{
			const FourOrbitalIndicesTuple indTuple = std::make_tuple(orbital1, orbital2, orbital3, orbital4);
			if (m_firstLevelIntegrals.find(indTuple) != m_firstLevelIntegrals.end())
				return m_firstLevelIntegrals.at(indTuple);

			// don't have it yet, compute it
			double result = 0;
			for (unsigned int s = 0; s < C.cols(); ++s)
				result += C(s, orbital4) * m_repo.getElectronElectron(orbital1, orbital2, orbital3, s);
			
			m_firstLevelIntegrals[indTuple] = result;
			
			return result;
		}

		inline double getElectronElectronSecondLevel(unsigned int orbital1, unsigned int orbital2, unsigned int orbital3, unsigned int orbital4, const Eigen::MatrixXd& C)
		{
			const FourOrbitalIndicesTuple indTuple = std::make_tuple(orbital1, orbital2, orbital3, orbital4);
			if (m_secondLevelIntegrals.find(indTuple) != m_secondLevelIntegrals.end())
				return m_secondLevelIntegrals.at(indTuple);

			// don't have it yet, compute it
			double result = 0;

			for (unsigned int l = 0; l < C.cols(); ++l)
				result += C(l, orbital3) * getElectronElectronFirstLevel(orbital1, orbital2, l, orbital4, C);
			
			m_secondLevelIntegrals[indTuple] = result;

			return result;
		}

		inline double getElectronElectronThirdLevel(unsigned int orbital1, unsigned int orbital2, unsigned int orbital3, unsigned int orbital4, const Eigen::MatrixXd& C)
		{
			const FourOrbitalIndicesTuple indTuple = std::make_tuple(orbital1, orbital2, orbital3, orbital4);
			if (m_thirdLevelIntegrals.find(indTuple) != m_thirdLevelIntegrals.end())
				return m_thirdLevelIntegrals.at(indTuple);

			// don't have it yet, compute it
			double result = 0;
			for (unsigned int nu = 0; nu < C.cols(); ++nu)
				result += C(nu, orbital2) * getElectronElectronSecondLevel(orbital1, nu, orbital3, orbital4, C);

			m_thirdLevelIntegrals[indTuple] = result;

			return result;
		}

		const IntegralsRepository& m_repo;
		std::unordered_map<FourOrbitalIndicesTuple, double, FourOrbitalIndicesTupleHash> m_firstLevelIntegrals;
		std::unordered_map<FourOrbitalIndicesTuple, double, FourOrbitalIndicesTupleHash> m_secondLevelIntegrals;
		std::unordered_map<FourOrbitalIndicesTuple, double, FourOrbitalIndicesTupleHash> m_thirdLevelIntegrals;
		std::unordered_map<FourOrbitalIndicesTuple, double, FourOrbitalIndicesTupleHash> m_fourthLevelIntegrals;
	};

}

