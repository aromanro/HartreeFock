#include "stdafx.h"
#include "IntegralsRepository.h"
#include "QuantumNumbers.h"

#include "BoysFunctions.h"

#include <psapi.h>

namespace GaussianIntegrals {


	IntegralsRepository::IntegralsRepository(Systems::Molecule *molecule)
		: m_Molecule(molecule), useLotsOfMemory(true)
	{
	}


	IntegralsRepository::~IntegralsRepository()
	{
	}

	void IntegralsRepository::Reset(Systems::Molecule* molecule)
	{
		ClearAllMaps();

		std::valarray<double> emptyV;
		electronElectronIntegrals.swap(emptyV);

		m_Molecule = molecule;
	}


	//************************************************************************************************************************************************************
	// OVERLAP integrals
	//************************************************************************************************************************************************************

	double IntegralsRepository::getOverlap(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2, bool extendForKinetic)
	{
		assert(m_Molecule);
		
		std::tuple<unsigned int, unsigned int, double, double > params(gaussian1.shellID, gaussian2.shellID, gaussian1.alpha, gaussian2.alpha);

		auto it = overlapIntegralsMap.find(params);
		if (overlapIntegralsMap.end() != it) return it->second.getOverlap(gaussian1.angularMomentum, gaussian2.angularMomentum);


		// unfortunately it's not yet calculated
		GaussianOverlap overlap;
		auto result = overlapIntegralsMap.insert(std::make_pair(params, overlap));



		Orbitals::QuantumNumbers::QuantumNumbers maxQN1(0, 0, 0), maxQN2(0, 0, 0);

		// now find out the maximum quantum numbers

		Systems::AtomWithShells *a1 = nullptr;
		Systems::AtomWithShells *a2 = nullptr;

		for (auto &atom : m_Molecule->atoms)
		{
			if (atom.position == gaussian1.center) a1 = &atom;
			if (atom.position == gaussian2.center) a2 = &atom;

			if (a1 && a2) break;
		}

		assert(a1);
		assert(a2);

		// now find the max quantum numbers

		a1->GetMaxQN(gaussian1.alpha, maxQN1);
		a2->GetMaxQN(gaussian2.alpha, maxQN2);

		if (extendForKinetic)
		{
			// calculating the kinetic integrals needs +1 quantum numbers for overlap integrals
			++maxQN1.l;
			++maxQN1.m;
			++maxQN1.n;

			++maxQN2.n;		
			++maxQN2.l;
			++maxQN2.m;
		}

		// calculate the integrals and that's about it

		result.first->second.Reset(gaussian1.alpha, gaussian2.alpha, gaussian1.center, gaussian2.center, maxQN1, maxQN2);

		return result.first->second.getOverlap(gaussian1.angularMomentum, gaussian2.angularMomentum);
	}


	double IntegralsRepository::getOverlap(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2, bool extendForKinetic)
	{
		double res = 0;
		
		for (auto &gaussian1 : orbital1.gaussianOrbitals)
			for (auto &gaussian2 : orbital2.gaussianOrbitals)
				res += gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient * getOverlap(gaussian1, gaussian2, extendForKinetic);

		return res;
	}



	//************************************************************************************************************************************************************
	// KINETIC integrals
	//************************************************************************************************************************************************************


	double IntegralsRepository::getKinetic(const Orbitals::ContractedGaussianOrbital& orbital1, const Orbitals::ContractedGaussianOrbital& orbital2)
	{
		double res = 0;

		for (auto &gaussian1 : orbital1.gaussianOrbitals)
			for (auto &gaussian2 : orbital2.gaussianOrbitals)
				res += gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient * getKinetic(gaussian1, gaussian2);

		return res;
	}

	double IntegralsRepository::getKinetic(const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2)
	{
		assert(m_Molecule);

		std::tuple<unsigned int, unsigned int, double, double > params(gaussian1.shellID, gaussian2.shellID, gaussian1.alpha, gaussian2.alpha);

		auto it = kineticIntegralsMap.find(params);
		if (kineticIntegralsMap.end() != it) return it->second.getKinetic(gaussian1.angularMomentum, gaussian2.angularMomentum);

		auto oit = overlapIntegralsMap.find(params);
		assert(oit != overlapIntegralsMap.end());

		// unfortunately it's not yet calculated
		GaussianKinetic kinetic(&gaussian1, &gaussian2, &oit->second);
		auto result = kineticIntegralsMap.insert(std::make_pair(params, kinetic));

		Orbitals::QuantumNumbers::QuantumNumbers maxQN1(0, 0, 0), maxQN2(0, 0, 0);

		// now find out the maximum quantum numbers

		Systems::AtomWithShells *a1 = nullptr;
		Systems::AtomWithShells *a2 = nullptr;

		for (auto &atom : m_Molecule->atoms)
		{
			if (atom.position == gaussian1.center) a1 = &atom;
			if (atom.position == gaussian2.center) a2 = &atom;

			if (a1 && a2) break;
		}

		assert(a1);
		assert(a2);

		// now find the max quantum numbers

		a1->GetMaxQN(gaussian1.alpha, maxQN1);
		a2->GetMaxQN(gaussian2.alpha, maxQN2);

		// calculate the integrals and that's about it

		result.first->second.Reset(gaussian1.alpha, gaussian2.alpha, maxQN1, maxQN2);

		return result.first->second.getKinetic(gaussian1.angularMomentum, gaussian2.angularMomentum);
	}


	//************************************************************************************************************************************************************
	// NUCLEAR integrals
	//************************************************************************************************************************************************************


	double IntegralsRepository::getNuclear(const Systems::Atom& nucleus, const Orbitals::ContractedGaussianOrbital* orbital1, const Orbitals::ContractedGaussianOrbital* orbital2)
	{
		if (orbital1->angularMomentum < orbital2->angularMomentum ||
			(orbital1->angularMomentum == orbital2->angularMomentum && orbital1->ID > orbital2->ID)) std::swap(orbital1, orbital2);

		
		std::tuple<unsigned int, unsigned int, unsigned int> params(nucleus.ID, orbital1->ID, orbital2->ID);
		auto it = nuclearIntegralsContractedMap.find(params);
		if (nuclearIntegralsContractedMap.end() != it) return it->second.getNuclear(orbital1->angularMomentum, orbital2->angularMomentum);
		
		const auto center1 = orbital1->getCenter();
		const auto center2 = orbital2->getCenter();

		Orbitals::QuantumNumbers::QuantumNumbers maxQN(0, 0, orbital1->angularMomentum + orbital2->angularMomentum);
	
		GaussianNuclear horizNuclear;
		auto result = nuclearIntegralsContractedMap.insert(std::make_pair(params, horizNuclear));
		result.first->second.matrixCalc = Eigen::MatrixXd::Zero(maxQN.GetTotalCanonicalIndex() + 1, 1);

		
		for (auto &gaussian1 : orbital1->gaussianOrbitals)
			for (auto &gaussian2 : orbital2->gaussianOrbitals)
			{
				const GaussianNuclear& nuclear = getNuclearVertical(nucleus, gaussian1, gaussian2);

				double factor = gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient;

				for (int row = 0; row < result.first->second.matrixCalc.rows(); ++row)
					result.first->second.matrixCalc(row, 0) += factor * nuclear.matrixCalc(row, 0);
			}			

		result.first->second.HorizontalRecursion(center1-center2, orbital1->angularMomentum, orbital2->angularMomentum);

		return result.first->second.getNuclear(orbital1->angularMomentum, orbital2->angularMomentum);
	}


	const GaussianNuclear& IntegralsRepository::getNuclearVertical(const Systems::Atom& nucleus, const Orbitals::GaussianOrbital& gaussian1, const Orbitals::GaussianOrbital& gaussian2)
	{
		assert(m_Molecule);

		std::tuple<unsigned int, unsigned int, unsigned int, double, double > params(nucleus.ID, gaussian1.shellID, gaussian2.shellID, gaussian1.alpha, gaussian2.alpha);
		auto it = nuclearVerticalIntegralsMap.find(params);
		if (nuclearVerticalIntegralsMap.end() != it) return it->second;

		// unfortunately it's not yet calculated
		GaussianNuclear nuclear;
		auto result = nuclearVerticalIntegralsMap.insert(std::make_pair(params, nuclear));

		// now find out the maximum quantum numbers

		Systems::AtomWithShells *a1 = nullptr;
		Systems::AtomWithShells *a2 = nullptr;

		for (auto &atom : m_Molecule->atoms)
		{
			if (atom.position == gaussian1.center) a1 = &atom;
			if (atom.position == gaussian2.center) a2 = &atom;

			if (a1 && a2) break;
		}

		assert(a1);
		assert(a2);

		// now find the max quantum numbers

		unsigned int maxL1 = a1->GetMaxAngularMomentum(gaussian1.alpha);
		unsigned int maxL2 = a2->GetMaxAngularMomentum(gaussian2.alpha);

		// calculate the integrals and that's about it

		result.first->second.Reset(this, gaussian1.alpha, gaussian2.alpha, nucleus.position, gaussian1.center, gaussian2.center, maxL1, maxL2, false);

		return result.first->second;
	}


	// **************************************************************************************************************************************************************
	// ELECTRON - ELECTRON from here
	// **************************************************************************************************************************************************************

	
	const BoysFunctions& IntegralsRepository::getBoysFunctions(unsigned int L, double T)
	{
		auto it = boysFunctions.find(T);		
		if (boysFunctions.end() != it) return it->second;

		BoysFunctions boys;
		auto result = boysFunctions.insert(std::make_pair(T, boys));

		unsigned int maxL = L + 1;

		Systems::Molecule *molecule = getMolecule();	
		if (molecule) maxL = 4 * molecule->GetMaxAngularMomentum() + 1;

		result.first->second.GenerateBoysFunctions(maxL, T);

		return result.first->second;
	}	

	template<class Orb> void IntegralsRepository::SwapOrbitals(Orb **orb1, Orb **orb2, Orb **orb3, Orb **orb4)
	{
		assert(orb1);
		assert(orb2);
		assert(orb3);
		assert(orb4);

		assert(*orb1);
		assert(*orb2);
		assert(*orb3);
		assert(*orb4);

		// the calculation algorithm requires a certain order of angular momenta
		// the L1 >= L2, L3 >= L4 and L1 + L2 >= L3 + L4 are needed by the algorithm
		// this also ensures that symmetry is used to speed it up
		
		if ((*orb1)->angularMomentum < (*orb2)->angularMomentum || ((*orb1)->angularMomentum == (*orb2)->angularMomentum && (*orb1)->ID < (*orb2)->ID)) std::swap(*orb1, *orb2);
		if ((*orb3)->angularMomentum < (*orb4)->angularMomentum || ((*orb3)->angularMomentum == (*orb4)->angularMomentum && (*orb3)->ID < (*orb4)->ID)) std::swap(*orb3, *orb4);

		if ((*orb1)->angularMomentum + (*orb2)->angularMomentum < (*orb3)->angularMomentum + (*orb4)->angularMomentum ||
			((*orb1)->angularMomentum + (*orb2)->angularMomentum == (*orb3)->angularMomentum + (*orb4)->angularMomentum &&
			((*orb1)->ID < (*orb3)->ID || ((*orb1)->ID == (*orb3)->ID && (*orb2)->ID < (*orb4)->ID))))
		{
			std::swap(*orb1, *orb3);
			std::swap(*orb2, *orb4);
		}
	}
	
	double IntegralsRepository::getElectronElectron(const Orbitals::ContractedGaussianOrbital* orbital1, const Orbitals::ContractedGaussianOrbital* orbital2, const Orbitals::ContractedGaussianOrbital* orbital3, const Orbitals::ContractedGaussianOrbital* orbital4)
	{
		SwapOrbitals(&orbital1, &orbital2, &orbital3, &orbital4);

		assert(orbital1->angularMomentum >= orbital2->angularMomentum);
		assert(orbital3->angularMomentum >= orbital4->angularMomentum);
		assert(orbital1->angularMomentum + orbital2->angularMomentum >= orbital3->angularMomentum + orbital4->angularMomentum);
		
		std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int> params(orbital1->angularMomentum, orbital2->angularMomentum, orbital3->angularMomentum, orbital4->angularMomentum, orbital1->shellID, orbital2->shellID, orbital3->shellID, orbital4->shellID);
		auto it = electronElectronIntegralsContractedMap.find(params);
		if (electronElectronIntegralsContractedMap.end() != it) return it->second.getValue(orbital1->angularMomentum, orbital2->angularMomentum, orbital3->angularMomentum, orbital4->angularMomentum);

		const unsigned int L1 = orbital1->angularMomentum;
		const unsigned int L2 = orbital2->angularMomentum;
		const unsigned int L3 = orbital3->angularMomentum;
		const unsigned int L4 = orbital4->angularMomentum;

		const unsigned int maxL12 = L1 + L2;
		const unsigned int maxL34 = L3 + L4;

		GaussianTwoElectrons contractedTwoElectrons;
		auto result = electronElectronIntegralsContractedMap.insert(std::make_pair(params, contractedTwoElectrons));

		// set up a zero filled matrix for the range L1 -> L1 + L2 on columns and L3 -> L3 + L4 on rows, the same range as the result from vertical and electron transfer relations
		result.first->second.matrixCalc = Eigen::MatrixXd::Zero(Orbitals::QuantumNumbers::QuantumNumbers(0, 0, maxL12).GetTotalCanonicalIndex() - Orbitals::QuantumNumbers::QuantumNumbers(L1, 0, 0).GetTotalCanonicalIndex() + 1, Orbitals::QuantumNumbers::QuantumNumbers(0, 0, maxL34).GetTotalCanonicalIndex() - Orbitals::QuantumNumbers::QuantumNumbers(L3, 0, 0).GetTotalCanonicalIndex() + 1);

		// now contract the results from the above mentioned two relations, the horizontal relations can be applied on the contracted results

		for (const auto &gaussian1 : orbital1->gaussianOrbitals)		
			for (const auto &gaussian2 : orbital2->gaussianOrbitals)
				for (const auto &gaussian3 : orbital3->gaussianOrbitals)
					for (const auto &gaussian4 : orbital4->gaussianOrbitals)
					{
						const double factor = gaussian1.normalizationFactor * gaussian2.normalizationFactor *  gaussian3.normalizationFactor * gaussian4.normalizationFactor * 
												gaussian1.coefficient * gaussian2.coefficient * gaussian3.coefficient * gaussian4.coefficient;

						bool swapped;
						const GaussianTwoElectrons& electronsVertical = getElectronElectronVerticalAndTransfer(&gaussian1, &gaussian2, &gaussian3, &gaussian4, swapped);

						if (swapped)
						{
							for (int i = 0; i < result.first->second.matrixCalc.rows(); ++i)
								for (int j = 0; j < result.first->second.matrixCalc.cols(); ++j)
									result.first->second.matrixCalc(i, j) += factor * electronsVertical.matrixCalc(j, i);
						}
						else
						{
							for (int i = 0; i < result.first->second.matrixCalc.rows(); ++i)
								for (int j = 0; j < result.first->second.matrixCalc.cols(); ++j)
									result.first->second.matrixCalc(i, j) += factor * electronsVertical.matrixCalc(i, j);
						}
					}

		// result.first->second.matrixCalc now holds the contraction of the results of vertical and electron transfer relations

		// now apply the two horizontal recurrence relations on it

		result.first->second.HorizontalRecursion1(orbital1->center - orbital2->center, orbital1->angularMomentum, orbital2->angularMomentum, orbital3->angularMomentum, orbital4->angularMomentum);
		result.first->second.HorizontalRecursion2(orbital3->center - orbital4->center, orbital1->angularMomentum, orbital2->angularMomentum, orbital3->angularMomentum, orbital4->angularMomentum);

		return result.first->second.getValue(orbital1->angularMomentum, orbital2->angularMomentum, orbital3->angularMomentum, orbital4->angularMomentum);
	}


	const GaussianTwoElectrons& IntegralsRepository::getElectronElectronVerticalAndTransfer(const Orbitals::GaussianOrbital* orbital1, const Orbitals::GaussianOrbital* orbital2, const Orbitals::GaussianOrbital* orbital3, const Orbitals::GaussianOrbital* orbital4, bool& swapped)
	{
		// the contracted gaussians were already swapped to have the angular momentum in order
		assert(orbital1->angularMomentum >= orbital2->angularMomentum);
		assert(orbital3->angularMomentum >= orbital4->angularMomentum);
		assert(orbital1->angularMomentum + orbital2->angularMomentum >= orbital3->angularMomentum + orbital4->angularMomentum);

		if (orbital1->angularMomentum == orbital2->angularMomentum && orbital1->ID < orbital2->ID) std::swap(orbital1, orbital2);
		if (orbital3->angularMomentum == orbital4->angularMomentum && orbital3->ID < orbital4->ID) std::swap(orbital3, orbital4);

		// the 'swapped' flag is needed to transpose the results matrix
		// the range is (L1 -> L1 + L2, s | L3 -> L3 + L4, s)
		// if 'swapped', it needs to be turned - by transposing - into:
		// (L3 -> L3 + L4, s | L1 -> L1 + L2, s)
		// since we're swapping here only if L1 + L2 == L3 + L4,
		// the 'swapping' flag only takes care of the lower bounds
		if (orbital1->angularMomentum + orbital2->angularMomentum == orbital3->angularMomentum + orbital4->angularMomentum &&
		       (orbital1->ID < orbital3->ID || (orbital1->ID == orbital3->ID && orbital2->ID < orbital4->ID)))
		{
			std::swap(orbital1, orbital3);
			std::swap(orbital2, orbital4);
			swapped = true;
		}
		else swapped = false;

		std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, double, double, double, double> params(orbital1->angularMomentum, orbital2->angularMomentum, orbital3->angularMomentum, orbital4->angularMomentum, 
																																											orbital1->shellID, orbital2->shellID, orbital3->shellID, orbital4->shellID, 
																																											orbital1->alpha, orbital2->alpha, orbital3->alpha, orbital4->alpha);
		auto it = electronElectronIntegralsVerticalAndTransferMap.find(params);
		if (electronElectronIntegralsVerticalAndTransferMap.end() != it) return it->second;
		
		// unfortunately it's not yet calculated
		GaussianTwoElectrons electronElectron;
		auto result = electronElectronIntegralsVerticalAndTransferMap.insert(std::make_pair(params, electronElectron));

		result.first->second.Reset(this, orbital1->alpha, orbital2->alpha, orbital3->alpha, orbital4->alpha, 
			                             orbital1->center, orbital2->center, orbital3->center, orbital4->center, 
			                             orbital1->angularMomentum, orbital2->angularMomentum, orbital3->angularMomentum, orbital4->angularMomentum);

		return result.first->second;
	}


	





	inline void IntegralsRepository::CalculateElectronElectronIntegrals4(int i, int j, int k, int ij, const Orbitals::ContractedGaussianOrbital& orb1, const Orbitals::ContractedGaussianOrbital& orb2, const Orbitals::ContractedGaussianOrbital& orb3)
	{
		int l = 0;
		for (const auto& atom4 : m_Molecule->atoms)
			for (const auto& shell4 : atom4.shells)
			{
				for (const auto& orb4 : shell4.basisFunctions)
				{
					const int kl = k * (k + 1) / 2 + l;
					
					if (ij <= kl) electronElectronIntegrals[GetElectronElectronIndex(i, j, k, l)] = getElectronElectron(&orb1, &orb2, &orb3, &orb4);

					if (++l > k) 
					{
						if (!useLotsOfMemory) ClearElectronElectronIntermediaries(); // if this is cleared, it more than doubles the execution time, but it uses a lot less memory

						return;
					}
				}

				if (!useLotsOfMemory) ClearElectronElectronIntermediaries(); // if this is cleared, it more than doubles the execution time, but it uses a lot less memory
			}
	}



	inline void IntegralsRepository::CalculateElectronElectronIntegrals23(int i, const Orbitals::ContractedGaussianOrbital& orb1)
	{
		int j = 0;
		for (const auto& atom2 : m_Molecule->atoms)
			for (const auto& shell2 : atom2.shells)
				for (const auto& orb2 : shell2.basisFunctions)
				{
					int ij = i * (i + 1) / 2 + j;

					int k = 0;
					for (const auto& atom3 : m_Molecule->atoms)
						for (const auto& shell3 : atom3.shells)
							for (const auto& orb3 : shell3.basisFunctions)
								CalculateElectronElectronIntegrals4(i, j, k++, ij, orb1, orb2, orb3);

					if (++j > i) return;
				}
	}



	void PrintMemoryInfo()
	{
		DWORD processID = GetProcessId(GetCurrentProcess());

		HANDLE hProcess;
		PROCESS_MEMORY_COUNTERS pmc;

		// Print information about the memory usage of the process.

		hProcess = OpenProcess(PROCESS_QUERY_INFORMATION |
			PROCESS_VM_READ,
			FALSE, processID);
		if (NULL == hProcess)
			return;

		if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
		{
			wchar_t buf[2048];

			swprintf(buf, sizeof(buf), L"\nProcess ID: %u\n\tPeakWorkingSetSize: %f\n\tWorkingSetSize: %f\n\tPagefileUsage: %f\n\tPeakPagefileUsage: %f\n",
				processID, pmc.PeakWorkingSetSize / (1024.*1024.), pmc.WorkingSetSize / (1024.*1024.), pmc.PagefileUsage / (1024.*1024.), pmc.PeakPagefileUsage / (1024.*1024.));

			AfxMessageBox(buf);
		}

		CloseHandle(hProcess);
	}






	void IntegralsRepository::CalculateElectronElectronIntegrals()
	{
		const int maxNr = m_Molecule->CountNumberOfContractedGaussians();
		const int maxIndex = GetElectronElectronIndex(maxNr, maxNr, maxNr, maxNr);

		electronElectronIntegrals.resize(maxIndex + 1);

		int i = 0;
		for (const auto& atom1 : m_Molecule->atoms)
			for (const auto& shell1 : atom1.shells)
				for (const auto& orb1 : shell1.basisFunctions)
					CalculateElectronElectronIntegrals23(i++, orb1);
					
		//PrintMemoryInfo();

		ClearElectronElectronMaps();

		//PrintMemoryInfo();
	}


}










