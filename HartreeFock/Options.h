#pragma once

class Options
{
public:
	Options();
	~Options();

	void Load();
	void Save();


	// Molecule

	CString m_atom1;
	CString m_atom2;
	bool twoAtom1;

	int basis; // 0 - STO3G, 1 - STO6G
	double bondAngle;

	int alphaElectrons;
	int betaElectrons;

	// also used for charting
	double XMaxBondLength;
	double XMinBondLength; // not really for chart, but for calculations

	// Hartree-Fock

	bool restricted;
	double alpha;
	double initialGuess;
	int iterations;
	bool addAsymmetry;
	double asymmetry;

	// Computation

	int nrThreads;
	bool useLotsOfMemory;
	int numberOfPoints;

	// Charts

	int YMaxEnergy; // eV
	int YMinEnergy; // eV
	unsigned int YBigTicksEnergy;
	unsigned int YSmallTicksEnergy;

	unsigned int XBigTicksBondLength;
	unsigned int XSmallTicksBondLength;

	bool useSplines;


	int DisplayHOMOEnergy; // 0 - groundstate energy, 1 - HOMO energy, 2 - binding energy
protected:
	static double GetDouble(LPCTSTR param, double defval);
};

