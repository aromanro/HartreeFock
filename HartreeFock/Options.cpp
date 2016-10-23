#include "stdafx.h"
#include "Options.h"

#include "HartreeFock.h"

Options::Options()
	:
	// Molecule
	m_atom1("H"),
	m_atom2("O"),
	twoAtom1(true),
	basis(1),
	bondAngle(104.48),
	alphaElectrons(5),
	betaElectrons(5),
	// also used for charting
	XMaxBondLength(3), // Angstroms
	XMinBondLength(0.5), // not really for chart, but for calculations

	// Hartree-Fock
	restricted(true),
	alpha(0.5),
	initialGuess(0.),
	iterations(3000),

	// Computation
	nrThreads(4),
	useLotsOfMemory(true),
	numberOfPoints(80),

	// Charts
	YMaxEnergy(-1700), //eV
	YMinEnergy(-2400), //eV
	YBigTicksEnergy(8),
	YSmallTicksEnergy(2),
	XBigTicksBondLength(10),
	XSmallTicksBondLength(2),
	useSplines(false)
{
}


Options::~Options()
{
}


void Options::Load()
{
	// Molecule
	m_atom1 = theApp.GetProfileString(L"options", L"Atom1", L"H");
	m_atom2 = theApp.GetProfileString(L"options", L"Atom2", L"O");
	twoAtom1 = (1 == theApp.GetProfileInt(L"options", L"TwoAtom1", 1) ? true : false);
	basis = theApp.GetProfileInt(L"options", L"Basis", 1);
	bondAngle = GetDouble(L"BondAngle", 104.474);
	alphaElectrons = theApp.GetProfileInt(L"options", L"AlphaElectrons", 5);
	betaElectrons = theApp.GetProfileInt(L"options", L"BetaElectrons", 5);

	XMaxBondLength = GetDouble(L"XMaxBondLength", 3);
	XMinBondLength = GetDouble(L"XMinBondLength", 0.3);

	// Hartree-Fock
	restricted = (1 == theApp.GetProfileInt(L"options", L"Restricted", 1) ? true : false);
	alpha = GetDouble(L"Alpha", 0.3);
	initialGuess = GetDouble(L"InitialGuess", 1.75);
	iterations = theApp.GetProfileInt(L"options", L"MaxIterations", 3000);

	// computations
	nrThreads = theApp.GetProfileInt(L"options", L"NrThreads", 4);
	useLotsOfMemory = (1 == theApp.GetProfileInt(L"options", L"UseLotsOfMemory", 1) ? true : false);
	numberOfPoints = theApp.GetProfileInt(L"options", L"NrPoints", 80);

	// charts
	YMaxEnergy = theApp.GetProfileInt(L"options", L"YMaxEnergy", -2000);
	YMinEnergy = theApp.GetProfileInt(L"options", L"YMinEnergy", -2070);
	YBigTicksEnergy = theApp.GetProfileInt(L"options", L"YBigTicksEnergy", 7);
	YSmallTicksEnergy = theApp.GetProfileInt(L"options", L"YSmallTicksEnergy", 2);
	XBigTicksBondLength = theApp.GetProfileInt(L"options", L"XBigTicksBondLengths", 6);
	XSmallTicksBondLength = theApp.GetProfileInt(L"options", L"XSmallTicksBondLength", 2);
	useSplines = (1 == theApp.GetProfileInt(L"options", L"UseSplines", 1) ? true : false);
}


void Options::Save()
{
	// Molecule
	theApp.WriteProfileString(L"options", L"Atom1", m_atom1);
	theApp.WriteProfileString(L"options", L"Atom2", m_atom2);
	theApp.WriteProfileInt(L"options", L"TwoAtom1", twoAtom1 ? 1 : 0);
	theApp.WriteProfileInt(L"options", L"Basis", basis);
	theApp.WriteProfileBinary(L"options", L"BondAngle", (LPBYTE)&bondAngle, sizeof(double));
	theApp.WriteProfileInt(L"options", L"AlphaElectrons", alphaElectrons);
	theApp.WriteProfileInt(L"options", L"BetaElectrons", betaElectrons);

	theApp.WriteProfileBinary(L"options", L"XMaxBondLength", (LPBYTE)&XMaxBondLength, sizeof(double));
	theApp.WriteProfileBinary(L"options", L"XMinBondLength", (LPBYTE)&XMinBondLength, sizeof(double));

	// Hartree-Fock
	theApp.WriteProfileInt(L"options", L"Restricted", restricted ? 1 : 0);
	theApp.WriteProfileBinary(L"options", L"Alpha", (LPBYTE)&alpha, sizeof(double));
	theApp.WriteProfileBinary(L"options", L"InitialGuess", (LPBYTE)&initialGuess, sizeof(double));
	theApp.WriteProfileInt(L"options", L"MaxIterations", iterations);

	// computations
	theApp.WriteProfileInt(L"options", L"NrThreads", nrThreads);
	theApp.WriteProfileInt(L"options", L"UseLotsOfMemory", useLotsOfMemory ? 1 : 0);
	theApp.WriteProfileInt(L"options", L"NrPoints", numberOfPoints);

	// charts
	theApp.WriteProfileInt(L"options", L"YMaxEnergy", YMaxEnergy);
	theApp.WriteProfileInt(L"options", L"YMinEnergy", YMinEnergy);
	theApp.WriteProfileInt(L"options", L"YBigTicksEnergy", YBigTicksEnergy);
	theApp.WriteProfileInt(L"options", L"YSmallTicksEnergy", YSmallTicksEnergy);
	theApp.WriteProfileInt(L"options", L"XBigTicksBondLengths", XBigTicksBondLength);
	theApp.WriteProfileInt(L"options", L"XSmallTicksBondLength", XSmallTicksBondLength);
	theApp.WriteProfileInt(L"options", L"UseSplines", useSplines ? 1 : 0);
}

double Options::GetDouble(LPCTSTR param, double defval)
{
	double val = defval;

	UINT sz = 0;
	LPBYTE buf = NULL;

	if (theApp.GetProfileBinary(L"options", param, &buf, &sz))
	{
		if (sizeof(double) == sz) val = *((double*)buf);
		delete[] buf;
	}

	return val;
}
