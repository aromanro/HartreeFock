#include "stdafx.h"
#include "Atom.h"

namespace Systems {

	Atom::Atom(unsigned int nrZ, int nrElectrons)
		: Z(nrZ), ID(0)
	{
		if (-1 == nrElectrons) electronsNumber = Z;
		else electronsNumber = nrElectrons;
	}

	unsigned int AtomWithShells::CountNumberOfContractedGaussians() const
	{
		unsigned int res = 0;

		for (auto const &shell : shells)
			res += shell.CountNumberOfContractedGaussians();

		return res;
	}

	unsigned int AtomWithShells::CountNumberOfGaussians() const
	{
		unsigned int res = 0;

		for (auto const &shell : shells)
			res += shell.CountNumberOfGaussians();

		return res;
	}


}