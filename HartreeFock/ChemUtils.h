#pragma once

#include <string>

namespace Chemistry {

	class ChemUtils
	{
	protected:
		static const char* atoms[];

	public:
		static unsigned int GetZForAtom(const std::string& name);
		static std::string GetAtomNameForZ(unsigned int Z);
	};

}