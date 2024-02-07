#pragma once

#include <string>

namespace Chemistry {

	class ChemUtils
	{
	public:
		static unsigned int GetZForAtom(const std::string& name);
		static std::string GetAtomNameForZ(unsigned int Z);

	private:
		static const char* atoms[];
	};

}