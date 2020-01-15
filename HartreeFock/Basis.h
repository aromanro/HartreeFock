#pragma once

#include <string>
#include <vector>

#include "Atom.h"

namespace Chemistry {

	class Basis
	{
	public:
		std::vector<Systems::AtomWithShells> atoms;

		void Load(const std::string& fileName);
		void Save(const std::string& name);
		void Normalize();
	};

}