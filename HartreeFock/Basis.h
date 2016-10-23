#pragma once

#include <string>
#include <vector>

#include "Atom.h"

namespace Chemistry {

	class Basis
	{
	public:
		std::vector<Systems::AtomWithShells> atoms;

		Basis();
		~Basis();

		void Load(std::string fileName);
		void Save(const std::string& name);
		void Normalize();
	};

}