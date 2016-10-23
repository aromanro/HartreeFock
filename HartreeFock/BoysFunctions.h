#pragma once

#include <vector>

namespace GaussianIntegrals {

	class BoysFunctions
	{
	public:
		BoysFunctions();
		~BoysFunctions();

		std::vector<double> functions;

		void GenerateBoysFunctions(int maxM, double T);
	};

}
