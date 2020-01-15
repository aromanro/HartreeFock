#pragma once

#include <vector>

namespace GaussianIntegrals {

	class BoysFunctions
	{
	public:
		std::vector<double> functions;

		void GenerateBoysFunctions(int maxM, double T);
	};

}
