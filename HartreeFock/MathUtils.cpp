#include "stdafx.h"
#include "MathUtils.h"



namespace GaussianIntegrals {

	template <long long int n> constexpr unsigned long long int double_factorial()
	{
		if constexpr (n > 0)
			return static_cast<unsigned long long int>(n) * double_factorial<n - 2>();

		return 1;
	}

	std::array<unsigned long long int, 21> MathUtils::doubleFactorialsTable {
			double_factorial<0>(),
			double_factorial<1>(),
			double_factorial<2>(),
			double_factorial<3>(),
			double_factorial<4>(),
			double_factorial<5>(),
			double_factorial<6>(),
			double_factorial<7>(),
			double_factorial<8>(),
			double_factorial<9>(),
			double_factorial<10>(),
			double_factorial<11>(),
			double_factorial<12>(),
			double_factorial<13>(),
			double_factorial<14>(),
			double_factorial<15>(),
			double_factorial<16>(),
			double_factorial<17>(),
			double_factorial<18>(),
			double_factorial<19>(),
			double_factorial<20>()
	};

}