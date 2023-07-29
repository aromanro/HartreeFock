#include "stdafx.h"
#include "Coefficients.h"

namespace CG
{
	template <unsigned long long int n> constexpr unsigned long long int factorial() 
	{
		if constexpr (n > 1)
			return n * factorial<n - 1>();
		
		return 1;
	}


	// filled at compile time
	std::array<unsigned long long int, 21> Coefficients::factorialsTable{
		factorial<0>(),
		factorial<1>(),
		factorial<2>(),
		factorial<3>(),
		factorial<4>(),
		factorial<5>(),
		factorial<6>(),
		factorial<7>(),
		factorial<8>(),
		factorial<9>(),
		factorial<10>(),
		factorial<11>(),
		factorial<12>(),
		factorial<13>(),
		factorial<14>(),
		factorial<15>(),
		factorial<16>(),
		factorial<17>(),
		factorial<18>(),
		factorial<19>(),
		factorial<20>()
	};

}
