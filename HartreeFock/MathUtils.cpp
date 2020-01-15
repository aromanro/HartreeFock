#include "stdafx.h"
#include "MathUtils.h"



namespace GaussianIntegrals {


	constexpr unsigned long long int factorial(long int n)
	{
		return n ? n * factorial(n - 1) : 1;
	}

	constexpr unsigned long long int double_factorial(long int n)
	{
		return (n > 1) ? n * double_factorial(n - 2) : 1;
	}

	// filled at compile time
	std::array<unsigned long long int, 21> MathUtils::factorialsTable{
		factorial(0),
		factorial(1),
		factorial(2),
		factorial(3),
		factorial(4),
		factorial(5),
		factorial(6),
		factorial(7),
		factorial(8),
		factorial(9),
		factorial(10),
		factorial(11),
		factorial(12),
		factorial(13),
		factorial(14),
		factorial(15),
		factorial(16),
		factorial(17),
		factorial(18),
		factorial(19),
		factorial(20)
	};


	std::array<unsigned long long int, 21> MathUtils::doubleFactorialsTable{
		double_factorial(0),
		double_factorial(1),
		double_factorial(2),
		double_factorial(3),
		double_factorial(4),
		double_factorial(5),
		double_factorial(6),
		double_factorial(7),
		double_factorial(8),
		double_factorial(9),
		double_factorial(10),
		double_factorial(11),
		double_factorial(12),
		double_factorial(13),
		double_factorial(14),
		double_factorial(15),
		double_factorial(16),
		double_factorial(17),
		double_factorial(18),
		double_factorial(19),
		double_factorial(20)
	};

}