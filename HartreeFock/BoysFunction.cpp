#include "stdafx.h"
#include "BoysFunction.h"


namespace GaussianIntegrals {

	BoysFunction::BoysFunction()
	{
	}


	BoysFunction::~BoysFunction()
	{
	}

	double BoysFunction::operator()(double m, double x) const
	{
		if (0. == m && 0. == x) return 1.;

		double res = 0;

		const double t = 1E-10;
		if (abs(x) < t) x = t;

		//BoysFunctor functor(m, x);
		//return MathUtils::ChebyshevIntegral(1E-15, 1000, functor) / 2.;

		if (x > 160) 
			return static_cast<double>(MathUtils::DoubleFactorial(static_cast<unsigned int>(2 * m - 1)) / pow(2, m + 1) * sqrt(M_PI / pow(x, 2 * m + 1)));

		MathUtils::IncompleteGamma(m + 1. / 2., x, res);

		return res / (2. * pow(x, m + 1. / 2.));
	}
}
