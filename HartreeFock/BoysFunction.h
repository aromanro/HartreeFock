#pragma once

#include "MathUtils.h"

namespace GaussianIntegrals {

	class BoysFunction
	{
	public:
		BoysFunction();
		~BoysFunction();

		class BoysFunctor : public MathUtils::FunctionFunctor
		{
		protected:
			double m_m;
			double m_x;
		public:
			BoysFunctor(double m, double x) : m_m(m), m_x(x) {}

			virtual double operator()(double t) const 
			{ 
				const double res = pow(t, 2. * m_m) * exp(-m_x * t * t);

				return res;  
			}
		};

		double operator()(double m, double x) const;
	};

}