#pragma once

#include "MathUtils.h"

namespace GaussianIntegrals {

	class BoysFunction
	{
	public:
		class BoysFunctor : public MathUtils::FunctionFunctor
		{
		public:
			BoysFunctor(double m, double x) : m_m(m), m_x(x) {}

			double operator()(double t) const noexcept override
			{ 
				const double res = pow(t, 2. * m_m) * exp(-m_x * t * t);

				return res;  
			}

		private:
			double m_m;
			double m_x;
		};

		double operator()(double m, double x) const noexcept;
	};

}