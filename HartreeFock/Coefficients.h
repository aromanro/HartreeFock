#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

#include <array>
#include <tuple>
#include <unordered_map>

namespace CG
{

	// Typically those are done by starting on top or bottom of the 'ladder' and by applying ladder operators 
	// you get all the needed ones recursively
	// since I need them for low L values for the purpose of the blog projects, this method should do
	// they will be calculated once and cached, so there is no big performance penalty, either

	// this one is here just to be able to use an unordered_map with a tuple of integers as a key
	template <typename... Tp> class TupleHash {
	public:
		size_t operator()(const std::tuple<Tp...>& t) const
		{
			size_t res = 1;
			std::apply([&res](auto&& ... args) { 
				auto compute = [&res](const auto& x) {
					res = 31 * res + x;
				};

				(compute(args), ...);
				}, t);

			return res;
		}
	};

	// not used yet, but I might use it in the future for various improvements (like computing the spectrum)
	// see the Wigner–Eckart theorem for details
	class Coefficients
	{
	public:
		static double Factorial(long long int n)
		{
			if (n <= 1) return 1;
			else if (n < 21) return static_cast<double>(factorialsTable[n]);

			double val = static_cast<double>(factorialsTable[20]);
			for (long int i = 21; i <= n; ++i)
				val *= i;

			return val;
		}

		double TriangleCoefficient(long long int a, long long int b, long long int c) const
		{
			if (a + b < c || a + c < b || b + c < a) return 0;

			return Factorial(a + b - c) * Factorial(a - b + c) * Factorial(-a + b + c) / Factorial(a + b + c + 1);
		}

		double CalculateClebschGordan(double j1, double j2, double j3, double m1, double m2, double m3) const
		{
			static const double eps = 1E-15;

			if (m1 < -j1 || m1 > j1 || m2 < -j2 || m2 > j2 || m3 < -j3 || m3 > j3) return 0;
			else if (j3 < abs(j1 - j2) || j3 > j1 + j2) return 0;
			else if (abs(m3 - m2 - m1) > eps) return 0;

			//if (m3 < 0) return pow(-1, j3 - j1 - j2) * CalculateClebschGordan(j1, j2, j3, -m1, -m2, -m3);
			//else if (j1 < j2) return pow(-1, j3 - j1 - j2) * CalculateClebschGordan(j2, j1, j3, m2, m1, m3);

			double val = 0;

			const double limMin = std::max<>(std::max<>(j2 - j3 - m1, j1 - j3 + m2), 0.);
			const double limMax = std::min<>(j2 + m2, std::min<>(j1 - m1, j1 + j2 - j3));
			for (long long int k = static_cast<long long int>(limMin); k <= limMax; ++k)
			{
				const double t1 = j1 + j2 - j3 - k;
				const double t2 = j1 - m1 - k;
				const double t3 = j2 + m2 - k;

				const double t4 = j3 - j2 + m1 + k;
				const double t5 = j3 - j1 - m2 + k;
				val += pow(-1, static_cast<double>(k)) / (Factorial(static_cast<unsigned long long int>(k)) * Factorial(static_cast<unsigned long long int>(t1)) * Factorial(static_cast<unsigned long long int>(t2)) *
					Factorial(static_cast<unsigned long long int>(t3)) * Factorial(static_cast<unsigned long long int>(t4)) * Factorial(static_cast<unsigned long long int>(t5)));
			}

			if (0 == val) return 0;

			return val * sqrt((2. * j3 + 1.) * Factorial(static_cast<unsigned long long int>(j3 + j1 - j2))* Factorial(static_cast<unsigned long long int>(j3 - j1 + j2)) *
				Factorial(static_cast<unsigned long long int>(j1 + j2 - j3)) / Factorial(static_cast<unsigned long long int>(j1 + j2 + j3 + 1.))) *
				sqrt(Factorial(static_cast<unsigned long long int>(j1 + m1)) * Factorial(static_cast<unsigned long long int>(j1 - m1)) *
					 Factorial(static_cast<unsigned long long int>(j2 + m2)) * Factorial(static_cast<unsigned long long int>(j2 - m2)) *
					 Factorial(static_cast<unsigned long long int>(j3 + m3)) * Factorial(static_cast<unsigned long long int>(j3 - m3)));
		}

		double CalculateWigner3j(double j1, double j2, double j3, double m1, double m2, double m3) const
		{
			static const double eps = 1E-15;

			if (m1 < -j1 || m1 > j1 || m2 < -j2 || m2 > j2 || m3 < -j3 || m3 > j3) return 0;
			else if (j3 < abs(j1 - j2) || j3 > j1 + j2) return 0;
			else if (abs(m1 + m2 + m3) > eps) return 0;

			if (abs(m1) < eps && abs(m2) < eps && abs(m3) < eps)
			{
				int intSum = static_cast<int>(j1 + j2 + j3);

				// the sum must be an even integer
				if (abs(j1 + j2 + j3 - intSum) > 0.1) return 0;
				else if (intSum % 2) return 0;
			}

			return pow(-1., j2 - j1 - m3) / sqrt(2. * j3 + 1.) * CalculateClebschGordan(j1, j2, j3, m1, m2, -m3);
		}

		double CalculateGaunt(double j1, double j2, double j3, double m1, double m2) const
		{
			int intSum = static_cast<int>(j1 + j2 + j3);
			// the sum must be an even integer (otherwise the first Wigner 3j term is zero)
			if (abs(j1 + j2 + j3 - intSum) > 0.1) return 0;
			else if (intSum % 2) return 0;

			static const double fourM_PI = 4. * M_PI;

			// either one would do, I prefer the Wigner 3j formula, although Wigner 3j is computed from Clebsch Gordan

			//return pow(-1., m2) * sqrt((2. * j1 + 1.) * (2. * j2 + 1.) / (fourM_PI * (2. * j3 + 1.))) *
			//	CalculateClebschGordan(j1, j2, j3, 0, 0, 0) * CalculateClebschGordan(j1, j2, j3, m1, -m2, m1 - m2);

			return pow(-1., m1) * sqrt((2. * j1 + 1.) * (2. * j2 + 1.) * (2. * j3 + 1.) / fourM_PI) *
				CalculateWigner3j(j1, j2, j3, 0, 0, 0) * CalculateWigner3j(j1, j2, j3, m1, -m2, m2 - m1);
		}


		void ClearCoefficientsCache()
		{
			coefficients.clear();
		}

		void PrecalculateCoefficients(int max_l)
		{
			static const double fourM_PI = 4. * M_PI;

			for (int j1 = 0; j1 <= max_l; ++j1)
				for (int j2 = 0; j2 <= max_l; ++j2)
					for (int j3 = abs(j1 - j2); j3 <= j1 + j2; j3 += 2)
					{
						//if ((j1 + j2 + j3) % 2) continue; // must be even to be different than 0, no need to check it, the for above ensures the sum is even

						for (int m1 = -j1; m1 <= j1; ++m1)
							for (int m2 = -j2; m2 <= j2; ++m2)
							{
								const double val = CalculateGaunt(j1, j2, j3, m1, m2);
								if (val == 0) continue;

								const auto ind = std::make_tuple(static_cast<int>(2. * j1), static_cast<int>(2. * j2), static_cast<int>(2. * j3), static_cast<int>(2. * m1), static_cast<int>(2. * m2));
								coefficients[ind] = val;
							}
					}
		}

		double getCoefficient(double j1, double j2, double j3, double m1, double m2) const
		{
			const auto ind = std::make_tuple(static_cast<int>(2. * j1), static_cast<int>(2. * j2), static_cast<int>(2. * j3), static_cast<int>(2. * m1), static_cast<int>(2. * m2));
			const auto res = coefficients.find(ind);
			if (res == coefficients.end())
				return 0;

			return res->second;
		}

	private:
		static std::array<unsigned long long int, 21> factorialsTable;

		std::unordered_map<std::tuple<int, int, int, int, int>, double, TupleHash<int, int, int, int, int>> coefficients;
	};

}

