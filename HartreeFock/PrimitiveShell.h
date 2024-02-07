#pragma once

#include <vector>
#include <string>
#include <cassert>

#include "GaussianOrbital.h"
#include "ContractedGaussianOrbital.h"
#include "Vector3D.h"


namespace Orbitals {

	// a shell has its center on the atom nucleus, it contains the basis functions with this center
	class Shell
	{
	public:
		unsigned int ID = 0;
		unsigned int centerID = 0;

		virtual ~Shell() = default;

		virtual const Vector3D<double>& getCenter() const = 0;
		virtual double operator()(const Vector3D<double>& r) const = 0;
		virtual Vector3D<double> getGradient(const Vector3D<double>& r) const = 0;
		virtual double getLaplacian(const Vector3D<double>& r) const = 0;
		virtual void Normalize() = 0;
	};

	// the basis functions inside have the same center and exponent
	// this is currently not used in the program
	class PrimitiveGaussianShell : public Shell
	{
	public:
		std::vector<GaussianOrbital> basisFunctions;
		
		const Vector3D<double>& getCenter() const override;

		double getAlpha() const;

		double operator()(const Vector3D<double>& r) const override;
		Vector3D<double> getGradient(const Vector3D<double>& r) const override;
		double getLaplacian(const Vector3D<double>& r) const override;

		void Normalize() override;
	};

	// all contracted gaussian orbitals that are contained share the same center and the same set of exponents
	// for example it might contain an s contracted gaussian orbital and the three ones for px, py, pz
	class ContractedGaussianShell : public Shell
	{
	public:
		std::vector<ContractedGaussianOrbital> basisFunctions;
		
		void AddOrbitals(char type);
		void AddGaussians(double exponent);

		const Vector3D<double>& getCenter() const override;
		std::string GetShellString() const;

		unsigned int CountOrbitals(char orbital) const;
		unsigned int CountContractedOrbitals(char orbital) const;

		unsigned int CountNumberOfContractedGaussians() const;
		unsigned int CountNumberOfGaussians() const;

		double operator()(const Vector3D<double>& r) const override;
		Vector3D<double> getGradient(const Vector3D<double>& r) const override;
		double getLaplacian(const Vector3D<double>& r) const override;

		void SetCenters(const Vector3D<double>& center);
		void Normalize() override;

	private:
		static unsigned int AdjustOrbitalsCount(char orbital, unsigned int res);
		void AddOrbitalsInCanonicalOrder(unsigned int L);
	};
}
