#pragma once
#include "RestrictedHartreeFock.h"

#include "CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository.h"

// implementation guided by this tutorial: https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2305

namespace HartreeFock {

    class RestrictedCCSD :
        public RestrictedHartreeFock
    {
    public:
        RestrictedCCSD(int iterations = 3000);
        virtual ~RestrictedCCSD();

        virtual void Init(Systems::Molecule* molecule) override;

    private:

        // this should correspond to Step #1: Preparing the Spin-Orbital Basis Integrals
        // implemented in a hurry, needs checking

        inline Eigen::MatrixXd getSpinOrbitalFockMatrix()
        {
            const int numberOfSpinOrbitals = 2 * numberOfOrbitals;
            Eigen::MatrixXd spinOrbitalFockMatrix(numberOfSpinOrbitals, numberOfSpinOrbitals);

            for (int p = 0; p < numberOfSpinOrbitals; ++p)
                for (int q = 0; q < numberOfSpinOrbitals; ++q)
                {
                    spinOrbitalFockMatrix(p, q) = (p % 2 == p % 2) * h(p / 2, q / 2);
                    for (int m = 0; m < numberOfSpinOrbitals; ++m)
                        spinOrbitalFockMatrix(p, q) += (*m_spinOrbitalBasisIntegrals)(p, m, q, m);
                }

            return spinOrbitalFockMatrix;
        }


        GaussianIntegrals::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository* m_spinOrbitalBasisIntegrals;
    };

}

