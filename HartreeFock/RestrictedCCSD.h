#pragma once
#include "RestrictedHartreeFock.h"

#include "CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository.h"

#include <unsupported/Eigen/CXX11/Tensor>

#include <eigen/eigen>

// implementation guided by this tutorial: https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2305

namespace HartreeFock {

    class RestrictedCCSD :
        public RestrictedHartreeFock
    {
    public:
        RestrictedCCSD(int iterations = 3000);
        virtual ~RestrictedCCSD();

        virtual void Init(Systems::Molecule* molecule) override;

        void InitCC()
        {
            m_spinOrbitalBasisIntegrals->Compute(integralsRepository, C);
            InitialGuessClusterAmplitudes();
            f = getSpinOrbitalFockMatrix();
        }

    private:

        // this should correspond to Step #1: Preparing the Spin-Orbital Basis Integrals
        // implemented in a hurry, needs checking

        inline Eigen::MatrixXd getSpinOrbitalFockMatrix()
        {
            Eigen::MatrixXd spinOrbitalFockMatrix(numberOfSpinOrbitals, numberOfSpinOrbitals);

            Eigen::MatrixXd FockMatrixMO = C.transpose() * h * C;

            for (int p = 0; p < numberOfSpinOrbitals; ++p)
                for (int q = 0; q < numberOfSpinOrbitals; ++q)
                {
                    spinOrbitalFockMatrix(p, q) = (p % 2 == q % 2) * FockMatrixMO(p / 2, q / 2);
                    for (int m = 0; m < numberOfSpinOrbitals; ++m)
                        spinOrbitalFockMatrix(p, q) += (*m_spinOrbitalBasisIntegrals)(p, m, q, m);
                }

            return spinOrbitalFockMatrix;
        }


        // Step #2: Build the Initial-Guess Cluster Amplitudes
        void InitialGuessClusterAmplitudes()
        {
            t2 = Eigen::MatrixXd::Zero(numberOfSpinOrbitals, numberOfSpinOrbitals);

            t4.resize(numberOfSpinOrbitals, numberOfSpinOrbitals, numberOfSpinOrbitals, numberOfSpinOrbitals);

            for (int i = 0; i < numberOfSpinOrbitals; ++i)
                for (int j = 0; j < numberOfSpinOrbitals; ++j)
                    for (int a = 0; a < numberOfSpinOrbitals; ++a)
                        for (int b = 0; b < numberOfSpinOrbitals; ++b)
                            t4(i, j, a, b) = (*m_spinOrbitalBasisIntegrals)(i, j, a, b) / (eigenvals(i) + eigenvals(j) - eigenvals(a) - eigenvals(b));
        }


    public:

        // just for checking against the MP2 energy
        double MP2EnergyFromt4() const
        {
            double result = 0;

            for (int i = 0; i < numberOfSpinOrbitals; ++i)
            {
                if (i >= occupied.size() || !occupied[i]) continue; // only occupied

                for (int j = 0; j < numberOfSpinOrbitals; ++j)
                {
                    if (j >= occupied.size() || !occupied[j]) continue; // only occupied

                    for (int a = 0; a < numberOfSpinOrbitals; ++a)
                    {
                        if (a < occupied.size() && occupied[a]) continue; // only unoccupied

                        for (int b = 0; b < numberOfSpinOrbitals; ++b)
                        {
                            if (b < occupied.size() && occupied[b]) continue;  // only unoccupied

                            result += (*m_spinOrbitalBasisIntegrals)(i, j, a, b) * t4(i, j, a, b);
                        }
                    }
                }
            }


            return 0.25 * result;
        }
        

    private:

        int numberOfSpinOrbitals;

        Eigen::MatrixXd f;

        Eigen::MatrixXd t2;
        
        // they should be better implemented than my own implementation
        Eigen::Tensor<double, 4> t4;
        
        GaussianIntegrals::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository* m_spinOrbitalBasisIntegrals;
    };

}

