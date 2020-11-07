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

        inline int delta(int i, int j)
        {
            return i == j ? 1 : 0;
        }

        // this should correspond to Step #1: Preparing the Spin-Orbital Basis Integrals
        // implemented in a hurry, needs checking

        inline Eigen::MatrixXd getSpinOrbitalFockMatrix()
        {
            Eigen::MatrixXd spinOrbitalFockMatrix(numberOfSpinOrbitals, numberOfSpinOrbitals);

            Eigen::MatrixXd FockMatrixMO = C.transpose() * h * C;

            for (int p = 0; p < numberOfSpinOrbitals; ++p)
            {
                const int hp = p / 2;

                for (int q = 0; q < numberOfSpinOrbitals; ++q)
                {
                    spinOrbitalFockMatrix(p, q) = (p % 2 == q % 2) * FockMatrixMO(hp, q / 2);
                    for (int m = 0; m < numberOfSpinOrbitals; ++m)
                    {
                        const int orb = m / 2;
                        if (orb >= occupied.size() || !occupied[orb]) continue; // only occupied

                        spinOrbitalFockMatrix(p, q) += (*m_spinOrbitalBasisIntegrals)(p, m, q, m);
                    }
                }
            }

            return spinOrbitalFockMatrix;
        }


        // Step #2: Build the Initial-Guess Cluster Amplitudes
        void InitialGuessClusterAmplitudes()
        {
            const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

            t2 = Eigen::MatrixXd::Zero(numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

            t4.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

            int indi = 0;
            for (int i = 0; i < numberOfSpinOrbitals; ++i)
            {
                const int hi = i / 2;
                if (hi >= occupied.size() || !occupied[hi]) continue; // only occupied

                int indj = 0;
                for (int j = 0; j < numberOfSpinOrbitals; ++j)
                {
                    const int hj = j / 2;
                    if (hj >= occupied.size() || !occupied[hj]) continue; // only occupied

                    int inda = 0;
                    for (int a = 0; a < numberOfSpinOrbitals; ++a)
                    {
                        const int ha = a / 2;
                        if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied

                        int indb = 0;
                        for (int b = 0; b < numberOfSpinOrbitals; ++b)
                        {
                            const int hb = b / 2;
                            if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

                            t4(indi, indj, inda, indb) = (*m_spinOrbitalBasisIntegrals)(i, j, a, b) / (eigenvals(hi) + eigenvals(hj) - eigenvals(ha) - eigenvals(hb));

                            ++indb;
                        }

                        ++inda;
                    }
                
                    ++indj;
                }

                ++indi;
            }
        }


        // Step #3: Calculate the CC Intermediates
        
        // for the formulae, see J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, J. Chem. Phys. volume 94, pp. 4334-4345 (1991)
        // "A direct product decomposition approach for symmetry exploitation in many-body methods. I. Energy calculations"

        // formulae 3 - 13

        // formulae 9 and 10 first, taus are needed in the following calculations

        void CalculateTaus()
        {
            const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

            tau.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);
            taut.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);

            int indi = 0;
            for (int i = 0; i < numberOfSpinOrbitals; ++i)
            {
                const int hi = i / 2;
                if (hi >= occupied.size() || !occupied[hi]) continue; // only occupied

                int indj = 0;
                for (int j = 0; j < numberOfSpinOrbitals; ++j)
                {
                    const int hj = j / 2;
                    if (hj >= occupied.size() || !occupied[hj]) continue; // only occupied

                    int inda = 0;
                    for (int a = 0; a < numberOfSpinOrbitals; ++a)
                    {
                        const int ha = a / 2;
                        if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied

                        int indb = 0;
                        for (int b = 0; b < numberOfSpinOrbitals; ++b)
                        {
                            const int hb = b / 2;
                            if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

                            const double term2 = t2(indi, inda) * t2(indj, indb) - t2(indi, indb) * t2(indj, inda);

                            tau(indi, indj, inda, indb) = t4(indi, indj, inda, indb) + 0.5 * term2;
                            taut(indi, indj, inda, indb) = t4(indi, indj, inda, indb) + term2;

                            ++indb;
                        }

                        ++inda;
                    }

                    ++indj;
                }

                ++indi;
            }

        }





        // end of computing intermediates



        // compute denominator arrays

        // formulae 12 and 13

        void CalculateDenominatorArrays()
        {
            const int numberOfUnoccupiedSpinOrbitals = numberOfSpinOrbitals - numberOfOccupiedSpinOrbitals;

            D2.resize(numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);
            D4.resize(numberOfOccupiedSpinOrbitals, numberOfOccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals, numberOfUnoccupiedSpinOrbitals);


            int indi = 0;
            for (int i = 0; i < numberOfSpinOrbitals; ++i)
            {
                const int hi = i / 2;
                if (hi >= occupied.size() || !occupied[hi]) continue; // only occupied

                int inda = 0;
                for (int a = 0; a < numberOfSpinOrbitals; ++a)
                {
                    const int ha = a / 2;
                    if (ha < occupied.size() && occupied[ha]) continue; // only unoccupied

                    D2(indi, inda) = f(i, i) - f(a, a);

                    int indj = 0;
                    for (int j = 0; j < numberOfSpinOrbitals; ++j)
                    {
                        const int hj = j / 2;
                        if (hj >= occupied.size() || !occupied[hj]) continue; // only occupied

                        int indb = 0;
                        for (int b = 0; b < numberOfSpinOrbitals; ++b)
                        {
                            const int hb = b / 2;
                            if (hb < occupied.size() && occupied[hb]) continue; // only unoccupied

                            D4(indi, indj, inda, indb) = f(i, i) + f(j, j) - f(a, a) - f(b, b);

                            ++indb;
                        }

                        ++indj;
                    }

                    ++inda;
                }

                ++indi;
            }

        }


    public:

        // just for checking against the MP2 energy
        double MP2EnergyFromt4() const
        {
            double result = 0;

            int indi = 0;
            for (int i = 0; i < numberOfSpinOrbitals; ++i)
            {
                const int orbi = i / 2;
                if (orbi >= occupied.size() || !occupied[orbi]) continue; // only occupied

                int indj = 0;
                for (int j = 0; j < numberOfSpinOrbitals; ++j)
                {
                    const int orbj = j / 2;
                    if (orbj >= occupied.size() || !occupied[orbj]) continue; // only occupied

                    int inda = 0;
                    for (int a = 0; a < numberOfSpinOrbitals; ++a)
                    {
                        const int orba = a / 2;
                        if (orba < occupied.size() && occupied[orba]) continue; // only unoccupied

                        int indb = 0;
                        for (int b = 0; b < numberOfSpinOrbitals; ++b)
                        {
                            const int orbb = b / 2;
                            if (orbb < occupied.size() && occupied[orbb]) continue;  // only unoccupied

                            result += (*m_spinOrbitalBasisIntegrals)(i, j, a, b) * t4(indi, indj, inda, indb);

                            ++indb;
                        }

                        ++inda;
                    }

                    ++indj;
                }

                ++indi;
            }


            return 0.25 * result;
        }
        

    private:

        int numberOfSpinOrbitals;
        int numberOfOccupiedSpinOrbitals;

        Eigen::MatrixXd f;


        Eigen::MatrixXd t2;
        
        // they should be better implemented than my own implementation
        Eigen::Tensor<double, 4> t4;
        
        GaussianIntegrals::CoupledClusterSpinOrbitalsElectronElectronIntegralsRepository* m_spinOrbitalBasisIntegrals;

        // intermediates
        Eigen::MatrixXd Fae; // unoccupied, unoccupied
        Eigen::MatrixXd Fmi; // occupied, occupied
        Eigen::MatrixXd Fme; // occupied, unoccupied

        Eigen::Tensor<double, 4> Wmnij; // occupied, occupied, occupied, occupied
        Eigen::Tensor<double, 4> Wabef; // unoccupied, unoccupied, unoccupied, unoccupied
        Eigen::Tensor<double, 4> Wmbej; // occupied, unoccupied, unoccupied, occupied

        // effective doubles (two particle excitation operators)
        Eigen::Tensor<double, 4> tau; // occupied, occupied, unoccupied, unoccupied
        Eigen::Tensor<double, 4> taut; // occupied, occupied, unoccupied, unoccupied

        // denominator arrays
        Eigen::MatrixXd D2;
        Eigen::Tensor<double, 4> D4;
    };

}

