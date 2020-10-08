#pragma once

#include <Eigen\eigen>

#include "QuantumNumbers.h"
#include "GaussianIntegral.h"

namespace GaussianIntegrals {


    // overlap integrals are a special case of moment integrals (with L3=0, corresponding to x^0 = 1, y^0 = 1 and z^0 = 1)
    // the recurrence relations are the same but there is one more for transfer of angular momentum from the third center
    // instead of two centers there are three
    // could be implemented in general, for multipoles, but I think I'll implement it only for dipole
    // this would require only three additional matrices (for x, y, z) compared with the overlap implementation
    // to transfer into them momentum 1 for the third center (see the HSERIlib paper for details - by the way, there might be some signs wrong there)
    // also because only x or y or z are required (no products, like xy, or xz and so on) no maxQN3 should be passed as a parameter, 1 can be implied for that angular momentum
    // also a third alpha will be missing, since the third center is not about a Gaussian

    // TODO: finish implementing it

    class GaussianMoment :
        public GaussianIntegral
    {
    public:
        Eigen::MatrixXd matrixX;
        Eigen::MatrixXd matrixY;
        Eigen::MatrixXd matrixZ;

        Eigen::MatrixXd matrixX1;
        Eigen::MatrixXd matrixY1;
        Eigen::MatrixXd matrixZ1;

        double factor;

        GaussianMoment();

        GaussianMoment(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2)
            : GaussianMoment(alpha1, alpha2, center1, center2, Vector3D<double>(0, 0, 0), maxQN1, maxQN2)
        {
        }

        GaussianMoment(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Vector3D<double>& center3, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2);


        void Reset(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2)
        {
            static const Vector3D<double> zeroV(0, 0, 0);
            Reset(alpha1, alpha2, center1, center2, zeroV, maxQN1, maxQN2);
        }

        void Reset(double alpha1, double alpha2, const Vector3D<double>& center1, const Vector3D<double>& center2, const Vector3D<double>& center3, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN1, const Orbitals::QuantumNumbers::QuantumNumbers& maxQN2);

        double getMomentX(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const
        {
            return getMoment(QN1, QN2, true, false, false);
        }

        double getMomentY(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const
        {
            return getMoment(QN1, QN2, false, true, false);
        }

        double getMomentZ(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2) const
        {
            return getMoment(QN1, QN2, false, false, true);
        }

        double getMoment(const Orbitals::QuantumNumbers::QuantumNumbers& QN1, const Orbitals::QuantumNumbers::QuantumNumbers& QN2, bool momentX, bool momentY, bool momentZ) const;

    protected:
        void CalculateMoment(Eigen::MatrixXd& matrix, Eigen::MatrixXd& matrix1, double alpha1, double alpha2, double center1, double center2, double center3, unsigned int maxQN1, unsigned int maxQN2);

    };

}
