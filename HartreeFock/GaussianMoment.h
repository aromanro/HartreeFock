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
    };

}
