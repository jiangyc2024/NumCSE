///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting gradit.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////


#include <Eigen/Dense>

using namespace Eigen;

VectorXd gradit(MatrixXd A, VectorXd b, VectorXd x, double rtol, double atol, unsigned int maxit) {
    // Residual 
    VectorXd r = b - A * x;

    for (int k = 0; k < maxit; k++) {
        VectorXd p = A * r;

        double numerator = (r.transpose() * r);
        double denominator = (r.transpose() * p);
        double ts = numerator/denominator; // cf. tmin
        
        x = x + ts * r;
        
        double cn = std::abs(ts) * r.norm(); // norm of correction

        if (cn < rtol * x.norm() || cn < atol) {
            break;
        }

        r = r - ts * p;
    }

    return x;
}