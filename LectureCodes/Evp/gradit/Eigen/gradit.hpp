///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting gradit.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson, Vivienne Langen
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>

using namespace Eigen;

// (any assertions for template parameters) 
template <typename Matrix, typename Vector> 
Vector gradit(Matrix A, Vector b, Vector x, double rtol, double atol, unsigned int maxit) {
    // Residual
    Vector r = b - A * x;

    for (unsigned int k = 0; k < maxit; k++) {
        Vector p = A * r;

        double numerator = (r.transpose() * r);
        double denominator = (r.transpose() * p);
        double ts = numerator / denominator;  // cf. tmin

        x += ts * r;

        double cn = std::abs(ts) * r.norm();  // norm of correction

        if (cn < rtol * x.norm() || cn < atol) {
            break;
            }
        
        r -= ts * p;
        }
    return x;
    }
