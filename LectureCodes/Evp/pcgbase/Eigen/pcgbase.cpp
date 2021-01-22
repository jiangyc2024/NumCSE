///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting pcgbase.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "pcgbase.hpp"

#include <Eigen/Dense>
#include <vector>

std::tuple<VectorXd, std::vector<double>, std::vector<VectorXd>> pcgbase(

    std::function<VectorXd(VectorXd)> evalA,
    VectorXd b, double tol, unsigned int maxit,
    std::function<VectorXd(VectorXd)> invB, VectorXd x) {

    // EvalA must pass a handle implementing A*x
    // invB is to be a handle to a function providing the action of the
    // preconditioning on a vector. The other arguments like for MATLAB's
    // pcg.

    std::vector<double> rn;
    std::vector<VectorXd> xk;

    VectorXd r = b - evalA(x);
    double rho = 1;

    xk.push_back(x);

    for (int i = 0; i < maxit; i++) {
        VectorXd y = invB(r);
        double rho_old = rho;
        rho = r.transpose() * y;
        rn.push_back(rho);

        VectorXd p(y.size());
        double rho0 = rho;
        if (i == 0) {
            p = y;
            rho0 = rho;
        }

        else if (rho < rho0 * tol) {
            break;
        }

        else {
            double beta = rho/rho_old;
            p = y + beta * p;
        }

        VectorXd q = evalA(p);
        double alpha = rho/(p.transpose() * q);
        x += alpha * p;
        r -= alpha * q;

        xk.push_back(x);
    }

    // Make output
    return std::make_tuple(x, rn, xk);
}
