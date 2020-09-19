#include "algcg.hpp"

#include <Eigen/Dense>

VectorXd cg(std::function<VectorXd(VectorXd)> evalA,
            VectorXd b, VectorXd x, double tol, unsigned int maxit) {
    // x supplies the initial guess,
    // maxit the maximal number of CG steps.
    // evalA must pass a handle to a function realising A*x

    VectorXd r = b - evalA(x);
    double rho = 1;
    double n0 = r.norm();
    VectorXd p;

    for (int i = 0; i < maxit; i++) {
        double rho_old = rho;
        rho = r.transpose() * r;

        if (i == 0) {
            p = r;
        }

        else {
            double beta = rho / rho_old;
            p = r + beta * p;
        }

        VectorXd q = evalA(p);
        double alpha = rho / (p.transpose() * q);

        x += alpha * p;  // Update approximate solution

        // Termination, see Rem. cgterm
        if ((b - evalA(x)).norm() <= tol * n0) {
            // ^ MAYBE the A * x should be evalA(x)? Since
            // evalA(x) should "realise A * x"?
            break;
        }

        r -= alpha * q;  // update of residual
    }

    return x;
}