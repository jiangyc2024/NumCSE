//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>


using namespace Eigen;

/* @brief Steffensen's method
 * @param[in] f Function handle
 * @param[in] x0 Initial guess
 * @param[out] x All estimations returned by Steffensen's iterations until convergence
 */
template <class function>
VectorXd steffensen(function&& f, double x0) {

    VectorXd x(1); x(0) = x0;

// TODO: Steffensen's method

    return x;
}

int main() {
    // Initialization
    // Different definitions of $f$ with the same zero:
    auto f = [] (double x) { return x*std::exp(x)-1; };
//    auto f = [] (double x) { return std::exp(x)-1/x; };
//    auto f = [] (double x) { return x-std::exp(-x); };
    double x0 = 1;
    double x_star = 0.567143290409784; // From Matlab: x_star = fzero(f,x0);

    // Steffensen's method
    VectorXd x = steffensen(f, x0);

    // Compute errors
    unsigned n = x.size();
    VectorXd residuals(n), errs(n), log_errs(n);
    for(unsigned i=0; i<n; ++i) {
        residuals(i) = f(x(i));
        errs(i) = std::abs(x(i)-x_star);
        log_errs(i) = std::log(errs(i));
    }
    VectorXd ratios = VectorXd::Zero(n);
    for(unsigned i=2; i<n; ++i) {
        ratios(i) = (log_errs(i)   - log_errs(i-1)) /
                    (log_errs(i-1) - log_errs(i-2));
    }

    // Print output
    std::cout << "x" << "\t" << "errors" <<" \t" << "ratios" << std::endl;
    for(unsigned i=0; i<n; ++i) {
        std::cout << x(i) << "\t" << errs(i) << "\t" << ratios(i) << std::endl;
    }

}
