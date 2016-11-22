//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>

/*!
 * \brief quadU Function implementing weighted Gauss quadrature
 * \param f Integrand, function handle
 * \param n Number of nodes
 * \return Value of integral
 */
template<typename Function>
double quadU(const Function& f, const unsigned n) {
    double q = 0, w, xi;
    for (unsigned j = 0; j < n; j++) {
        w = M_PI/(n+1) * std::pow(std::sin((j+1)*M_PI/(n+1)), 2);
        xi = std::cos( (j+1.)/(n+1)*M_PI );
        q += w*f(xi);
    }
    return q;
}
int main(){
    // Integrand
    auto f = [](double t) { return 1. / (2+std::exp(3*t)); };

    // "Exact" value of integral
    const double exact = 0.483296828976607;

    // Store error
    std::vector<double> e(25);
    // Header
    std::cout << std::setw(5) << "n"
              << std::setw(20) << "Error"
              << std::setw(20) << "q"
              << std::endl;
    for (unsigned int n = 0; n < 25; n++) {
        // Compute error
        e[n] = std::abs(quadU(f,n+1)-exact);
        if (n > 1)
            // Print table
            std::cout << std::setw(5) << n
                      << std::setw(20) << e[n]
                      << std::setw(20) << e[n]/e[n-1]
                      << std::endl;
    }
    // After the 18th iteration only roundoff error is present.

    return 0;
}
