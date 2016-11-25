//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iomanip>
#include <iostream>
#include <cmath>

#include "gaussquad.hpp"

//! @brief Compute $\int_a^b f(x) dx \approx \sum w_i f(x_i)$ (with scaling of w and x)
//! @tparam func template type for function handle f (e.g. lambda func.)
//! @param[in] a left boundary in [a,b]
//! @param[in] b right boundary in [a,b]
//! @param[in] f integrand
//! @param[in] N number of quadrature points
//! @return Approximation of integral $\int_a^b f(x) dx$
template <class Function>
double evalgaussquad(const double a, const double b,
                     const Function& f,
                     const QuadRule& Q) {
    double I = 0;
    // TODO: Compute the qudrature and scale the interal according to interval $[a,b]$
    return I;
}

//! @brief Compute double integral $\int_\Delta f(x,b) dx dy$.
//! Use nested Gauss quadrature.
//! @tparam func Template type for function handle f (e.g. lambda func.),
//! having operator (double x, double y) -> double
//! @param[in] f integrand, f(x,y) must be defined
//! @param[in] N number of quadrature points (in each direction)
//! @return Approximation of integral $\int_\Delta f(x,b) dx dy$
template <class Function>
double gaussquadtriangle(const Function& f, const unsigned N) {
    // Get nodes/weights for integral over dx and dy
    QuadRule Q;
    gaussquad(N, Q);

    // TODO: Compute double integral
    return 0;
}

int main() {
    // Parameters
    const double alpha = 1, p = 0, q = 0;

    // Laser beam intensity
    auto I = [alpha, p, q] (double x, double y) { 
      return std::exp(- alpha * ( (x-p)*(x-p) + (y-q)*(y-q) ) ); 
    };
    
    // Max num of Gauss points to use (in each direction)
    const unsigned max_N = 12;
    // "Exact" integral
    const double I_ex = 0.366046550000405;

    // Observed: exponential convergence (as exepcted)
    std::cout << std::setw(3) << "N"
              << std::setw(15) << "I_approx"
              << std::setw(15) << "error"
              << std::endl;
    for(unsigned N = 1; N < max_N; ++N) {

        double I_approx = gaussquadtriangle(I, N);

        std::cout << std::setw(3) << N
                  << std::setw(15) << I_approx
                  << std::setw(15) << std::abs(I_ex - I_approx)
                  << std::endl;
    }

    return 0;
}
