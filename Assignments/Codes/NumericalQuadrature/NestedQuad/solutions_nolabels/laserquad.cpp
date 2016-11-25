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
    // Loop over all nodes/weights pairs
    for(int i = 0; i < Q.weights.size(); ++i) {
        I += f( (Q.nodes(i) + 1) * (b - a) / 2 + a ) * Q.weights(i);
    }
    return I * (b - a) / 2.;
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

//#if SOLUTION
    // We define an auxiliary function of x defined
    // as f\_y := [\&y] (double x) \{ return f(x,y) \}, where we fix y.
    // We define the function $g$ as the function of $y$
    // $g(y) := \int_0^{1-y} f_y(x) dx$, which is $= \int_0^{1-y} I(x,y) dx$
    auto g = [&f, &Q] (double y) { 
      return evalgaussquad(0, 1-y, [&f, &y] (double x) { return f(x,y); }, Q); 
    };
    // We integrate the function g over y from 0 to 1
    return evalgaussquad(0, 1, g, Q);

    /* EQUIVALENT: Loop based, copy-and-paste implementation
    // Integration over y from 0 to 1 of $g(y) := \int_0^{1-y} I(x,y) dx$
    double I = 0;
    double a = 0., b  = 1.;
    for(int i = 0; i < Q.weights.size(); ++i) {
        // Find out the y at which we are
        double y = (qr.x(i) + 1) * (b - a) / 2 + a;
        // Define $f_y(x)$ (y is fixed and f\_y is a function of x)
        auto f_y = [&f, &y] (double x) { return f(x,y); };
        // Compute g(y) as \int_0^{1-y} I(x,y) dx
        I += evalgaussquad(0, 1-y, f_y, Q) * Q.weights(i);
    }
    // Rescale interval
    return I * (b - a) / 2.;
    */

