//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <Eigen/Dense>

#include <figure.hpp>

#include "trigpolyvalequid.hpp"

using namespace Eigen;

/*!
 * \brief plot_basis Plot the shifted basis polynomials.
 * \param n $2*n+1$ will be the number of basis polynomials.
 */
void plot_basis(int n) {
    // Mesh size
    const int M = 1e3;

    // TODO: plot the basis $b_n$.
}

/*!
 * \brief trigIpL Compute $\lambda(n)$.
 *
 * \param[in] n $2*n+1$ will be the number of basis polynomials.
 * \return Value $\lambda(n)$.
 */
double trigIpL(std::size_t n) {
    double ret = 0;

    // TODO: implement the function returning $\lambda(n)$
}

int main() {

    /// PART 1
    const int n = 5;
    plot_basis(n);

    /// PART 2
    const int s = 11;

    std::cout << std::setw(s) << "2^k"
              << std::setw(s) << "lambda(k)" << std::endl;

    for(unsigned int i = 1 << 2; i < (1 << 15); i = i << 1) {
        double l = trigIpL(i);
        std::cout << std::setprecision(3)
                  << std::setw(s) << i
                  << std::setw(s) << l << std::endl;

    }
}
