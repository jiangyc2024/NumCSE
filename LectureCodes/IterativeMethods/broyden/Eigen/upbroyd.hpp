///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#ifndef UPBROYD_HPP
#define UPBROYD_HPP

#include <Eigen/Dense>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

// convenience typedef
template <typename T, int N>
using Vector = Eigen::Matrix<T, N, 1>;

/**
 * \brief Good Broyden rank-1-update quasi-Newton method
 * Straightforward implementation for small problems
 * \param F Non-linear mapping in n dimensions
 * \param x initial guess
 * \param J initial guess for Jacobi matrix at x0
 * \param tol tolerance for termination
 */

#endif
