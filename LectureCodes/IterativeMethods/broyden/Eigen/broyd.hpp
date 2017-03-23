///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#ifndef BROYD_HPP
#define BROYD_HPP

#include <Eigen/Dense>
#include "void_cb.hpp"

using namespace Eigen;

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
 * \param callback to be run in every iteration step
 */
template <typename FuncType, typename JacType, typename Scalar=double,
          int N=Dynamic, typename CB=void_cb>
Vector<Scalar, N> broyd(const FuncType F, Vector<Scalar, N> x, JacType J, const Scalar tol, 
                        const unsigned maxit=20, CB callback=nullptr)
{
    // calculate LU factorization
    auto fac = J.lu();

    unsigned k = 1;
    Vector2d s = fac.solve(F(x));
    x -= s;
    auto sn = s.squaredNorm();
    auto f = F(x);

    while ((sn > tol*tol) && (k < 10)) {
        J = J - f*s.transpose()/sn;
        fac = J.lu();
        s = fac.solve(f);
        x -= s;
        f = F(x); sn=s.squaredNorm();

        if (callback != nullptr)
            callback(k, x, f, s);
        k+=1; 
    }
    return x;
}


#endif
