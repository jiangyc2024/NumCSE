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
#include <vector>
#include <iostream>
#include <tuple>
#include <utility>

using namespace Eigen;

// convenience typedef
template <typename T, int N>
using Vector = Eigen::Matrix<T, N, 1>;

template <typename Scalar>
using upbroyd_history_t = std::vector<std::tuple<unsigned, Scalar, Scalar, Scalar>>;

/**
 * \brief Good Broyden rank-1-update quasi-Newton method
 * Straightforward implementation for small problems
 * \param F Non-linear mapping in n dimensions
 * \param x initial guess
 * \param J initial guess for Jacobi matrix at x0
 * \param tol tolerance for termination
 */
template <typename FuncType, typename JacType, typename Scalar=double, int N=Dynamic, typename CB=void(*)()>
Vector<Scalar, N> upbroyd(const FuncType &F, Vector<Scalar, N> x, JacType J, 
                          const Scalar tol, const unsigned maxit=20, CB callback=nullptr)
{
    // calculate LU factorization
    auto fac = J.lu();

    unsigned k = 1;
    Vector2d s = fac.solve(F(x));
    x -= s;
    auto sn = s.squaredNorm();

    // update vector
    std::vector<VectorXd> dx = {s};
    std::vector<Scalar> dxn = {sn};

    // callback once before we start the algorithm
    callback(k, x, F(x), s, Vector<Scalar, N>::Zero(), dxn);

    while ((sn > tol*tol) && (k < maxit)) {
        Vector<Scalar, N> w = fac.solve(F(x));
        for (unsigned l=1; l<k; ++l) {
            w = w + dx[l]*(dx[l-1].dot(w))/dxn[l-1];
        }
        auto z = s.dot(w);
        s = (1+z/(sn-z))*w;
        sn = s.squaredNorm();
        dx.push_back(s);
        dxn.push_back(sn);
        x -= s;
        if (callback != nullptr)
            callback(k, x, F(x), s, w, dxn);
        ++k;
    }
    return x;
}

#endif
