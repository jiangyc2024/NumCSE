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
#include <vector>
#include <iostream>
#include <tuple>
#include <utility>

using namespace Eigen;

// convenience typedef
template <typename T, int N>
using Vector = Eigen::Matrix<T, N, 1>;

template <typename Scalar>
using broyd_history_t = std::vector<std::tuple<unsigned, Scalar, Scalar>>;

/**
 * \brief Good Broyden rank-1-update quasi-Newton method
 * Straightforward implementation for small problems
 * \param F Non-linear mapping in n dimensions
 * \param x initial guess
 * \param J initial guess for Jacobi matrix at x0
 * \param tol tolerance for termination
 */
template <typename FuncType, typename JacType, typename Scalar=double, int N=Dynamic>
std::pair<Vector<Scalar, N>, broyd_history_t<Scalar>> broyd(const FuncType &F, Vector<Scalar, N> x, JacType J, 
                                                            const Scalar tol, const unsigned maxit=20)
{
    // calculate LU factorization
    auto fac = J.lu();

    unsigned k = 1;
    Vector2d s = fac.solve(F(x));
    x -= s;
    auto sn = s.squaredNorm();
    auto f = F(x);

    // keeping a record of the convergence history
    std::vector<std::tuple<unsigned, Scalar, Scalar>> history;
    history.emplace_back(k, s.norm(), f.norm());

    while ((sn > tol*tol) && (k < 10)) {
        J = J - f*s.transpose()/sn;
        fac = J.lu();
        s = fac.solve(f);
        x -= s;
        f = F(x); sn=s.squaredNorm();
        std::cout << "Iteration " << k << ": |s| " << s.norm() << " " << f.norm() << std::endl;
        history.emplace_back(k, std::sqrt(sn), f.norm());
        k+=1; 
    }
    return std::make_pair(x, history);
}


#endif
