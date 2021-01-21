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

// convenience typedef
template <typename T, int N>
using Vector = Eigen::Matrix<T, N, 1>;

/**
 * \brief Good Broyden rank-1-update quasi-Newton method
 * Naive and inefficient implementation for small problems
 * \param F Non-linear mapping in n dimensions
 * \param x initial guess
 * \param J initial guess for Jacobi matrix at x0
 * \param tol tolerance for termination
 * \param callback to be run in every iteration step
 */
template <typename FuncType, typename JacType, typename Scalar,
          int N = Eigen::Dynamic, typename CB = void_cb>
Vector<Scalar, N> broyd(FuncType&& F, Vector<Scalar, N> x, JacType J,
                        const Scalar tol, const unsigned maxit = 20,
                        CB callback = nullptr) {
  // Compute first quasi-Newton update
  Vector<Scalar, N> s = -J.lu().solve(F(x));
  // Compute $\Vx^{(1)}$
  x += s;
  auto f = F(x);
  for (unsigned int k = 1; ((s.norm() > tol * x.norm()) && (k < maxit)); ++k) {
    // Rank-1 update of Jacobian \eqref{eq:broydenqn}
    J += f * s.transpose() / s.squaredNorm();
    // Next quasi-Newton correction
    s = -J.lu().solve(f);
    // Update of iterate
    x += s;
    // Evaluate F at next iterate
    f = F(x);
    // Optional output 
    if (callback != nullptr) {
      callback(k, x, f, s);
    }
  }
  return x;
}

/* SAM_LISTING_BEGIN_1 */
template <typename FuncType, typename JacType, typename Scalar = double,
          int N = Eigen::Dynamic, typename CB = void_cb>
Vector<Scalar, N> upbroyd(const FuncType &F, Vector<Scalar, N> x, JacType J,
                          const Scalar tol, const unsigned maxit = 20,
                          CB callback = nullptr) {
  // Calculate LU factorization of initial Jacobian once
  auto fac = J.lu();
  // First quasi-Newton correction $\cob{\Delta\Vx^{(0)} := -\VJ_0^{-1}F(\Vx^{(0)})}$
  Vector<Scalar, N> s = -fac.solve(F(x));
  // Store the first quasi-Newton correction $\cob{\Delta\Vx^{(0)}}$
  std::vector<Vector<Scalar, N>> dx{s};
  // First update of iterates: $\cob{\Vx^{(1)} := \Vx^{(0)} + \Delta\Vx^{(0)}}$
  x += s;
  auto f = F(x); // Here $\cob{ = F(\Vx^{(1)}}$
  // Empty sequence of simplified quasi-Newton corrections 
  std::vector<Vector<Scalar, N>> dxs{};
  // Store denominators $\cob{\N{\Vx^{(\ell)}}^{2}+(\Delta\Vx^{(\ell)})^{\top}\cop{\Delta\overline{\Vx}^{(\ell+1)}}}$
  std::vector<Scalar> den{};
  // callback once before we start the algorithm
  callback(0, x, f, s);
  for (unsigned int k = 1; ((s.norm() > tol * x.norm()) && (k < maxit)); ++k) {
    // Compute $\cob{\VJ_{0}^{-1}F(\Vx^{(k)})}$
    s = fac.solve(f);
    // Compute next simplified quasi-Newton correction recursively 
    Vector<Scalar, N> ss = s;
    for (unsigned int l = 1; l < k; ++l) {
      ss -= dxs[l-1] * (dx[l - 1].dot(ss)) / den[l - 1];
    }
    // Store next denominator
    den.push_back(dx[k - 1].squaredNorm() + dx[k - 1].dot(ss));
    // Store current simplified quasi-Newton correction
    dxs.push_back(ss);
    // Compute next quasi-Newton correction recursively 
    for (unsigned int l = 0; l < k; ++l) {
      s -= dxs[l] * (dx[l].dot(s)) / den[l];
    }
    s *= (-1.0);
    dx.push_back(s);
    // Compute next iterate
    x += s;
    // Evaluation F at next iterate
    f = F(x);
    if (callback != nullptr) { // NOLINT
      callback(k, x, f, s);
    }
  }
  return x;
}
/* SAM_LISTING_END_1 */
#endif
