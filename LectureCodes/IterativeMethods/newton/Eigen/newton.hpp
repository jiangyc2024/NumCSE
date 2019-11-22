///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair <hiptmair@sam.math.ethz.ch>,
//             Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#ifndef NEWTON_HPP
#define NEWTON_HPP

#include "void_cb.hpp"
#include <Eigen/Dense>

using namespace Eigen;

// convenience typedef
template <typename T, int N> using Vector = Eigen::Matrix<T, N, 1>;

/* SAM_LISTING_BEGIN_0 */
template <typename FuncType, typename JacType, typename VecType>
VecType newton(FuncType &&F, JacType &&DFinv, VecType x, const double rtol,
               const double atol) {
  // Note that the vector x passes both the initial guess and also
  // contains the iterates
  VecType s(x.size()); // Vector for Newton corrections
  // Main loop
  do {
    s = DFinv(x, F(x)); // compute Newton correction
    x -= s;             // compute next iterate
  }
  // correction based termination (relative and absolute)
  while ((s.norm() > rtol * x.norm()) && (s.norm() > atol));
  return x;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename FuncType, typename JacType, typename Scalar, int N = Dynamic,
          typename CB = void_cb>
Vector<Scalar, N> newton_ext(FuncType &F, JacType &DFinv,
                             const Vector<Scalar, N> x0, const double rtol,
                             const double atol, CB callback = nullptr) {
  Vector<Scalar, N> x = x0;
  Vector<Scalar, N> s;

  do {
    s = DFinv(x, F(x)); // compute Newton correction
    x -= s;             // compute next iterate

    if (callback != nullptr)
      callback(x, s);
  }
  // correction based termination (relative and absolute)
  while ((s.norm() > rtol * x.norm()) && (s.norm() > atol));

  return x;
}
/* SAM_LISTING_END_1 */
#endif
