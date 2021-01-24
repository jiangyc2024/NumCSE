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
#include <functional>

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
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTION, typename JACOBIAN, typename SCALAR,
          int N = Eigen::Dynamic,
          typename MONITOR =
              std::function<void(unsigned int, Vector<SCALAR, N>,
                                 Vector<SCALAR, N>, Vector<SCALAR, N>)>>
Vector<SCALAR, N> broyd(
    FUNCTION &&F, Vector<SCALAR, N> x, JACOBIAN J, SCALAR reltol,
    SCALAR abstol, unsigned int maxit = 20,
    MONITOR &&monitor = [](unsigned int /*itnum*/,
                           const Vector<SCALAR, N> & /*x*/,
                           const Vector<SCALAR, N> & /*fx*/,
                           const Vector<SCALAR, N> & /*dx*/) {}) {
  // Compute first quasi-Newton update
  Vector<SCALAR, N> s = -J.lu().solve(F(x));
  // Compute $\Vx^{(1)}$
  x += s;
  auto f = F(x);
  monitor(0, x, f, s);
  for (unsigned int k = 1;
       ((s.norm() >= reltol * x.norm()) && (s.norm() >= abstol) && (k < maxit));
       ++k) {
    // Rank-1 update of Jacobian \eqref{eq:broydenqn}
    J += f * s.transpose() / s.squaredNorm();
    // Next quasi-Newton correction
    s = -J.lu().solve(f);
    // Update of iterate
    x += s;
    // Evaluate F at next iterate
    f = F(x);
    // Optional output
    monitor(k, x, f, s);
  }
  return x;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTION, typename JACOBIAN, typename SCALAR,
          int N = Eigen::Dynamic,
          typename MONITOR =
              std::function<void(unsigned int, Vector<SCALAR, N>,
                                 Vector<SCALAR, N>, Vector<SCALAR, N>)>>
Vector<SCALAR, N> upbroyd(
    FUNCTION &&F, Vector<SCALAR, N> x, JACOBIAN &&J, SCALAR reltol,
    SCALAR abstol, unsigned int maxit = 20,
    MONITOR &&monitor = [](unsigned int /*itnum*/,
                           const Vector<SCALAR, N> & /*x*/,
                           const Vector<SCALAR, N> & /*fx*/,
                           const Vector<SCALAR, N> & /*dx*/) {}) {
  // Calculate LU factorization of initial Jacobian once, cf. \cref{rem:seqsolvelse}
  auto fac = J.lu();
  // First quasi-Newton correction $\cob{\Delta\Vx^{(0)} := -\VJ_0^{-1}F(\Vx^{(0)})}$
  Vector<SCALAR, N> s = -fac.solve(F(x));
  // Store the first quasi-Newton correction $\cob{\Delta\Vx^{(0)}}$
  std::vector<Vector<SCALAR, N>> dx{s};
  x += s; // $\cob{\Vx^{(1)} := \Vx^{(0)} + \Delta\Vx^{(0)}}$
  auto f = F(x); // Here $\cob{ = F(\Vx^{(1)})}$
  // Array storing simplified quasi-Newton corrections $\cob{\Delta\overline{\Vx}^{(\ell)}}$
  std::vector<Vector<SCALAR, N>> dxs{};
  // Array of denominators $\cob{\N{\Vx^{(\ell)}}_2^{2}+(\Delta\Vx^{(\ell)})^{\top}{\Delta\overline{\Vx}^{(\ell+1)}}}$
  std::vector<SCALAR> den{};
  monitor(0, x, f, s); // Record start of iteration 
  // Main loop with correction based termination control
  for (unsigned int k = 1;
       ((s.norm() >= reltol * x.norm()) && (s.norm() >= abstol) && (k < maxit));
       ++k) {
    // Compute $\cob{\VJ_{0}^{-1}F(\Vx^{(k)})}$, needed for both recursions
    s = fac.solve(f);
    // \eqref{eq:qncrec}: recursion for next simplified quasi-Newton correction 
    Vector<SCALAR, N> ss = s;
    for (unsigned int l = 1; l < k; ++l) {
      ss -= dxs[l - 1] * (dx[l - 1].dot(ss)) / den[l - 1];
    }
    // Store next denominator $\cob{\N{\Vx^{(k-1)}}_{2}^{2}+(\Delta\Vx^{(k-1)})^{\top}{\Delta\overline{\Vx}^{(k)}}}$
    den.push_back(dx[k - 1].squaredNorm() + dx[k - 1].dot(ss));
    // Store current simplified quasi-Newton correction $\cob{\Delta\overline{\Vx}^{(k)}}$
    dxs.push_back(ss);
    // \eqref{eq:broyrec}: Compute next quasi-Newton correction recursively
    for (unsigned int l = 0; l < k; ++l) {
      s -= dxs[l] * (dx[l].dot(s)) / den[l];
    }
    s *= (-1.0); // Comply with sign convention
    dx.push_back(s);
    // Compute next iterate $\cob{\Vx^{(k+1)}}$ and $\cob{F(\Vx^{(k+1)})}$
    x += s;
    f = F(x);
    monitor(k, x, f, s); // Record progress
  }
  return x;
}
/* SAM_LISTING_END_1 */
#endif
