///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): J. Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#ifndef INTPOLYVAL_HPP
#define INTPOLYVAL_HPP

#include <Eigen/Dense>

using namespace Eigen;
// IN:  \texttt{t}: vector of nodes \Blue{$t_0, \ldots, t_n$}
//      \texttt{y}: vector of data \Blue{$y_0, \ldots, y_n$}
//      \texttt{x}: vector of evaluation points \Blue{$x_1, \ldots, x_N$}
// OUT: \texttt{p}: interpolant evaluated at x
/* SAM_LISTING_BEGIN_0 */
VectorXd intpolyval(const VectorXd& t, const VectorXd& y, const VectorXd& x) {
  const unsigned
      n = t.size(),  // no. of interpolation nodes = deg. of polynomial $-1$
      N = x.size();  // no. of evaluation points

  VectorXd p = VectorXd(N);

  // Precompute the weights \Blue{$\lambda_i$} with effort \Blue{$O(n^2)$}
  VectorXd lambda(n);
  for (unsigned k = 0; k < n; ++k) {
    // little workaround: cannot subtract a vector from a scalar
    // -> multiply scalar by vector of ones
    lambda(k) =
        1. / ((t(k) * VectorXd::Ones(k) - t.head(k)).prod() *
              (t(k) * VectorXd::Ones(n - k - 1) - t.tail(n - k - 1)).prod());
  }
  // Compute quotient of weighted sums  of \Blue{$\frac{\lambda_i}{t - t_i}$},
  // effort \Blue{$O(n)$}
  for (unsigned i = 0; i < N; ++i) {
    VectorXd z = (x(i) * VectorXd::Ones(n) - t);

    // check if we want to evaluate at a node <-> avoid division by zero
    double* ptr = std::find(z.data(), z.data() + n, 0.0);
    if (ptr != z.data() + n) {  // if ptr = z.data + n = z.end no zero was found
      p(i) = y(ptr - z.data());  // ptr - z.data gives the position of the zero
    } else {
      VectorXd mu = lambda.cwiseQuotient(z);
      p(i) = (mu.cwiseProduct(y)).sum() / mu.sum();
    }
  }
  return p;
}
/* SAM_LISTING_END_0 */

#endif
