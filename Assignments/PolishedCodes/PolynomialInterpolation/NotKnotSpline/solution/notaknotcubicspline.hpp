#ifndef NOTAKNOTCUBICSPLINE_HPP
#define NOTAKNOTCUBICSPLINE_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

class NotAKnotCubicSpline {
public:
  NotAKnotCubicSpline(Eigen::VectorXd t, Eigen::VectorXd y);
  ~NotAKnotCubicSpline() = default;

  [[nodiscard]] double eval(double t) const;
  [[nodiscard]] double evalDerivative(double t) const;

private:
  [[nodiscard]] bool checkC2() const;
  [[nodiscard]] int getPtIdx(double t) const;
  Eigen::VectorXd t_; // knot vector
  Eigen::VectorXd y_; // data values
  Eigen::VectorXd c_; // slopes
};

// Compute the slopes of the natural cubic spline interpolant
// The knot vector t is supposed to be sorted!
/* SAM_LISTING_BEGIN_1 */
NotAKnotCubicSpline::NotAKnotCubicSpline(Eigen::VectorXd t, Eigen::VectorXd y)
    : t_(std::move(t)), y_(std::move(y)), c_(y_.size()) {
  const long n = t_.size() - 1;
  assert(n > 2);
  assert(n == y_.size() - 1);
  // Sequence $\cob{\left(h_i\right)_{i=1}^{n}}$, knots have to be sorted
  auto h = (t_.tail(n) - t_.head(n)).array();
  auto hs = h.pow(2);                             // $\cob{h_i^2}$
  auto b = 1.0 / h;                               // $\cob{b_i := h_i^{-1}}$
  auto a = 2.0 * (b.head(n - 1) + b.tail(n - 1)); // $\cob{a_i}$
  // Initialize sparse coefficient matrix for the linear system of equations
  // needed to compute the unknown slopes. Compared to \lref{eq:NCSILSE} only
  // the first and last row has to be adjusted, see \prbeqref{eq:naklse}.
  Eigen::SparseMatrix<double> M(n + 1, n + 1);
  M.reserve(Eigen::RowVectorXi::Constant(n + 1, 3));
  // Matrix row arising from the continuity of
  M.insert(0, 0) = 6.0 / hs[0];
  M.insert(0, 1) = (6.0 / hs[0] - 6.0 / hs[1]);
  M.insert(0, 2) = -6.0 / hs[1];
  for (int i = 1; i < n; ++i) {
    M.insert(i, i - 1) = b[i - 1];
    M.insert(i, i) = a[i - 1];
    M.insert(i, i + 1) = b[i];
  }
  // Matrix due to continuity of $\cob{s^{(3)}}$ at  $\cob{t_{n-1}}$
  M.insert(n, n - 2) = 6.0 / hs[n - 2];
  M.insert(n, n - 1) = 6.0 / hs[n - 2] - 6.0 / hs[n - 1];
  M.insert(n, n) = -6.0 / hs[n - 1];
  // Initialize right-hand-side vector, see \prbeqref{eq:naklse}
  Eigen::VectorXd rhs(n + 1);
  rhs[0] = -12.0 / (hs[0] * h[0]) * (y_[0] - y_[1]) +
           12.0 / (hs[1] * h[1]) * (y_[1] - y_[2]);
  for (int i = 1; i < n; ++i) {
    rhs[i] =
        3.0 * ((y_[i] - y_[i - 1]) / hs[i - 1] + (y_[i + 1] - y_[i]) / hs[i]);
  }
  rhs[n] = -12.0 / (hs[n - 2] * h[n - 2]) * (y_[n - 2] - y_[n - 1]) +
           12.0 / (hs[n - 1] * h[n - 1]) * (y_[n - 1] - y_[n]);
  // Compute slopes by solving an (almost) tridiagonal linear system
  // with asymptotic effort $\cob{O(n)}$ for $\cob{n\to\infty}$.
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("LU factorization failed");
  }
  c_ = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Elimination failed");
  }
  if (checkC2()) {
    std::cout << "C2 continuity satisfied" << std::endl;
  }
}
/* SAM_LISTING_END_1 */

bool NotAKnotCubicSpline::checkC2() const {
  const long n = t_.size() - 1;
  // Vector of h_i's , t has to be sorted
  auto h = (t_.tail(n) - t_.head(n)).array();
  // Vector h_i^2
  auto hs = h.pow(2);

  double spp; // second derivative
  // Visit all knot intervals
  bool error = false;
  for (int j = 1; j <= n; ++j) {
    const double spp_left = 6.0 / hs[j - 1] * (y_[j] - y_[j - 1]) -
                            2.0 / h[j - 1] * (2 * c_[j - 1] + c_[j]);
    const double spp_right = 6.0 / hs[j - 1] * (y_[j - 1] - y_[j]) +
                             2.0 / h[j - 1] * (c_[j - 1] + 2 * c_[j]);
    if (j > 1) {
      if (std::abs(spp - spp_left) >
          1.0E-10 * std::max(std::abs(spp_left), std::abs(spp_right))) {
        std::cout << "Mismatch of 2nd derivative at knot " << j - 1
                  << std::endl;
        error = true;
      }
    }
    spp = spp_right;
  }
  return !error;
}

// Determines the index i of the knot interval [t_{i-1},t_{i}] containing the
// point t

int NotAKnotCubicSpline::getPtIdx(double t) const {
  if ((t < t_[0]) || (t > t_[t_.size() - 1])) {
    throw std::runtime_error("t out of range");
  }
  // Binary search for first knot to the right of t
  const auto *const last = t_.data() + t_.size();
  const auto *const ptr = std::lower_bound(t_.data(), last, t);
  assert(ptr != last);
  const int idx = std::distance(t_.data(), ptr);
  assert(idx > 0);
  return idx;
}

/* SAM_LISTING_BEGIN_2 */
double NotAKnotCubicSpline::eval(double t) const {
  const int idx = getPtIdx(t);
  const double h = t_[idx] - t_[idx - 1];
  const double tau = (t - t_[idx - 1]) / h;
  const double a0 = y_[idx - 1];
  const double a1 = h * c_[idx - 1];
  const double a2 = -3.0 * y_[idx - 1] + 3.0 * y_[idx] - 2.0 * a1 - h * c_[idx];
  const double a3 = 2.0 * y_[idx - 1] - 2.0 * y_[idx] + a1 + h * c_[idx];
  return (a0 + tau * (a1 + tau * (a2 + tau * a3)));
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
double NotAKnotCubicSpline::evalDerivative(double t) const {
  const int idx = getPtIdx(t);
  const double h = t_[idx] - t_[idx - 1];
  const double tau = (t - t_[idx - 1]) / h;
  const double a0 = c_[idx - 1];
  const double a1 =
      6.0 / h * (y_[idx] - y_[idx - 1]) - 4.0 * c_[idx - 1] - 2.0 * c_[idx];
  const double a2 =
      6.0 / h * (y_[idx - 1] - y_[idx]) + 3.0 * c_[idx - 1] + 3.0 * c_[idx];
  return a0 + tau * (a1 + tau * a2);
}
/* SAM_LISTING_END_3 */

#endif
