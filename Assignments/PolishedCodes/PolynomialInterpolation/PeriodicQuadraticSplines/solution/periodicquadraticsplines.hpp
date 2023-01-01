#ifndef PERIODICQUADRATICSPLINESHPP
#define PERIODICQUADRATICSPLINESHPP

/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

/* SAM_LISTING_BEGIN_1 */
class ClosedQuadraticSplineCurve {
 public:
  // Constructor
  explicit ClosedQuadraticSplineCurve(const Eigen::Matrix2Xd &p);
  ~ClosedQuadraticSplineCurve() = default;
  ClosedQuadraticSplineCurve() = delete;
  ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &) = delete;
  ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &&) = delete;
  ClosedQuadraticSplineCurve &operator=(const ClosedQuadraticSplineCurve &) =
      delete;
  // Point evaluation operator for sorted parameter arguments
  [[nodiscard]] Eigen::Matrix2Xd curve_points(const Eigen::VectorXd &v) const;
  // Curvature evaluation operator for sorted parameter arguments
  [[nodiscard]] Eigen::VectorXd local_curvatures(
      const Eigen::VectorXd &v) const;
  // Approximate computation of the length of the curve
  [[nodiscard]] double length(double rtol = 1E-6) const;

 private:
  [[nodiscard]] bool checkC1(double tol = 1.0E-4) const;
  // Number of points to be interpolated
  unsigned int n_;
  // Coordinates of interpolated points, duplicate: $\cob{\Vp^0=\Vp^n}$
  Eigen::Matrix2Xd p_;
  // Knot sequence
  Eigen::VectorXd t_;
  // Coefficients in local representation
  Eigen::Matrix2Xd x_;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
ClosedQuadraticSplineCurve::ClosedQuadraticSplineCurve(
    const Eigen::Matrix2Xd &p)
    : n_(p.cols()), p_(2, p.cols() + 1), t_(p.cols() + 1), x_(2, p.cols()) {
  assert((n_ > 2) && "At least three points have to be supplied");
  assert((n_ % 2 == 1) &&
         "Number of points must be odd!");  // \Label[line]{pqs:1}
  // Save point coordinates, duplicate endpoints
  p_.col(0) = p.col(n_ - 1);
  p_.rightCols(n_) = p;
  // TODO: (5-15.c) Write the constructor.
  // START
  // Compute knot set by accumulating length of polygon segments
  Eigen::VectorXd lensum{n_};
  lensum[0] = (p_.col(1) - p_.col(0)).norm();
  for (unsigned int l = 1; l < n_; ++l) {
    lensum[l] = lensum[l - 1] + (p_.col(l + 1) - p_.col(l)).norm();
  }
  t_[0] = 0.0;
  for (unsigned int l = 1; l <= n_; ++l) {
    t_[l] = lensum[l - 1] / lensum[n_ - 1];
  }
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
// Compute vector coefficients $\cob{\Vx^j}$
#ifndef SMART_LSE_SOLUTION
  // Initialize sparse matrix efficiently, see
  Eigen::SparseMatrix<double> A(n_, n_);
  A.reserve(Eigen::RowVectorXi::Constant(n_, 2));
  A.insert(0, 0) = 1.0;
  A.insert(0, n_ - 1) = 1.0;
  for (unsigned int k = 1; k < n_; ++k) {
    A.insert(k, k) = 1.0;
    A.insert(k, k - 1) = 1.0;
  }
  A.makeCompressed();
  // Compute LU-decomposition
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Matrix factorization failed!");
  }
#endif
  // Right-hand side vector for linear system of equations
  Eigen::VectorXd b(n_);
  // Treat first and second component
  for (unsigned int d = 0; d <= 1; ++d) {
    // Compile right-hand side vector
    const auto c{p_.row(d)};
    b[0] = (c[n_] - c[n_ - 1]) / h[n_ - 1] - (c[1] - c[0]) / h[0];
    for (unsigned int k = 1; k < n_; ++k) {
      b[k] = (c[k] - c[k - 1]) / h[k - 1] - (c[k + 1] - c[k]) / h[k];
    }
#ifdef SMART_LSE_SOLUTION
    // Forward elimination in-situ in b
    for (unsigned int j = 1; j < n_; ++j) {
      b[j] -= b[j - 1];
    }
    // Backward substitution in-situ in b
    b[n_ - 1] *= 0.5;
    int sign = 1;
    for (unsigned int j = 0; j < n_ - 1; ++j) {
      b[j] -= (sign * b[n_ - 1]);
      sign *= -1;
    }
    x_.row(d) = b.transpose();
#else
    // Solve sparse linear system
    x_.row(d) = solver.solve(b).transpose();
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Elimination failed");
    }
#endif
  }
  // END
}
/* SAM_LISTING_END_2 */

bool ClosedQuadraticSplineCurve::checkC1(double tol) const {
  bool ok = true;
  // Lengths of knot intervals
  Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  double diff = 1.0E-8;
  // Last point
  double tl = 1.0 - diff / h[n_ - 1];
  Eigen::Vector2d left_val = (1.0 - tl) * p_.col(n_ - 1) +
                             h[n_ - 1] * x_.col(n_ - 1) * tl * (1.0 - tl) +
                             tl * p_.col(n_);
  double tr = diff / h[0];
  Eigen::Vector2d right_val = (1.0 - tr) * p_.col(0) +
                              h[0] * x_.col(0) * tr * (1.0 - tr) +
                              tr * p_.col(1);
  Eigen::Vector2d left_slope = (left_val - p_.col(0)) / diff;
  Eigen::Vector2d right_slope = (p_.col(n_) - right_val) / diff;
  std::cout << "CHECK: at p[0]=p[n]: slope diff = "
            << (left_slope - right_slope).norm() << std::endl;
  // First point
  tl = 1.0 - diff / h[0];
  tr = diff / h[0];
  left_val = (1.0 - tl) * p_.col(0) + h[0] * x_.col(0) * tl * (1.0 - tl) +
             tl * p_.col(1);
  right_val = (1.0 - tr) * p_.col(1) + h[1] * x_.col(1) * tr * (1.0 - tr) +
              tr * p_.col(2);
  left_slope = (left_val - p_.col(1)) / diff;
  right_slope = (p_.col(1) - right_val) / diff;
  std::cout << "CHECK: at p[0]=p[n]: slope diff = "
            << (left_slope - right_slope).norm() << std::endl;
  // Internal points
  for (unsigned int k = 1; k < n_ - 1; ++k) {
    tl = 1.0 - diff / h[k];
    tr = diff / h[k + 1];
    left_val = (1.0 - tl) * p_.col(k) + h[k] * x_.col(k) * tl * (1.0 - tl) +
               tl * p_.col(k + 1);
    right_val = (1.0 - tr) * p_.col(k + 1) +
                h[k + 1] * x_.col(k + 1) * tr * (1.0 - tr) + tr * p_.col(k + 2);
    left_slope = (left_val - p_.col(k + 1)) / diff;
    right_slope = (p_.col(k + 1) - right_val) / diff;
    std::cout << "CHECK: at p[0]=p[n]: slope diff = "
              << (left_slope - right_slope).norm() << std::endl;
  }
  return ok;
}

// Computation of points on the curve for many parameter values
// Those are assumed to be sorted
/* SAM_LISTING_BEGIN_3 */
Eigen::Matrix2Xd ClosedQuadraticSplineCurve::curve_points(
    const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Matrix containing points to be computed
  Eigen::Matrix2Xd s(2, N);
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  // TODO: (5-15.d) Compute the points on the curve for the sorted parameter
  // values v.
  // START

  // Run through all the provided parameter values
  unsigned int knot_idx = 0;
  for (unsigned int k = 0; k < N; ++k) {
    assert(((v[k] >= 0.0) && (v[k] < 1.0)) && "Out of range!");
    if (k > 0) {
      assert(v[k] >= v[k - 1] && "Parameter values not sorted!");
    }
    // Find right knot interval: knot\_idx stores index of knot to the right of
    // current parameter value
    while ((knot_idx < n_) && (v[k] >= t_[knot_idx])) {
      knot_idx++;
    }
    assert(knot_idx > 0);
    const double tau = (v[k] - t_[knot_idx - 1]) / h[knot_idx - 1];
    s.col(k) = ((1.0 - tau) * p_.col(knot_idx - 1)) +
               (h[knot_idx - 1] * tau * (1 - tau)) * x_.col(knot_idx - 1) +
               (tau * p_.col(knot_idx));
  }
  // END
  return s;
}
/* SAM_LISTING_END_3 */

// Computation of local curvatures of the curve for many parameter values
// Those are assumed to be sorted
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd ClosedQuadraticSplineCurve::local_curvatures(
    const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Vector for storing the computed curvatures
  Eigen::VectorXd c(N);
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  // TODO: (5-15.e) Compute the local curvature at the parameter values v.
  // START

  // Same as in curve\_points(): Run through all the
  // provided parameter values
  unsigned int knot_idx = 0;
  for (unsigned int k = 0; k < N; ++k) {
    assert(((v[k] >= 0.0) && (v[k] < 1.0)) && "Out of range!");
    if (k > 0) {
      assert(v[k] >= v[k - 1] && "Parameter values not sorted!");
    }
    // Find right knot interval: knot\_idx stores index of knot to the right of
    // current parameter value
    while ((v[k] >= t_[knot_idx]) && (knot_idx < n_)) {
      knot_idx++;
    }
    assert(knot_idx > 0);
    const double tau = (v[k] - t_[knot_idx - 1]) / h[knot_idx - 1];
    // Derivative vector
    const Eigen::Vector2d ds =
        (p_.col(knot_idx) - p_.col(knot_idx - 1)) / h[knot_idx - 1] +
        (1.0 - 2 * tau) * x_.col(knot_idx - 1);
    // Second derivative vector
    const Eigen::Vector2d dds = -2.0 * x_.col(knot_idx - 1);
    // Evaluate formula for signed curvature
    const double determinant = ds[0] * dds[1] - ds[1] * dds[0];
    c[k] = determinant / ds.norm();
  }
  // END
  return c;
}
/* SAM_LISTING_END_4 */

// Adaptive computation of the length of the curve by means
// of polygonal approximation up to a given relative tolerance
double ClosedQuadraticSplineCurve::length(double rtol) const {
  double length = 0.0;

  // TODO: (5-15.f) Return an approximation of the length of the curve.
  // START

  // END

  return length;
}

#endif
