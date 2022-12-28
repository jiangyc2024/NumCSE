#ifndef QSPLINES_HPP
#define QSPLINES_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * \brief Computes length of knot intervals
 *
 * \param t Vector of size n-1
 * \return std::pair<Eigen::VectorXd, Eigen::VectorXd> of size n+2 and n+1,
 * respectively
 */
/* SAM_LISTING_BEGIN_0 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> increments(
    const Eigen::VectorXd &t) {
  Eigen::VectorXd dt, ddt;
  // TODO: (5-6.g) compute the increments t_j - t_{j-1}
  // for j = 0,...,n+1 and t_{j+1} - t_{j-1} for j = 0,...,n
  // Assume that t is sorted and that t does not include the endpoints.
  // Hint: use periodic conditions to define t_{-1}, t_{n+1} for vectorization.

  // START
  // Create extended vector using the periodicity:
  const unsigned int n = t.size() + 1;  // number of intervals
  Eigen::VectorXd ext_t(n + 3);
  ext_t << t(n - 2) - 1, 0, t, 1, 1 + t(0);
  // Increments in t
  dt = ext_t.tail(n + 2).array() - ext_t.head(n + 2).array();
  ddt = ext_t.tail(n + 1).array() - ext_t.head(n + 1).array();
  // END

  return std::make_pair(dt, ddt);
}
/* SAM_LISTING_END_0 */

/**
 * \brief Computes spline coefficients c
 *
 * \param t Vector of size n-1
 * \param y Vector of size n
 * \return Eigen::VectorXd Vector of size n
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd compute_c(const Eigen::VectorXd &t, const Eigen::VectorXd &y) {
  const unsigned int n = y.size();  // number of intervals
  Eigen::VectorXd c(n);
  // TODO: (5-6.g) Build the (sparse) matrix for spline interpolation
  // at midpoints. Then compute the coefficients c for the data y.

  // START
  assert(t.size() == n - 1 && "number of intervals mismatch");
  assert(t.minCoeff() > 0 && t.maxCoeff() < 1 && "mesh nodes out of range");

  std::pair<Eigen::VectorXd, Eigen::VectorXd> p = increments(t);
  Eigen::VectorXd dt = p.first;
  Eigen::VectorXd ddt = p.second;

  // Note: A_0 = A_n and C_1 = C_{n+1}
  Eigen::VectorXd A = dt.tail(n).cwiseQuotient(2 * ddt.tail(n));
  Eigen::VectorXd C = dt.head(n).cwiseQuotient(2 * ddt.head(n));

  std::vector<Eigen::Triplet<double>> triplets(3 * n);

  // Build matrix as triplets
  for (unsigned int j = 0; j < n - 1; ++j) {
    triplets.push_back(Eigen::Triplet<double>(j, j, A(j) + C(j) + 1));
    triplets.push_back(Eigen::Triplet<double>(j, j + 1, C(j + 1)));
    triplets.push_back(Eigen::Triplet<double>(j + 1, j, A(j)));
  }
  triplets.push_back(
      Eigen::Triplet<double>(n - 1, n - 1, A(n - 1) + C(n - 1) + 1));
  triplets.push_back(Eigen::Triplet<double>(n - 1, 0, C(0)));
  triplets.push_back(Eigen::Triplet<double>(0, n - 1, A(n - 1)));

  Eigen::SparseMatrix<double> M(n, n);
  M.setFromTriplets(triplets.begin(), triplets.end());
  M.makeCompressed();
  // solve linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  c = solver.solve(y);
  // END
  return c;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Computes spline coefficients d with known c
 *
 * \param c Vector of size n
 * \param t Vector of size n-1
 * \return Eigen::VectorXd Vector of size n+1
 */
/* SAM_LISTING_BEGIN_9 */
Eigen::VectorXd compute_d(const Eigen::VectorXd &c, const Eigen::VectorXd &t) {
  const unsigned int n = c.size();  // number of intervals
  Eigen::VectorXd d_ext(n + 1);
  // TODO: (5-6.h) compute coefficients d_j for j = 0...n
  // Hint: periodic conditions give d_0 = d_n
  // START
  std::pair<Eigen::VectorXd, Eigen::VectorXd> p = increments(t);
  Eigen::VectorXd dt = p.first;
  Eigen::VectorXd ddt = p.second;

  // extend c periodically for vectorization
  Eigen::VectorXd c_ext(c.size() + 1);
  c_ext << c, c(0);

  // coefficients d_j for j = 1...n
  Eigen::VectorXd d =
      c.cwiseProduct(dt.tail(n)) + c_ext.tail(n).cwiseProduct(dt.segment(1, n));
  d = 2 * d.cwiseQuotient(ddt.tail(n));
  // extend d by periodic condition
  d_ext << d(n - 1), d;
  // END
  return d_ext;
}
/* SAM_LISTING_END_9 */

/**
 * \brief Computes values of interpolating quadratic spline at x
 *
 * \param t Vector of size n-1
 * \param y Vector of size n
 * \param x Vector of size N
 * \return Eigen::VectorXd Vector of size N
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd quadspline(const Eigen::VectorXd &t, const Eigen::VectorXd &y,
                           const Eigen::VectorXd &x) {
  Eigen::VectorXd fval(x.size());

  // TODO: (5-6.i) evaluate the spline at the points defined in x.
  // Assume that x is sorted.
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 &&
         "evaluation samples out of range");

  // START
  Eigen::VectorXd c = compute_c(t, y);
  const unsigned int n = y.size();  // number of intervals

  std::pair<Eigen::VectorXd, Eigen::VectorXd> p = increments(t);
  Eigen::VectorXd dt = p.first;
  Eigen::VectorXd ddt = p.second;

  Eigen::VectorXd d_ext = compute_d(c, t);
  // extend t for vectorization
  Eigen::VectorXd t_ext(n + 1);
  t_ext << 0, t, 1;

  unsigned int i = 0;
  // loop over intervals
  for (unsigned int j = 1; j <= n; ++j) {
    // loop over sample points in the interval
    while (i < x.size() && x(i) <= t_ext(j)) {
      const double tau = (x(i) - t_ext(j - 1)) / dt(j);
      const double uat = 1 - tau;
      fval(i) = d_ext(j) * tau * tau + 4 * c(j - 1) * tau * uat +
                d_ext(j - 1) * uat * uat;
      ++i;
    }
    if (i == x.size()) {
      break;
    }
  }
  // END
  return fval;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void plotquadspline(const std::string &filename) {
  Eigen::VectorXd mesh = Eigen::VectorXd::LinSpaced(9, .1, .9);
  auto f = [](Eigen::VectorXd t) { return (2 * M_PI * t).array().sin().exp(); };

  plt::figure();
  // TODO: (5-6.j) plot the quadratic spline for the function f based on
  // the intervals defined in t. Plot also the data (interpolation) points.

  // START
  // The interpolation points are the midpoints of the intervals
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(10, .05, .95);
  Eigen::VectorXd y = f(t);
  // plot data points
  plt::plot(t, y, "o", {{"label", "data"}});

  // std::cout << y <<std::endl;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(200, 0, 1);
  Eigen::VectorXd spline_val = quadspline(mesh, y, x);
  // std::cout << spline_val <<std::endl;

  plt::plot(x, spline_val, {{"label", "spline"}});
  plt::legend();
  plt::title("Quadratic spline");
  plt::xlabel("t");
  plt::ylabel("s(t)");
  // END
  plt::savefig("./cx_out/" + filename + ".png");
}
/* SAM_LISTING_END_3 */

/* [input] q, integer
   [output] Err vector of size n = 2^q
*/
/**
 * \brief Computes error norms for n = 2, 4, ..., $2^q$
 *
 * \param q
 * \return std::vector<double> error normes
 */
/* SAM_LISTING_BEGIN_4 */
std::vector<double> qsp_error(unsigned int q) {
  std::vector<double> Err;
  constexpr unsigned int N = 10000;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, 0, 1);
  auto f = [](Eigen::VectorXd t) { return (2 * M_PI * t).array().sin().exp(); };
  // TODO: (5-6.k) compute L^infty errors for all n = 2, 4, ..., 2^q
  // START
  const unsigned int n_max = std::pow(2, q);
  Eigen::VectorXd f_exact = f(x);

  for (unsigned int n = 2; n <= n_max; n <<= 1) {
    Eigen::VectorXd t =
        Eigen::VectorXd::LinSpaced(n - 1, 1.0 / n, 1.0 - 1.0 / n);
    // The interpolation points are the midpoints of the intervals
    Eigen::VectorXd interp =
        Eigen::VectorXd::LinSpaced(n, .5 / n, 1.0 - .5 / n);
    Eigen::VectorXd y = f(interp);
    Eigen::VectorXd f_spline = quadspline(t, y, x);
    // compute max error
    Err.push_back((f_spline - f_exact).cwiseAbs().maxCoeff());
  }
  // END
  return Err;
}
/* SAM_LISTING_END_4 */

#endif
