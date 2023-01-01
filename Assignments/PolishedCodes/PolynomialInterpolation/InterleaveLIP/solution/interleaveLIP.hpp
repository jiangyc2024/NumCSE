#ifndef INTERLEAVELIP_HPP
#define INTERLEAVELIP_HPP

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <string>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * @brief "argsort": find indices sorting a vector
 *
 * @param values Array for which we want to find the sorting indices
 * @return std::vector<std::size_t> Permutation of indices resulting in sorting
 * of values
 */
std::vector<std::size_t> order(const Eigen::VectorXd& values) {
  std::vector<std::size_t> indices(values.size());
  std::iota(begin(indices), end(indices), static_cast<std::size_t>(0));
  std::sort(begin(indices), end(indices),
            [&](size_t a, size_t b) { return values[a] < values[b]; });
  return indices;
}

/**
 * @brief Interpolator class
 *
 */
/* SAM_LISTING_BEGIN_1 */
class PwLinIP {
 public:
  /**
   * \brief Construct a new PwLinIP object
   *
   * \param x Vector of knots
   * \param t Vector of nodes
   * \param y Vector of values of interpolant in nodes
   */
  PwLinIP(const Eigen::VectorXd& x, const Eigen::VectorXd& t,
          const Eigen::VectorXd& y);

  /**
   * \brief evaluate interpolant at $arg$.
   *
   * \param arg argument to evaluate at
   * \return double
   */
  double operator()(double arg) const;

 private:
  Eigen::VectorXd x_;
  Eigen::VectorXd t_;
  Eigen::VectorXd y_;
  Eigen::VectorXd s_;
};
/* SAM_LISTING_END_1 */

/**
 * @brief Compute values of interpolant in knots $\mathbf{x}$ from $(t_i,y_i)$
 *
 * @param x Vector of knots
 * @param t Vector of nodes
 * @param y Vector of values of interpolant in nodes $\Vt$
 * @return Eigen::VectorXd Vector of values of interpolant in knots $\Vx$
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd tentBasCoeff(const Eigen::VectorXd& x, const Eigen::VectorXd& t,
                             const Eigen::VectorXd& y) {
  // Initialization
  const std::size_t n = t.size();
  Eigen::VectorXd s = Eigen::VectorXd::Zero(n);

  // TODO: (6-2.d) Return the interpolant evaluated at the knots x.
  // START

  auto x_indices = order(x);
  auto t_indices = order(t);
  // You can also implement a solution which does not need
  // sorted vectors and e.g.\ for each knot $x_j$ looks
  // for the closest node $t_{i1}$ and the next closest node $t_{i2}$.
  // However, such solution will not become more efficient
  // if you give as input already sorted vectors: for each knot $x_j$
  // you will always have to iterate along the sorted vector $t$
  // to find the included node $t_i$.

  // Check condition of subproblem c
  std::size_t i = 0;
  std::size_t k = 0;
  for (std::size_t j = 0; j < (n - 1); ++j) {
    bool intervalOK = false;
    while (i < n) {
      const bool inInterval = (x(x_indices[j]) < t(t_indices[i])) &&
                              (t(t_indices[i]) < x(x_indices[j + 1]));

      if (inInterval) {
        intervalOK = true;
        if (i == j) {  // Index of interval which contains 2 nodes
                       // $t_i$ and $t_{i+1}$. After that, we have
                       // $i > j$...
          k = j;
        }
        break;
      } else {
        ++i;
      }
    }
    if (!intervalOK) {
      std::exit(EXIT_FAILURE);
    }
  }

  // 1. Find slope $\gamma$ and intercept $\beta$
  // in interval with 2 nodes $k$
  // 2. Find $s_k$ and $s_{k+1}$
  double gamma = (y(t_indices[k + 1]) - y(t_indices[k])) /
                 (t(t_indices[k + 1]) - t(t_indices[k]));
  double beta = y(t_indices[k]) - gamma * t(t_indices[k]);

  s(x_indices[k]) = gamma * x(x_indices[k]) + beta;
  s(x_indices[k + 1]) = gamma * x(x_indices[k + 1]) + beta;

  // Find intercept, slope and value $s$ at the lower bound $x$
  // for all intervals on the left of interval with 2 nodes $k$
  for (int j = k - 1; j >= 0; --j) {
    gamma = (s(x_indices[j + 1]) - y(t_indices[j])) /
            (x(x_indices[j + 1]) - t(t_indices[j]));
    beta = y(t_indices[j]) - gamma * t(t_indices[j]);

    s(x_indices[j]) = gamma * x(x_indices[j]) + beta;
  }

  // Find intercept, slope and value $s$ at the upper bound $x$
  // for all intervals on the right of interval with 2 nodes $k$
  for (int j = k + 2; j < n; ++j) {
    gamma = (y(t_indices[j]) - s(x_indices[j - 1])) /
            (t(t_indices[j]) - x(x_indices[j - 1]));
    beta = s(x_indices[j - 1]) - gamma * x(x_indices[j - 1]);

    s(x_indices[j]) = gamma * x(x_indices[j]) + beta;
  }
  // END

  return s;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_2 */
PwLinIP::PwLinIP(const Eigen::VectorXd& x, const Eigen::VectorXd& t,
                 const Eigen::VectorXd& y) {
  assert(t.size() == y.size() && t.size() == x.size() &&
         "x, t, y must have same size!");

  const std::size_t n = t.size();
  x_.resize(n);
  t_.resize(n);
  y_.resize(n);

  // TODO: (6-2.e) Fill the member variables of the class.
  // START
  auto x_indices = order(x);
  for (std::size_t i = 0; i < n; ++i) {
    x_(i) = x[x_indices[i]];
  }

  auto t_indices = order(t);
  for (std::size_t i = 0; i < n; ++i) {
    t_(i) = t[t_indices[i]];
    y_(i) = y[t_indices[i]];
  }

  s_ = tentBasCoeff(x_, t_, y_);
  // END
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
double PwLinIP::operator()(double arg) const {
  double ret_val = 0.;

  // TODO: (6-2.e) Implement the evaluation operator efficiently!
  // START
  if (arg < x_(0) || arg > x_(x_.size() - 1)) {
    ret_val = 0;
  } else {
    std::vector<double> x(
        x_.size() -
        1);  // Exclude $x_0$ from search of interval including $arg$
    Eigen::VectorXd::Map(&x.front(), x_.size() - 1) = x_.tail(x_.size() - 1);
    const std::size_t j = std::lower_bound(x.begin(), x.end(), arg) -
                          x.begin() + 1;  // Binary search
    // '+1' is needed to restore indexing from $x_0$

    const double gamma = (s_(j) - s_(j - 1)) / (x_(j) - x_(j - 1));
    const double beta = s_(j - 1) - gamma * x_(j - 1);

    ret_val = gamma * arg + beta;
  }
  // END

  return ret_val;
}
/* SAM_LISTING_END_4 */

/**
 * @brief Plots the cardinal basis functions as described in the problem
 * statement.
 *
 */
/* SAM_LISTING_BEGIN_3 */
void plotCardinalBasisFunctions() {
  // Initialization
  constexpr std::size_t n = 11;

  // Nodes
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0, 10);
  Eigen::VectorXd t(n);
  t(0) = 0;
  t.tail(n - 1).setLinSpaced(n - 1, 0.5, 9.5);

  plt::figure();
  // Plot colors
  std::vector<std::string> C = {"b", "r", "m"};

  // TODO: (6-2.f) Plot cardinal basis functions for interpolation.
  // START

  // Loop over basis functions
  for (unsigned int i = 0; i < 3; ++i) {
    // Basis function $e_i$
    Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
    y(i) = 1;

    // Create interpoland with values $y$
    PwLinIP cardinalBasis(x, t, y);

    // Plotting values
    constexpr unsigned int N = 1000;
    Eigen::VectorXd xval = Eigen::VectorXd::LinSpaced(N, 0, 10);
    Eigen::VectorXd s(N);
    for (unsigned int j = 0; j < N; ++j) {
      s(j) = cardinalBasis(xval(j));
    }

    plt::plot(xval, s, C[i], {{"label", "Basis k = " + std::to_string(i)}});
  }

  plt::xlabel("t");
  plt::ylabel("y");
  plt::legend();
  plt::title("Tent basis functions");
  // END

  plt::savefig("cx_out/tent_basis_functions.png");
}
/* SAM_LISTING_END_3 */

#endif