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
   * @brief Construct a new PwLinIP object
   *
   * @param x Vector of knots
   * @param t Vector of nodes
   * @param y Vector of values of interpolant in nodes
   */
  PwLinIP(const Eigen::VectorXd& x, const Eigen::VectorXd& t,
          const Eigen::VectorXd& y);

  /**
   * @brief evaluate interpolant at $arg$.
   *
   * @param arg argument to evaluate at
   * @return double
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

  // END
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
double PwLinIP::operator()(double arg) const {
  double ret_val = 0.;

  // TODO: (6-2.e) Implement the evaluation operator efficiently!
  // START

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

  // END

  plt::savefig("cx_out/tent_basis_functions.png");
}
/* SAM_LISTING_END_3 */

#endif