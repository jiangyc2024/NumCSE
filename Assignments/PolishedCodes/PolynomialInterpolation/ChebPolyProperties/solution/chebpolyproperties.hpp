#ifndef CHEBPOLYPROPERTIES_HPP
#define CHEBPOLYPROPERTIES_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "chebpolmult.hpp"

/**
 * @brief Check the orthogonality of Chebychev polynomials
 *
 * @param n maximal order to check
 * @return true if polynomials are orthogonal
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_1 */
bool checkDiscreteOrthogonality(unsigned int n) {
  bool return_val = false;
  constexpr double tol =
      std::numeric_limits<double>::epsilon();  // machine epsilon

  // TODO: (6-4.d) Check the orthogonality property of the Chebychev polynomials
  // up to order n.
  // START
  Eigen::MatrixXd V;
  // generate Chebychev nodes
  Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(n + 1, 0, n);
  x = x.unaryExpr(
      [n](double x) { return std::cos(M_PI * (2. * x + 1.) / 2. / (n + 1.)); });
  chebpolmult(n, x, V);  // compute values of polynomials

  // check for orthogonality
  // note that the second loop starts at k + 1 because this inner product is
  // (anti-)symmetric
  return_val = true;
  for (unsigned int k = 0; k < n + 1; ++k) {
    for (unsigned int l = k + 1; l < n + 1; ++l) {
      return_val &= (std::abs(V.row(k).dot(V.row(l))) < tol * 100.);
    }
  }
  // END

  return return_val;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Calculates the vector of coefficients $\alpha_j, \quad j = 0,\dots,m$,
 * given f
 *
 * @tparam Function of signature std::function<double(double)>
 * @param f the function
 * @param n number of Chebychev nodes
 * @param m degree of polynomial
 * @return Eigen::VectorXd vector of coefficients
 */
/* SAM_LISTING_BEGIN_0 */
template <typename Function>
Eigen::VectorXd bestpolchebnodes(const Function& f, unsigned int n,
                                 unsigned int m) {
  Eigen::VectorXd alpha = Eigen::VectorXd::Zero(m + 1);

  // TODO: (6-4.h) Compute the vector of coefficients such that the Chebychev
  // polynomials of up to order m interpolate the function f.
  // START

  // generate Chebychev nodes
  Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(n + 1, 0, n);
  x = x.unaryExpr(
      [n](double x) { return std::cos(M_PI * (2. * x + 1.) / 2. / (n + 1.)); });

  // evaluate function
  Eigen::VectorXd f_eval = x.unaryExpr(f);

  // evaluate Chebychev polynomials
  Eigen::MatrixXd V;
  chebpolmult(m, x, V);

  // find coefficients using previously determined formula in subproblem (f)
  alpha = 2. * V * f_eval / (n + 1.);
  alpha(0) /= 2.;
  // END

  return alpha;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Tests bestpolchebnodes and tabulates errors.
 *
 * @param n maximal degree of Chebychev polynomials
 * @return std::vector<double> of errors
 */
/* SAM_LISTING_BEGIN_2 */
std::vector<double> testBestPolyChebNodes(unsigned int n) {
  std::vector<double> errors(n + 1, 0.);
  auto f = [](double x) {
    return 1. / (std::pow(5. * x, 2) + 1.);
  };  // function to test on
  const Eigen::VectorXd X =
      Eigen::VectorXd::LinSpaced(1e6, -1, 1);  // evaluation grid

  // TODO: (6-4.i) Test the function from the previous subproblem by tabulating
  // the error norms.
  // START
  std::cout << " m  L^inf error" << std::endl;

  for (unsigned int m = 0; m < n + 1; ++m) {
    Eigen::VectorXd alpha = bestpolchebnodes(f, m, m);

    Eigen::MatrixXd V;
    chebpolmult(m, X, V);

    Eigen::VectorXd qm = V.transpose() * alpha;
    Eigen::VectorXd f_eval = X.unaryExpr(f);

    double err_max = (f_eval - qm).lpNorm<Eigen::Infinity>();

    // tabulate errors
    std::cout << std::setw(2) << m << std::setw(2) << " "
              << std::setprecision(10) << err_max << std::endl;
  }
  // END

  return errors;
}
/* SAM_LISTING_END_2 */

#endif