#ifndef CHEBPOLYPROPERTIES_HPP
#define CHEBPOLYPROPERTIES_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * @brief Computes the values of the Chebychev polynomials
 * \Blue{$T_{0},\ldots,T_{d}$}, \Blue{$d\geq 2$} using the 3-term recursion.
 *
 * @param d maximal order of polynomials
 * @param x evaluation points
 * @param V values \Blue{$T_k(x_j)$} are returned in row \Blue{$k+1$}
 */
/* SAM_LISTING_BEGIN_3 */
void chebpolmult(const unsigned int d, const Eigen::RowVectorXd& x,
                 Eigen::MatrixXd& V) {
  const unsigned int n = x.size();
  V = Eigen::MatrixXd::Ones(d + 1, n);  // \Blue{$T_0 \equiv 1$}
  if (d > 0) {
    V.block(1, 0, 1, n) = x;  // \Blue{$T_1(x) = x$}
    for (unsigned int k = 1; k < d; ++k) {
      const Eigen::RowVectorXd p = V.block(k, 0, 1, n);      // $p = T_{k}$
      const Eigen::RowVectorXd q = V.block(k - 1, 0, 1, n);  // $q = T_{k-1}$
      V.block(k + 1, 0, 1, n) =
          2 * x.cwiseProduct(p) - q;  // \Magenta{3-term recursion}
    }
  }
}
/* SAM_LISTING_END_3 */

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

  // END

  return errors;
}
/* SAM_LISTING_END_2 */

#endif