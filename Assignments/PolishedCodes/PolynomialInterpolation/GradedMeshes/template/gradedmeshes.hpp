#ifndef GRADEDMESHES_HPP
#define GRADEDMESHES_HPP

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * @brief Least squares polynomial fit to the data (x,y).
 * @param x, y vectors of same size
 * @param order, the degree of the fitted polynomial
 * @param coeff, vector of size order+1 containing the monomial
 * coefficients of the fitted polynomial p(x) = coeff[0] + coeff[1]*x + ... +
 * coeff[order]*x^order.
 */
/* SAM_LISTING_BEGIN_9 */
Eigen::VectorXd polyfit(const Eigen::VectorXd &x, const Eigen::VectorXd &y,
                        std::size_t order) {
  Eigen::MatrixXd A(x.size(), order + 1);

  assert(x.size() == y.size());
  assert(x.size() >= order + 1);

  // Create matrix
  for (std::size_t i = 0; i < x.size(); ++i) {
    for (std::size_t j = 0; j < order + 1; ++j) {
      A(i, j) = pow(x(i), j);
    }
  }

  // Solve for linear least squares fit
  Eigen::VectorXd coeff = A.householderQr().solve(y);
  coeff.conservativeResize(order + 1);

  return coeff;
}
/* SAM_LISTING_END_9 */

/**
 * @brief Compute the $L^\infty$ error of piecewise linear interpolation.
 * @param f function with evaluation operator double
 * operator()(double)const.
 * @param t sorted mesh with t[0] = 0 and t[t.size()-1] = 1.
 * @param $\|f - I_{\mathcal{M}}f\|_{L^\infty}$
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FUNCTION>
double pwlintpMaxError(FUNCTION &&f, const Eigen::VectorXd &t) {
  constexpr unsigned int N =
      1e5;              // Sampling resolution for approximating L-infinty norm.
  double maxerr = 0.0;  // L-infinity error so far.

  // TODO: (6-1.c) Approximate the maximum norm of (f - If) on [0,1],
  // where If is the piecewise linear interpolation of f on the mesh t.
  // START

  // END
  return maxerr;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Plots error norms of piecewise linear interpolation of
 * $f(t)=t^\alpha$ using equidistant meshes on [0,1].
 * @param alpha, a vector of positive values not equal to 0 or 1.
 */
/* SAM_LISTING_BEGIN_5 */
void cvgplotEquidistantMesh(const Eigen::VectorXd &alpha) {
  const unsigned int n_alphas = alpha.size();
  plt::figure();
  // TODO: (6-1.d) Create log-log plots of the maximum norm
  // errors (obtained by pwlintMaxError()) with
  // number of mesh intervals 32, 64, ..., 4096.
  // START

  // END
  plt::savefig("./cx_out/cvgplotequidistant.png");
}
/* SAM_LISTING_END_5 */

/**
 * @brief Estimates convergence rates of piecewise linear interpolation of
 * $f(t)=t^\alpha$ using equidistant meshes on [0,1].
 * @param alpha, a vector of positive values not equal to 0 or 1.
 * @param Rates, a vector with Rates[j] = convergence rate of piecewise
 * linear interpolation of t^alpha[j] using equidistant meshes on [0,1].
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd cvgrateEquidistantMesh(const Eigen::VectorXd &alpha) {
  const unsigned int n_alphas = alpha.size();
  Eigen::VectorXd Rates(n_alphas);

  // TODO: (6-1.e) Fill in the entries of Rates.
  // Hint: For each alpha(j), use polyfit() to estimate convergence
  // rates of the maximum norm (obtained by pwlintMaxError()) with
  // number of mesh intervals 32, 64, ..., 4096.
  // START

  // END
  return Rates;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Tabulates and plots convergence rate of piecewise linear interpolation
 * of t^alpha on equidistant meshes.
 */
/* SAM_LISTING_BEGIN_2 */
void testcvgEquidistantMesh() {
  plt::figure();
  // TODO: (6-1.f) Plot and tabulate the convergence rates of pw.
  // lin. intp. of $t^\alpha$ using equidistant meshes on [0,1],
  // for alpha = 0.05, 0.15, 0.25 ..., 2.95.
  // START

  // END
  plt::savefig("./cx_out/cvgRateEquidistant.png");
}
/* SAM_LISTING_END_2 */

/**
 * @brief Estimates convergence rates of piecewise linear interpolation of
 * $f(t)=t^\alpha$ on [0,1].
 * @param alpha, a vector of positive values not equal to 0 or 1.
 * @param beta, a vector of positive values.
 * @param Rates, a matrix with Rates(k,j) = convergence rate of piecewise
 * linear interpolation of t^alpha[j] using a beta[k]-graded mesh on [0,1].
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd cvgrateGradedMesh(const Eigen::VectorXd &alpha,
                                  const Eigen::VectorXd &beta) {
  const unsigned int n_alphas = alpha.size();
  const unsigned int n_betas = beta.size();
  Eigen::MatrixXd Rates(n_betas, n_alphas);

  // TODO: (6-1.i) Fill in the entries of Rates.
  // START

  // END
  return Rates;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Plots convergence rates of piecewise linear interpolation of t^alpha
 * on graded meshes.
 */
/* SAM_LISTING_BEGIN_4 */
void testcvgGradedMesh() {
  plt::figure();
  // TODO: (6-1.j) Plot the convergence rates from
  // cvgrateGradedMesh() using alpha = 0.05, 0.15, ..., 2.95
  // and beta = 0.1, 0.2, ..., 2.0.
  // Note: Running this code may take a while, so start by
  // using fewer values for alpha and beta.
  // Hint: You can use plt::plot_surface(X,Y,Z) where X, Y, Z
  // are matrices that all have the same dimensions.
  // You can use x.replicate() to create X and Y
  // for an appropriate vector x.
  // START

  // END
  plt::savefig("./cx_out/alphabeta.png");
}
/* SAM_LISTING_END_4 */

#endif