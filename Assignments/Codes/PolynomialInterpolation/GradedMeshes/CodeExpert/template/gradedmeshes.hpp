#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "matplotlibcpp.h"


using namespace Eigen;
namespace plt = matplotlibcpp;

/* @brief Least squares polynomial fit to the data (x,y).
 * @param[in] x, y vectors of same size
 * @param[in] order, the degree of the fitted polynomial
 * @param[out] coeff, vector of size order+1 containing the monomial
 * coefficients of the fitted polynomial p(x) = coeff[0] + coeff[1]*x + ... +
 * coeff[order]*x^order.
 */
/* SAM_LISTING_BEGIN_9 */
VectorXd polyfit(const VectorXd &x, const VectorXd &y, size_t order) {
  Eigen::MatrixXd A(x.size(), order + 1);

  assert(x.size() == y.size());
  assert(x.size() >= order + 1);

  // Create matrix
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < order + 1; ++j) {
      A(i, j) = pow(x(i), j);
    }
  }

  // Solve for linear least squares fit
  VectorXd coeff = A.householderQr().solve(y);
  coeff.conservativeResize(order + 1);

  return coeff;
}
/* SAM_LISTING_END_9 */

/* @brief Compute the $L^\infty$ error of piecewise linear interpolation.
 * @param[in] f function with evaluation operator double
 * operator()(double)const.
 * @param[in] t sorted mesh with t[0] = 0 and t[t.size()-1] = 1.
 * @param[out] $\|f - I_{\mathcal{M}}f\|_{L^\infty}$
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FUNCTION>
double pwlintpMaxError(FUNCTION &&f, const Eigen::VectorXd &t) {
  int n = t.size() - 1;
  int N = 1e5;         // Sampling resolution for approximating L-infinty norm.
  double maxerr = 0.0; // L-infinity error so far.

  // TO DO (7-1.c): Approximate the maximum norm of (f - If) on [0,1],
  // where If is the piecewise linear interpolation of f on the mesh t.
  // START
  
  // END
  return maxerr;
}
/* SAM_LISTING_END_0 */


/* @brief Plots error norms of piecewise linear interpolation of
 * $f(t)=t^\alpha$ using equidistant meshes on [0,1].
 * @param[in] alpha, a vector of positive values not equal to 0 or 1.
 */
/* SAM_LISTING_BEGIN_5 */
void cvgplotEquidistantMesh(const Eigen::VectorXd &alpha) {
  int n_alphas = alpha.size();
  plt::figure();
  // TO DO (7-1.d): Create log-log plots of the maximum norm
  // errors (obtained by pwlintMaxError()) with
  // number of mesh intervals 32, 64, ..., 4096.
  // START
  
  // END
  plt::savefig("./cx_out/cvgplotEquidistant.png");
}
/* SAM_LISTING_END_5 */

/* @brief Estimates convergence rates of piecewise linear interpolation of
 * $f(t)=t^\alpha$ using equidistant meshes on [0,1].
 * @param[in] alpha, a vector of positive values not equal to 0 or 1.
 * @param[out] Rates, a vector with Rates[j] = convergence rate of piecewise
 * linear interpolation of t^alpha[j] using equidistant meshes on [0,1].
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd cvgrateEquidistantMesh(const Eigen::VectorXd &alpha) {
  int n_alphas = alpha.size();
  VectorXd Rates(n_alphas);

  // TO DO (7-1.e): Fill in the entries of Rates.
  // Hint: For each alpha(j), use polyfit() to estimate convergence
  // rates of the maximum norm (obtained by pwlintMaxError()) with
  // number of mesh intervals 32, 64, ..., 4096.
  // START
  
  // END
  return Rates;
}
/* SAM_LISTING_END_1 */

/* @brief Tabulates and plots convergence rate of piecewise linear interpolation
 * of t^alpha on equidistant meshes.
 */
/* SAM_LISTING_BEGIN_2 */
void testcvgEquidistantMesh(void) {
  // TO DO (7-1.f): Plot and tabulate the convergence rates of pw.
  // lin. intp. of $t^\alpha$ using equidistant meshes on [0,1],
  // for alpha = 0.05, 0.15, 0.25 ..., 2.95.
  // START
  
  // END
}
/* SAM_LISTING_END_2 */

/* @brief Estimates convergence rates of piecewise linear interpolation of
 * $f(t)=t^\alpha$ on [0,1].
 * @param[in] alpha, a vector of positive values not equal to 0 or 1.
 * @param[in] beta, a vector of positive values.
 * @param[out] Rates, a matrix with Rates(k,j) = convergence rate of piecewise
 * linear interpolation of t^alpha[j] using a beta[k]-graded mesh on [0,1].
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd cvgrateGradedMesh(const Eigen::VectorXd &alpha, const Eigen::VectorXd &beta) {
  int n_alphas = alpha.size();
  int n_betas = beta.size();
  MatrixXd Rates(n_betas, n_alphas);

  // TO DO (7-1.i): Fill in the entries of Rates.
  // START
  
  // END
  return Rates;
}
/* SAM_LISTING_END_3 */

/* @brief Plots convergence rates of piecewise linear interpolation of t^alpha
 * on graded meshes.
 */
/* SAM_LISTING_BEGIN_4 */
void testcvgGradedMesh(void) {
  // TO DO (7-1.j): Plot the convergence rates from
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
}
/* SAM_LISTING_END_4 */
