#pragma once

#include <Eigen/SVD>
#include <iostream>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
VectorXd symRankOneBestApproxSym(const MatrixXd &M) {
  // ensure that M is symmetric
  assert(M.cols() == M.rows() && M == M.transpose() && "M is not symmetric");
  // Vector to be used to return result
  VectorXd z(M.cols());
  // TODO: (9-14.a)
  //  solve equ. (9.14.1) for a symmetric matrix with the aid of singular value
  //  decomposition (SVD)
  // START

  // END
  return z;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
VectorXd computeKronProdVecMult(const VectorXd &v, const VectorXd &b) {
  const unsigned int n = v.size();
  assert(b.size() == n * n && "Size mismatch");
  // Vector to be used to return result
  VectorXd res;
  // TODO: (9-14.e)
  //  evaluate the expression (9.14.14) efficiently with the aid of Eigen::map
  // START

  // END
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
VectorXd symmRankOneApprox(const MatrixXd &M, const double rtol = 1e-6,
                           const double atol = 1e-8) {
  // ensure M is square
  assert(M.cols() == M.rows() && "Matrix must be square");

  const unsigned int n = M.cols();
  // Initial guess: Column of symmetric part of M with largest norm
  auto Msym = 0.5 * (M + M.transpose());
  Eigen::Index j, dummy;
  Msym.colwise().norm().maxCoeff(&dummy, &j);
  // Set initial guess. This vector is also used to return the result
  VectorXd z{Msym.col(j)};
  // TODO: (9-14.g)
  //  solve equ. (9.14.1) using Gauss-Newton iteration with tolerances rtol and
  //  atol
  // START

  // END
  return z;
}
/* SAM_LISTING_END_2 */
