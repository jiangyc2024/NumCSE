#pragma once

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
VectorXd symRankOneBestApproxSym(const MatrixXd &M) {
  // ensure that M is symmetric
  assert(M.cols() == M.rows() && M == M.transpose() && "M is not symmetric");
  // Vector to be used to return result
  VectorXd z(M.cols());
  // To do: (0-3.a)
  // START

  // END
  return z;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
VectorXd computeKronProdVecMult(const VectorXd &v, const VectorXd &b) {
  unsigned int n = v.size();
  assert(b.size() == n * n && "Size mismatch");
  // Vector to be used to return result
  VectorXd res;
  // To do: (0-3.d)
  // START

  // END
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
VectorXd symmRankOneApprox(const MatrixXd &M, double rtol = 1e-6,
                           double atol = 1e-8) {
  // ensure M is square
  assert(M.cols() == M.rows() && "Matrix must be square");

  const unsigned int n = M.cols();
  // Initial guess: Column of symmetric part of M with largest norm
  auto Msym = 0.5 * (M + M.transpose());
  Eigen::Index j, dummy;
  Msym.colwise().norm().maxCoeff(&dummy, &j);
  // Set initial guess. This vector is also used to return the result
  VectorXd z{Msym.col(j)};
  // To do: (0-3.e)
  // START

  // END
  return z;
}
/* SAM_LISTING_END_2 */
