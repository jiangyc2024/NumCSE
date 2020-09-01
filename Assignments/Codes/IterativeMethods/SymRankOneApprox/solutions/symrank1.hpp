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
  // Compute economical SVD of M
  Eigen::JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
  // Extract first column of u and first singular value
  VectorXd first_U = svd.matrixU().col(0);
  VectorXd sv = svd.singularValues();
  // Scaling
  z = first_U * std::sqrt(sv[0]);
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
  // Raw data array for argument vector b
  const double *data = b.data();
  // Twice reinterpretation of data, no copy
  res = (MatrixXd::Map(data, n, n) + MatrixXd::Map(data, n, n).transpose()) * v;
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
  // Gauss-Newton correction vector
  VectorXd s(n);
  // Direct access to the array storing the entries of M
  const double *dataM = M.data();

  // Main Gauss-Newton loop with correction-based termination
  do {
    // Compute b = vec(M) - kron(z,z)
    const MatrixXd Z = z * z.transpose();
    const double *dataz = Z.data();
    const VectorXd b =
        MatrixXd::Map(dataM, n * n, 1) - MatrixXd::Map(dataz, n * n, 1);
    // Solve normal equations
    double normalize = 2 * z.squaredNorm();
    VectorXd kronMult = computeKronProdVecMult(z, b);
    s = (kronMult - z.dot(kronMult) / normalize * z) / normalize;
    // Apply Gauss-Newton update
    z = z + s;
    std::cout << "|s| = " << s.norm() << std::endl;
  } while (s.norm() > std::max(atol, rtol * z.norm()));
  // END
  return z;
}
/* SAM_LISTING_END_2 */
