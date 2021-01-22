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
  // Compute economical SVD of M
  JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
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
  const unsigned int n = v.size();
  assert(b.size() == n * n && "Size mismatch");
  // Vector to be used to return result
  VectorXd res;
  // TODO: (9-14.f)
  //  evaluate the expression (9.14.14) efficiently with the aid of Eigen::map
  // START
  // Raw data array for argument vector b
  const double *data = b.data();
  // Twice reinterpretation of data, no copy
  res =
      MatrixXd::Map(data, n, n) * v + MatrixXd::Map(data, n, n).transpose() * v;
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
  // Gauss-Newton correction vector
  VectorXd s(n);
  // Direct access to the array storing the entries of M
  const double *dataM = M.data();

  std::cout << " Correction s to current newton step:\n";
  // Main Gauss-Newton loop with correction-based termination
  do {
    // Compute b = vec(M) - kron(z,z)
    const MatrixXd Z = z * z.transpose();
    const double *dataz = Z.data();
    const VectorXd b =
        MatrixXd::Map(dataM, n * n, 1) - MatrixXd::Map(dataz, n * n, 1);
    // Solve normal equations
    const double normalize = 2 * z.squaredNorm();
    VectorXd kronMult = computeKronProdVecMult(z, b);
    s = (kronMult - z.dot(kronMult) / normalize * z) / normalize;
    // Apply Gauss-Newton update
    z = z + s;
    std::cout << "|s| = " << s.norm() << std::endl;
  } while (s.norm() > std::max(atol, rtol * z.norm()));
  std::cout << "\n";
  // END
  return z;
}
/* SAM_LISTING_END_2 */
