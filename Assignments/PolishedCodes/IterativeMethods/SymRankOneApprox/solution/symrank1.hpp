#pragma once

#include <Eigen/SVD>
#include <iostream>

/**
 * @brief Computes a symmetric rank-1 matrix closest in Frobenius norm to M.
 *
 * @param M symmetric input matrix
 * @return Eigen::VectorXd z s.t. zz^T is the rank-1 matrix that is closest in
 * Frobenius norm to M
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd symRankOneBestApproxSym(const Eigen::MatrixXd &M) {
  // ensure that M is symmetric
  assert(M.cols() == M.rows() && M == M.transpose() && "M is not symmetric");
  // Vector to be used to return result
  Eigen::VectorXd z = Eigen::VectorXd::Zero(M.cols());
  // TODO: (8-14.a)
  //  solve equ. (8.14.1) for a symmetric matrix with the aid of singular value
  //  decomposition (SVD)
  // START
  // Compute economical SVD of M
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // Extract first column of u and first singular value
  Eigen::VectorXd first_U = svd.matrixU().col(0);
  Eigen::VectorXd sv = svd.singularValues();
  // Scaling
  z = first_U * std::sqrt(sv[0]);
  // END
  return z;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Evaluates the expression (8.14.14).
 *
 * @param v
 * @param b
 * @return Eigen::VectorXd
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd computeKronProdVecMult(const Eigen::VectorXd &v,
                                       const Eigen::VectorXd &b) {
  const unsigned int n = v.size();
  assert(b.size() == n * n && "Size mismatch");
  // Vector to be used to return result
  Eigen::VectorXd res = Eigen::VectorXd::Zero(n);
  // TODO: (8-14.f)
  //  evaluate the expression (8.14.14) efficiently with the aid of Eigen::map
  // START
  // Raw data array for argument vector b
  const double *data = b.data();
  // Twice reinterpretation of data, no copy
  res = Eigen::MatrixXd::Map(data, n, n) * v +
        Eigen::MatrixXd::Map(data, n, n).transpose() * v;
  // END
  return res;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Computes a symmetric rank-1 matrix closest in Frobenius norm to M.
 *
 * @param M matrix M
 * @param rtol relative tolerance for stopping the iteration
 * @param atol absolute tolerance for stopping the iteration
 * @return Eigen::VectorXd
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd symmRankOneApprox(const Eigen::MatrixXd &M,
                                  const double rtol = 1e-6,
                                  const double atol = 1e-8) {
  // ensure M is square
  assert(M.cols() == M.rows() && "Matrix must be square");

  const unsigned int n = M.cols();
  // Initial guess: Column of symmetric part of M with largest norm
  auto Msym = 0.5 * (M + M.transpose());
  Eigen::Index j, dummy;
  Msym.colwise().norm().maxCoeff(&dummy, &j);
  // Set initial guess. This vector is also used to return the result
  Eigen::VectorXd z{Msym.col(j)};
  // TODO: (8-14.g)
  //  solve equ. (8.14.1) using Gauss-Newton iteration with tolerances rtol and
  //  atol
  // START
  // Gauss-Newton correction vector
  Eigen::VectorXd s(n);
  // Direct access to the array storing the entries of M
  const double *dataM = M.data();

  std::cout << " Correction s to current newton step:\n";
  // Main Gauss-Newton loop with correction-based termination
  do {
    // Compute b = vec(M) - kron(z,z)
    const Eigen::MatrixXd Z = z * z.transpose();
    const double *dataz = Z.data();
    const Eigen::VectorXd b = Eigen::MatrixXd::Map(dataM, n * n, 1) -
                              Eigen::MatrixXd::Map(dataz, n * n, 1);
    // Solve normal equations
    const double normalize = 2 * z.squaredNorm();
    Eigen::VectorXd kronMult = computeKronProdVecMult(z, b);
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
