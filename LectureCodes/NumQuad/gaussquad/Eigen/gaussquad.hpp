#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/* Golub-Welsch algorithm for computing Gauss points */

/* SAM_LISTING_BEGIN_0 */
struct QuadRule {
  Eigen::VectorXd nodes_, weights_;
} __attribute__((aligned(32)));

inline QuadRule gaussquad(const unsigned int n) {
  QuadRule qr;
  // Symmetric matrix whose eigenvalues provide Gauss points
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  for (unsigned int i = 1; i < n; ++i) {  // \Label[line]{gw:3}
    const double b = i / std::sqrt(4. * i * i - 1.);
    M(i, i - 1) = M(i - 1, i) = b;
  }  // \Label[line]{gw:3x}
  // using Eigen's built-in solver for symmetric eigenvalue problems
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(M);

  qr.nodes_ = eig.eigenvalues();  //
  qr.weights_ =
      2 * eig.eigenvectors().topRows<1>().array().pow(2);  // \Label[line]{gw:4}
  return qr;
}
/* SAM_LISTING_END_0 */
