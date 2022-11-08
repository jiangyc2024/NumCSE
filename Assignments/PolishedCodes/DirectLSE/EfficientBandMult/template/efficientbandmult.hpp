#ifndef EFFICIENTBANDMULT_HPP
#define EFFICIENTBANDMULT_HPP

#include <Eigen/Sparse>
#include <iostream>

/**
 * @brief Compute $y = A*x$ with A banded matrix with diagonal structure
 *
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param b An $(n-2)$-dimensional vector for the second lower diagonal
 * @param x An $n$-dimensional vector for $Ax = y$
 * @return The $n$-dimensional vector $y = Ax$
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd multAx(Eigen::VectorXd& a, const Eigen::VectorXd& b,
                       const Eigen::VectorXd& x) {
  const unsigned int n = x.size();
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return y;
  }

  // TODO: (2-6.a) Multiply A*x and exploit the diagonal structure of the
  // matrix.
  // START

  // END
  return y;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Solve $r = A*x$ with $A$ banded matrix with upper triangular sparse
 * structure
 *
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param r An $n$-dimensional vector for $Ax = r$
 * @return The $n$-dimensional vector from $Ax = r$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solvelseAupper(const Eigen::VectorXd& a,
                               const Eigen::VectorXd& r) {
  // Set up dimensions
  const unsigned int n = r.size();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  if (a.size() < n - 1) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return x;
  }

  // TODO: (2-6.c) Solve the LSE A*x = r for x by exploiting the diagonal
  // structure of A and the fact that b is zero.
  // START

  // END
  return x;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solve $r = A*x$ with $A$ banded matrix using Gaussian elimination (no
 * pivot)
 *
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param b An $(n-2)$-dimensional vector for the second lower diagonal
 * @param r An $n$-dimensional vector for $Ax = r$
 * @return The $n$-dimensional vector from $Ax = r$
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solvelseA(const Eigen::VectorXd& a, const Eigen::VectorXd& b,
                          const Eigen::VectorXd& r) {
  // Set up dimensions
  const unsigned int n = r.size();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return x;
  }

  // TODO: (2-6.d) Compute the solution to A*x=r by Gaussian elimination and
  // without the help of Eigen's solvers.
  // START

  // END
  return x;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Solve $r = A*x$ with $A$ banded matrix using Eigen::SparseLU
 *
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param b An $(n-2)$-dimensional vector for the second lower diagonal
 * @param r An $n$-dimensional vector for $Ax = r$
 * @return The $n$-dimensional vector from $Ax = r$
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd solvelseASparse(const Eigen::VectorXd& a,
                                const Eigen::VectorXd& b,
                                const Eigen::VectorXd& r) {
  // Set up dimensions
  const unsigned int n = r.size();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return x;
  }

  // TODO: (2-6.f) Solve A*x=r using Eigen's sparse solver.
  // START

  // END
  return x;
}
/* SAM_LISTING_END_3 */

#endif
