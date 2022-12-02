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
Eigen::VectorXd multAx(Eigen::VectorXd &a, const Eigen::VectorXd &b,
                       const Eigen::VectorXd &x) {
  const unsigned int n = x.size();
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return y;
  }

  // TODO: (2-6.a) Multiply A*x and exploit the diagonal structure of the
  // matrix.
  // START
  // Handle first two rows
  if (n > 1) {
    y(0) = 2 * x(0) + a(0) * x(1);
  } else {
    y(0) = 2 * x(0);
    return y;
  }
  if (n > 2) {
    y(1) = 2 * x(1) + a(1) * x(2);
  } else {
    y(1) = 2 * x(1);
    return y;
  }

  // Row if $n > 3$ and without last row
  for (unsigned int i = 2; i < n - 1; ++i) {
    y(i) = 2 * x(i) + b(i - 2) * x(i - 2) + a(i) * x(i + 1);
  }

  // Last row special case
  if (n > 2)
    y(n - 1) = 2 * x(n - 1) + b(n - 3) * x(n - 3);
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
Eigen::VectorXd solvelseAupper(const Eigen::VectorXd &a,
                               const Eigen::VectorXd &r) {
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
  // Backward substitution
  x(n - 1) = 0.5 * r(n - 1);
  for (int j = n - 2; j >= 0; --j) {
    x(j) = 0.5 * (r(j) - a(j) * x(j + 1));
  }
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
Eigen::VectorXd solvelseA(const Eigen::VectorXd &a, const Eigen::VectorXd &b,
                          const Eigen::VectorXd &r) {
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
  // Initializing vectors for the diagonal and lower off-diagonal
  Eigen::VectorXd d = Eigen::VectorXd::Constant(n, 2);
  Eigen::VectorXd c = Eigen::VectorXd::Constant(n - 1, 0);

  // A modifiable RHS vector
  Eigen::VectorXd R(r);

  // Performing Gaussian elimination
  // We only need to update the vectors c and d!
  for (unsigned i = 0; i < n - 2; ++i) {
    /*

      ... d(i)  a(i)     0       0      ...
      ... c(i)  d(i+1)   a(i+1)  0      ...
      ... b(i)  c(i+1)   d(i+2)  a(i+2) 0

                    ===>

      ... d(i)  a(i)                    0       0      ...
      ... 0     d(i+1)-c(i)/d(i)*a(i)   a(i+1)  0      ...
      ... b(i)  c(i+1)                  d(i+2)  a(i+2) 0
    */
    d(i + 1) -= c(i) / d(i) * a(i);
    // Modifying the RHS accordingly
    R(i + 1) -= c(i) / d(i) * R(i);

    /*
      ... d(i)  a(i)                    0       0      ...
      ... 0     d(i+1)-c(i)/d(i)*a(i)   a(i+1)  0      ...
      ... b(i)  c(i+1)                  d(i+2)  a(i+2) 0

                    ===>

      ... d(i)  a(i)                    0       0      ...
      ... 0     d(i+1)-c(i)/d(i)*a(i)   a(i+1)  0      ...
      ... 0     c(i+1)-b(i)/d(i)*a(i)   d(i+2)  a(i+2) 0
    */
    c(i + 1) -= b(i) / d(i) * a(i);
    // Modifying the RHS accordingly
    R(i + 2) -= b(i) / d(i) * R(i);
  }
  // Last row
  d(n - 1) -= c(n - 2) / d(n - 2) * a(n - 2);
  R(n - 1) -= c(n - 2) / d(n - 2) * R(n - 2);

  // Backward substitution
  x(n - 1) = R(n - 1) / d(n - 1);
  for (int i = n - 2; i >= 0; --i) {
    x(i) = (R(i) - a(i) * x(i + 1)) / d(i);
  }
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
Eigen::VectorXd solvelseASparse(const Eigen::VectorXd &a,
                                const Eigen::VectorXd &b,
                                const Eigen::VectorXd &r) {
  // Set up dimensions
  const unsigned int n = r.size();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return x;
  }

  // TODO: (2-6.f) Solve A*x=r using Eigen's sparse solver.
  // START
  // Fill in matrix:
  // We reserve three nonzero entries per row for Gaussian fill-in
  Eigen::SparseMatrix<double> A(n, n);
  A.reserve(3*n);
  for (unsigned int i = 0; i < n; ++i) {
    A.insert(i, i) = 2;
    if (i < n - 1)
      A.insert(i, i + 1) = a(i);
    if (i >= 2)
      A.insert(i, i - 2) = b(i - 2);
  }
  A.makeCompressed();

  // Call SparseLU
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(A);
  solver.compute(A);
  x = solver.solve(r);
  // END
  return x;
}
/* SAM_LISTING_END_3 */

#endif
