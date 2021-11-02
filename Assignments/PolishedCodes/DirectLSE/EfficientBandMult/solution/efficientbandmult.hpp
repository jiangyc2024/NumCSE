#ifndef EFFICIENTBANDMULT_HPP
#define EFFICIENTBANDMULT_HPP

#include <Eigen/Sparse>
#include <iostream>

/**
 * @brief Compute $y = A*x$ with A banded matrix with diagonal structure
 *
 * @tparam Vector random-access container, must provide size(), operator(),
 * resize()
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param b An $(n-2)$-dimensional vector for the second lower diagonal
 * @param x An $n$-dimensional vector for $Ax = y$
 * @param y The $n$-dimensional vector $y = Ax$
 */
/* SAM_LISTING_BEGIN_0 */
template <class Vector>
void multAx(const Vector& a, const Vector& b, const Vector& x, Vector& y) {
  const unsigned int n = x.size();
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return;
  }
  y.resize(n);

  // TODO: (2-6.a) Multiply A*x and exploit the diagonal structure of the
  // matrix. START Handle first two rows
  if (n > 1) {
    y(0) = 2 * x(0) + a(0) * x(1);
  } else {
    y(0) = 2 * x(0);
    return;
  }
  if (n > 2) {
    y(1) = 2 * x(1) + a(1) * x(2);
  } else {
    y(1) = 2 * x(1);
    return;
  }

  // Row if $n > 3$ and without last row
  for (unsigned int i = 2; i < n - 1; ++i) {
    y(i) = 2 * x(i) + b(i - 2) * x(i - 2) + a(i) * x(i + 1);
  }

  // Last row special case
  if (n > 2) y(n - 1) = 2 * x(n - 1) + b(n - 3) * x(n - 3);
  // END
}
/* SAM_LISTING_END_0 */

/**
 * @brief Solve $r = A*x$ with $A$ banded matrix with upper triangular sparse
 * structure
 *
 * @tparam Vector random-access container, must provide size(), operator(),
 * resize()
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param r An $n$-dimensional vector for $Ax = r$
 * @param x The $n$-dimensional vector from $Ax = r$
 */
/* SAM_LISTING_BEGIN_1 */
template <class Vector>
void solvelseAupper(const Vector& a, const Vector& r, Vector& x) {
  // Set up dimensions
  const unsigned int n = r.size();
  if (a.size() < n - 1) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return;
  }
  x.resize(n);

  // TODO: (2-6.c) Solve the LSE A*x = r for x by exploiting the diagonal
  // structure of A and the fact that b is zero.
  // START
  // Backward substitution
  x(n - 1) = 0.5 * r(n - 1);
  for (int j = n - 2; j >= 0; --j) {
    x(j) = 0.5 * (r(j) - a(j) * x(j + 1));
  }
  // END
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solve $r = A*x$ with $A$ banded matrix using Gaussian elimination (no
 * pivot)
 *
 * @tparam Vector random-access container, must provide size(), operator(),
 * resize()
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param b An $(n-2)$-dimensional vector for the second lower diagonal
 * @param r An $n$-dimensional vector for $Ax = r$
 * @param x The $n$-dimensional vector from $Ax = r$
 */
/* SAM_LISTING_BEGIN_2 */
template <class Vector>
void solvelseA(const Vector& a, const Vector& b, const Vector& r, Vector& x) {
  // Set up dimensions
  const unsigned int n = r.size();
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return;
  }

  // TODO: (2-6.d) Compute the solution to A*x=r by Gaussian elimination and
  // without the help of Eigen's solvers.
  // START
  Vector c, d, y;
  c.resize(n - 1);
  d.resize(n);
  for (unsigned int i = 0; i < n - 1; ++i) {
    c(i) = 0;
    d(i) = 2;
  }
  d(n - 1) = 2;
  x = r;
  y = r;

  // Plain vectors are enough
  // Fill-in is confined to a single lower off-diagonal, which can be held in
  // vector $c$
  for (unsigned int i = 0; i < n - 2; ++i) {
    c(i + 1) = -b(i) / d(i) * a(i);
    d(i + 1) -= c(i) / d(i) * a(i);
    y(i + 1) -= c(i) / d(i) * y(i);
    y(i + 2) -= b(i) / d(i) * y(i);
  }

  x(n - 1) = y(n - 1) / d(n - 1);
  for (int i = n - 2; i >= 0; --i) {
    x(i) = (y(i) - a(i) * x(i + 1)) / d(i);
  }
  // END
}
/* SAM_LISTING_END_2 */

/**
 * @brief Solve $r = A*x$ with $A$ banded matrix using Eigen::SparseLU
 *
 * @tparam Vector random-access container, must provide size(), operator(),
 * resize()
 * @param a An $(n-1)$-dimensional vector for the upper diagonal
 * @param b An $(n-2)$-dimensional vector for the second lower diagonal
 * @param r An $n$-dimensional vector for $Ax = r$
 * @param x The $n$-dimensional vector from $Ax = r$
 */
/* SAM_LISTING_BEGIN_3 */
template <class Vector>
void solvelseAEigen(const Vector& a, const Vector& b, const Vector& r,
                    Vector& x) {
  // Set up dimensions
  typedef typename Vector::Scalar Scalar;
  const unsigned int n = r.size();
  if (a.size() < n - 1 || b.size() < n - 2) {
    std::cerr << "Error: size mismatch!" << std::endl;
    return;
  }

  // TODO: (2-6.f) Solve A*x=r using Eigen's sparse solver.
  // START
  // Fill in matrix:
  // We reserve three nonzero entries per row for Gaussian fill-in
  Eigen::SparseMatrix<Scalar> A(n, n);
  A.reserve(3);
  for (unsigned int i = 0; i < n; ++i) {
    A.insert(i, i) = 2;
    if (i < n - 1) A.insert(i, i + 1) = a(i);
    if (i >= 2) A.insert(i, i - 2) = b(i - 2);
  }
  A.makeCompressed();

  // Call SparseLU
  Eigen::SparseLU<Eigen::SparseMatrix<Scalar>> solver;
  solver.analyzePattern(A);
  solver.compute(A);
  x = solver.solve(r);
  // END
}
/* SAM_LISTING_END_3 */

#endif
