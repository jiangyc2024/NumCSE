#ifndef GRIDFUN_HPP
#define GRIDFUN_HPP

//
// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
// Contributors: tille, jgacon, dcasati
// This file is part of the NumCSE repository.
//

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <vector>

// We use this as type of each index (will be int or similar)
using index_t = Eigen::Index;
// Use this to contain the "size" of our matrix
using shape_t = std::array<index_t, 2>;

/**
 * @brief Matrix filling function
 *
 * @param n number of rows
 * @param m number of columns
 * @return std::function<double(index_t, index_t)> matrix filler function
 */
/* SAM_LISTING_BEGIN_0 */
std::function<double(index_t, index_t)> define_f(unsigned int n,
                                                 unsigned int m) {
  // TODO: (2-16.c) Define a lambda function f that takes as input n,m
  // and outputs a double
  // START
  auto f = [n, m](index_t i, index_t j) -> double {
    if (i > n / 4 && i < n * 3 / 4 && j > m / 4 && j < m * 3 / 4)
      return 1.;
    else
      return 0.;
  };
  return f;
  // END
}
/* SAM_LISTING_END_0 */

/**
 * @brief Evaluate function $f$ at the indices of $\mathbf{X}$.
 * Computes $(\mathbf{X})_{i,j} = f(i,j)$
 * ($+1$ added because C++-indices start at $0$).
 *
 * @param X Matrix $\mathbf{X}$, must have correct size.
 * @param f A function $\mathbf{N}^2 \rightarrow \mathbf{R}$.
 */
/* SAM_LISTING_BEGIN_2 */
void eval(Eigen::MatrixXd& X, std::function<double(index_t, index_t)> f) {
  // TODO: (2-16.c) fill the matrix X with values from f.
  // Hint: do not forget the shift of indices in C++
  // START
  for (index_t i = 0; i < static_cast<index_t>(X.rows()); ++i) {
    for (index_t j = 0; j < static_cast<index_t>(X.cols()); ++j) {
      // Notice that we shift indices due to C++ indexing.
      X(i, j) = f(i + 1, j + 1);
    }
  }
  // END
}
/* SAM_LISTING_END_2 */

/**
 * @brief Given "grid" indices i and j, convert to the vectorized index
 * Return $I$ s.t. $(\mathbf{X})_{i,j} = vec(\mathbf{X})_I$.
 *
 * @param i Coordinate in $i$ direction
 * @param j Coordinate in $j$ direction
 * @param size The size of the "grid" matrix (tuple $(n,m)$)
 * @return index_t $I$ of entry $(i,j)$ on the vectorized matrix
 */
/* SAM_LISTING_BEGIN_5 */
inline index_t to_vector_index(index_t i, index_t j, const shape_t& size) {
  index_t index = 0;
  // TODO: (hint in (2-16.d)) for a given pair of indices (i,j)
  // return the index of vec(X) corresponding to (X)_{i,j}
  // START
  index = size[1] * i + j;
  // END
  return index;
}
/* SAM_LISTING_END_5 */

/**
 * @brief Build sparse matrix constructed from stencil matrix $\mathbf{S}$.
 * The matrix is constructed from triplets, then compressed and returned.
 *
 * @param S The stencil matrix: i.e. an $3 \times 3$ matrix
 * @param size The size of the grid, i.e. the tuple $(n,m)$.
 * @return Eigen::SparseMatrix<double> $\mathbf{A}$ representing linear operator
 * $L$.
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> build_matrix(const Eigen::Matrix3d& S,
                                         const shape_t& size) {
  // Will be used to store triplets
  std::vector<Eigen::Triplet<double>> triplets;
  // The matrix $\mathbf{A}$ will have size $N \times N$.
  // N = n*m
  const index_t N = size[0] * size[1];

  // Now, construct $\mathbf{A}$.
  Eigen::SparseMatrix<double> A(N, N);

  // TODO: (2-16.d) store the matrix in COO format in the vector "triplets"
  // then convert triplets to CRS/CCS format using the functions
  // setFromTriplets() available in Eigen.
  // Hint: use the inline function to_vector_index to find quickly
  // the indices of each triplet.
  // Hint 2: to improve efficiency, you can allocate first a sufficient
  // number of triplets with the function triplets.reserve()

  // START

  // Used to preallocate size of vector "triplets":
  // predicted maximum number of triplets
  const index_t ntriplets =
      (size[0] - 2) * (size[1] - 2) * 3 *
          3                    // elements not at boundary:
                               // those have 3*3 entries per element
      + 2 * 3 * (size[0] - 2)  // top/bottom boundary rows (minus corners)
      + 2 * 3 * (size[1] - 2)  // left/right columns (minus corners)
      + 2 * 2 * 4;             // corners

  // Reserve enough space (more efficient)
  triplets.reserve(ntriplets);

  // Loop over entire grid
  for (index_t i = 0; i < size[0]; ++i) {
    for (index_t j = 0; j < size[1]; ++j) {
      // Loop over the entire stencil matrix $\mathbf{S}$.
      for (index_t k = 0; k < 3; ++k) {
        for (index_t l = 0; l < 3; ++l) {
          if (i + k - 1 >= 0 && i + k - 1 <= size[0] - 1 && j + l - 1 >= 0 &&
              j + l - 1 <= size[1] - 1) {
            // The index of $vec(\mathbf{X})$ w.r.t
            // the indices $(i,j)$ of $\mathbf{X}$
            const index_t I = to_vector_index(i, j, size);
            // Index on $vec(\mathbf{X})$ relative to index
            // $(i+k-2,j+l-2)$
            // NOTE: C++ indexing starts from $0$, so we change the
            // offset.
            const index_t J = to_vector_index(i + k - 1, j + l - 1, size);
            // Add triplet of corresponding stencil matrix entry
            triplets.push_back(Eigen::Triplet<double>(J, I, S(k, l)));
          }
        }
      }
    }
  }

  // Build from triplets and compress matrix (no need to compress through
  // makeCompressed(); already included in setFromTriplets()).
  A.setFromTriplets(triplets.begin(), triplets.end());
  // END

  return A;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute $vec(\mathbf{Y}) := \mathbf{A} * vec(\mathbf{X})$
 * Multiply the sparse matrix by $vec(\mathbf{X})$
 * and store result in $vec(\mathbf{Y})$.
 *
 * @param A A $nm \times nm$ sparse matrix.
 * @param X A $n \times m$ "grid" matrix.
 * @param Y A "grid" matrix s.t.
 * $vec(\mathbf{Y}) := \mathbf{A} * vec(\mathbf{X})$.
 */
/* SAM_LISTING_BEGIN_3 */
void mult(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& X,
          Eigen::MatrixXd& Y) {
  // Size check
  assert(A.rows() == X.rows() * X.cols() && A.cols() == X.rows() * X.cols() &&
         "Inconsistent size of A!");
  // Ensure correct size of $\mathbf{Y}$.
  Y.resizeLike(X);
  // TODO: (2-16.e) compute vec(Y) = A*vec(X)
  // Hint: Use Eigen::Map to reshape X into vec(X).

  // START
  // The second "const" here is needed to "promise" X will not be modified.
  const Eigen::Map<const Eigen::VectorXd> x(X.data(), X.size());
  // Notice how "Map" can also be written onto.
  Eigen::Map<Eigen::VectorXd>(Y.data(), X.size()) = A * x;
  // END
}
/* SAM_LISTING_END_3 */

/**
 * @brief Solve $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$ for
 * $\mathbf{X}$. Solve the linear system given by the sparse matrix $\mathbf{A}$
 * and r.h.s $\mathbf{Y}$.
 *
 * @param A A $nm \times nm$ sparse matrix.
 * @param Y A $n \times m$ "grid" matrix.
 * @param X A "grid" matrix s.t.
 * $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$.
 */
/* SAM_LISTING_BEGIN_4 */
void solve(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& Y,
           Eigen::MatrixXd& X) {
  // Size check
  assert(A.rows() == Y.rows() * Y.cols() && A.cols() == Y.rows() * Y.cols() &&
         "Inconsistent size of A!");
  // Ensure correct size of $\mathbf{X}$.
  X.resizeLike(Y);
  // TODO: (2-16.f) compute X such that vec(Y) = A * vec(X) using a
  // sparse solver in Eigen.
  // START
  // Sparse LU decomposition object
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  // Prepare LU factorization
  solver.analyzePattern(A);
  solver.factorize(A);

  // Check if factorization failed.
  if (solver.info() != Eigen::Success) {
    std::cerr << "Fatal: LU decomposition of A failed!!!" << std::endl
              << "Is the matrix invertible?" << std::endl;
    return;
  }

  // Solve system after factorization.
  // Notice how you can assign to a Map object.
  Eigen::Map<Eigen::VectorXd>(X.data(), Y.size()) =
      solver.solve(Eigen::Map<const Eigen::VectorXd>(Y.data(), Y.size()));
  // END
}
/* SAM_LISTING_END_4 */

#endif