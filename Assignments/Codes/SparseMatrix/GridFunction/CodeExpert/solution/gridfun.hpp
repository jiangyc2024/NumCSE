#ifndef GRIDFUN_HPP
#define GRIDFUN_HPP

////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <array>
#include <functional>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

using namespace Eigen;

// We use this as type of each index (will be int or similar)
using index_t = std::size_t;
// Use this to contain the "size" of our matrix
using shape_t = std::array<index_t, 2>;

/* SAM_LISTING_BEGIN_0 */
std::function<double(index_t, index_t)> define_f(int n, int m) {
  // TO DO (3-16.c): Define a lambda function f that takes as input n,m
  // and outputs a double
  // START
  auto f = [n, m](index_t i, index_t j) {
    if (i > n / 4 && i < n * 3 / 4 && j > m / 4 && j < m * 3 / 4)
      return 1.;
    else
      return 0.;
  };
  // END
  return f;
}
/* SAM_LISTING_END_0 */

/* \brief Evaluate function $f$ at the indices of $\mathbf{X}$.
 * Computes $(\mathbf{X})_{i,j} = f(i,j)$
 * ($+1$ added because C++-indices start at $0$).
 * \param[in,out] X Matrix $\mathbf{X}$, must have correct size.
 * \param[in] f A function $\mathbf{N}^2 \rightarrow \mathbf{R}$.
 */
/* SAM_LISTING_BEGIN_2 */
void eval(MatrixXd &X, std::function<double(index_t, index_t)> f) {
  // TO DO (3-16.c): fill the matrix X with values from f.
  // Hint: do not forget the shift of indices in C++
  // START
  for (index_t i = 0; i < (index_t)X.rows(); ++i) {
    for (index_t j = 0; j < (index_t)X.cols(); ++j) {
      // Notice that we shift indices due to C++ indexing.
      X(i, j) = f(i + 1, j + 1);
    }
  }
  // END
}
/* SAM_LISTING_END_2 */

/* \brief Given "grid" indices i and j, convert to the vectoried index
 * Return $I$ s.t. $(\mathbf{X})_{i,j} = vec(\mathbf{X})_I$.
 * \param[in] i Coordinate in $i$ direction
 * \param[in] j Coordinate in $j$ direction
 * \param[in] size The size of the "grid" matrix (tuple $(n,m)$)
 * \return Index $I$ of entry $(i,j)$ on the vectorized matrix
 */
inline index_t to_vector_index(index_t i, index_t j, const shape_t &size) {
  // TO DO, hint in (3-16.d) : for a given pair of indices (i,j)
  // return the index of vec(X) corresponding to (X)_{i,j}
  // START
  return size[1] * i + j;
  // END
}

/* \brief Build sparse matrix constructed from stencil matrix $\mathbf{S}$.
 * The matrix is constructed from triplets, then compressed and returned.
 * \param[in] S The stencil matrix: i.e. an $3 \times 3$ matrix
 * \param[in] size The size of the grid, i.e. the tuple $(n,m)$.
 * \return Sparse matrix $\mathbf{A}$ representing linear operator $L$.
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> build_matrix(const Matrix3d &S, const shape_t &size) {
  // Will be used to store triplets
  std::vector<Triplet<double>> triplets;
  // The matrix $\mathbf{A}$ will have size $N \times N$.
  // N = n*m
  const index_t N = size[0] * size[1];

  // Now, construct $\mathbf{A}$.
  SparseMatrix<double> A(N, N);

  // TO DO: (3-16.d) store the matrix in COO format in the vector "triplets"
  // then convert triplets to CRS/CCS format using the functions
  // setFromTriplets() and makeCompressed() available in Eigen.
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
      + 2 * 3 * (size[0] - 2); // top/bottom boundary rows (minus corners)
  +2 * 3 * (size[1] - 2);      // left/right columns (minus corners)
  +2 * 2 * 4;                  // corners

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
            triplets.push_back(Triplet<double>(J, I, S(k, l)));
          }
        }
      }
    }
  }

  // Build from triplets and compress matrix.
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();

  // END
  return A;
}
/* SAM_LISTING_END_1 */

/* \brief Compute $vec(\mathbf{Y}) := \mathbf{A} * vec(\mathbf{X})$
 * Multiply the sparse matrix by $vec(\mathbf{X})$
 * and store result in $vec(\mathbf{Y})$.
 * \param[in] A A $nm \times nm$ sparse marix.
 * \param[in] X A $n \times m$ "grid" matrix.
 * \param[out] Y A "grid" matrix s.t.
 * $vec(\mathbf{Y}) := \mathbf{A} * vec(\mathbf{X})$.
 */
/* SAM_LISTING_BEGIN_3 */
void mult(const SparseMatrix<double> &A, const MatrixXd &X, MatrixXd &Y) {
  // Size check
  assert(A.rows() == X.rows() * X.cols() && A.cols() == X.rows() * X.cols() &&
         "Inconsistent size of A!");
  // Ensure correct size of $\mathbf{Y}$.
  Y.resizeLike(X);
  // TO DO (3-16.e): compute vec(Y) = A*vec(X)
  // Hint: Use Eigen::Map to reshape X into vec(X).

  // START
  // The second "const" here is needed to "promise" X will not be modified.
  const Map<const VectorXd> x(X.data(), X.size());
  // Notice how "Map" can also be written onto.
  Map<VectorXd>(Y.data(), X.size()) = A * x;

  // END
}
/* SAM_LISTING_END_3 */

/* \brief Solve $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$ for
 * $\mathbf{X}$. Solve the linear system given by the sparse matrix $\mathbf{A}$
 * and r.h.s $\mathbf{Y}$.
 * \param[in] A A $nm \times nm$ sparse marix.
 * \param[in] Y A $n \times m$ "grid" matrix.
 * \param[out] X A "grid" matrix s.t.
 * $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$.
 */
/* SAM_LISTING_BEGIN_4 */
void solve(const SparseMatrix<double> &A, const MatrixXd &Y, MatrixXd &X) {
  // Size check
  assert(A.rows() == Y.rows() * Y.cols() && A.cols() == Y.rows() * Y.cols() &&
         "Inconsistent size of A!");
  // Ensure correct size of $\mathbf{X}$.
  X.resizeLike(Y);
  // TO DO (3-16.f): compute X such that vec(Y) = A * vec(X) using a
  // sparse solver in Eigen.
  // START
  // Sparse LU decomposition object
  SparseLU<SparseMatrix<double>> solver;

  // Prepare LU factorization
  solver.analyzePattern(A);
  solver.factorize(A);

  // Check if factorization failed.
  if (solver.info() != Success) {
    std::cerr << "Fatal: LU decomposition of A failed!!!" << std::endl
              << "Is the matrix invertible?" << std::endl;
    return;
  }

  // Solve system after factorization.
  // Notice how you can assign to a Map object.
  Map<VectorXd>(X.data(), Y.size()) =
      solver.solve(Map<const VectorXd>(Y.data(), Y.size()));
  // END
}
/* SAM_LISTING_END_4 */

#endif
