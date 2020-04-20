#ifndef GRIDFUN_HPP
#define GRIDFUN_HPP

//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <vector>
#include <array>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

using namespace Eigen;



// We use this as type of each index (will be int or similar)
using index_t = std::size_t;
// Use this to contain the "size" of our matrix
using shape_t =  std::array<index_t, 2>;




/* SAM_LISTING_BEGIN_0 */
std::function<double (index_t,index_t)> define_f(int n, int m) {
    // TO DO (3-16.c): Define a lambda function f that takes as input n,m
    // and outputs a double 
    //START

    return [](index_t x,index_t y){ return NAN; };
    //END
}
/* SAM_LISTING_END_0 */

/* \brief Evaluate function $f$ at the indices of $\mathbf{X}$.
 * Computes $(\mathbf{X})_{i,j} = f(i,j)$
 * ($+1$ added because C++-indices start at $0$).
 * \param[in,out] X Matrix $\mathbf{X}$, must have correct size.
 * \param[in] f A function $\mathbf{N}^2 \rightarrow \mathbf{R}$.
 */
/* SAM_LISTING_BEGIN_2 */
void eval(MatrixXd & X, std::function<double(index_t,index_t)> f) {
  //TO DO (3-16.c): fill the matrix X with values from f. 
  // Hint: do not forget the shift of indices in C++ 
  //START
  
  //END
}
/* SAM_LISTING_END_2 */

/* \brief Given "grid" indices i and j, convert to the vectoried index
 * Return $I$ s.t. $(\mathbf{X})_{i,j} = vec(\mathbf{X})_I$.
 * \param[in] i Coordinate in $i$ direction
 * \param[in] j Coordinate in $j$ direction
 * \param[in] size The size of the "grid" matrix (tuple $(n,m)$)
 * \return Index $I$ of entry $(i,j)$ on the vectorized matrix
 */
inline index_t to_vector_index(index_t i,
                               index_t j,
                               const shape_t & size) {
    //TO DO, hint in (3-16.d) : for a given pair of indices (i,j) 
    // return the index of vec(X) corresponding to (X)_{i,j}                         
    //START                             
    return  0;
    //END
}

/* \brief Build sparse matrix constructed from stencil matrix $\mathbf{S}$.
 * The matrix is constructed from triplets, then compressed and returned.
 * \param[in] S The stencil matrix: i.e. an $3 \times 3$ matrix
 * \param[in] size The size of the grid, i.e. the tuple $(n,m)$.
 * \return Sparse matrix $\mathbf{A}$ representing linear operator $L$.
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> build_matrix(const Matrix3d & S,
                                 const shape_t & size) {
    // Will be used to store triplets
    std::vector<Triplet<double> > triplets;
    // The matrix $\mathbf{A}$ will have size $N \times N$.
    // N = n*m
    const index_t N = size[0]*size[1];

    // Now, construct $\mathbf{A}$.
    SparseMatrix<double> A(N,N);
    
    // TO DO: (3-16.d) store the matrix in COO format in the vector "triplets"
    // then convert triplets to CRS/CCS format using the functions 
    // setFromTriplets() and makeCompressed() available in Eigen.
    // Hint: use the inline function to_vector_index to find quickly 
    // the indices of each triplet.
    // Hint 2: to improve efficiency, you can allocate first a sufficient 
    // number of triplets with the function triplets.reserve()
    
    //START
    
    //END
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
void mult(const SparseMatrix<double> & A,
          const MatrixXd & X, MatrixXd & Y) {
    // Size check
    assert( A.rows() == X.rows()*X.cols() &&
            A.cols() == X.rows()*X.cols() &&
            "Inconsistent size of A!");
    // Ensure correct size of $\mathbf{Y}$.
    Y.resizeLike(X);
    //TO DO (3-16.e): compute vec(Y) = A*vec(X)
    // Hint: Use Eigen::Map to reshape X into vec(X).
    
    //START
    
    //END
}
/* SAM_LISTING_END_3 */

/* \brief Solve $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$ for $\mathbf{X}$.
 * Solve the linear system given by the sparse matrix $\mathbf{A}$
 * and r.h.s $\mathbf{Y}$.
 * \param[in] A A $nm \times nm$ sparse marix.
 * \param[in] Y A $n \times m$ "grid" matrix.
 * \param[out] X A "grid" matrix s.t.
 * $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$.
 */
/* SAM_LISTING_BEGIN_4 */
void solve(const SparseMatrix<double> & A,
           const MatrixXd & Y, MatrixXd & X) {
    // Size check
    assert( A.rows() == Y.rows()*Y.cols() &&
            A.cols() == Y.rows()*Y.cols() &&
            "Inconsistent size of A!");
    // Ensure correct size of $\mathbf{X}$.
    X.resizeLike(Y);
    //TO DO (3-16.f): compute X such that vec(Y) = A * vec(X) using a
    // sparse solver in Eigen.
    //START
    
    //END
}
/* SAM_LISTING_END_4 */

#endif