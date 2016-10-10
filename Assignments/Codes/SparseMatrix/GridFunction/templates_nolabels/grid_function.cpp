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


#include <mgl2/mgl.h>

using namespace Eigen;

// We use this as type of each index (will be int or similar)
using index_t = std::size_t;
// Use this to contain the "shape" of our matrix
using shape_t =  std::array<index_t, 2>;

/* \brief Given "grid" indices i and j, convert to the vectorie index
 *
 * \param[in] i Coordinate in i direction
 * \param[in] j Coordinate in j direction
 * \param[in] size The shape of the "grid" matrix
 * \return Index I of point (i,j) on the vectorized matrix
 */
inline index_t to_vector_index(index_t i,
                               index_t j,
                               const shape_t & size) {
    // TODO: convert grid index (i,j) to index I and return
    return 0;
}

/* \brief BUild a sparse matrix constructed from the stencil S
 * The matrix is constructed from triplets, compressed and returned.
 * \tparam n Size of stencil (odd number)
 * \param[in] S The stencil: a $n \times n$ matrix
 * \param[in] size The size if the grid
 * \return Sparse matrix representing linear operator $L$.
 */
template <int s>
SparseMatrix<double> build_matrix(const Matrix<double, s, s> & S,
                                 const shape_t & size) {
    static_assert( s % 2 == 1, "Fatal: n must be odd!");

    // Will be used to store triplets
    std::vector<Triplet<double>> triplets;

    // The matrix A wil have size $N \times N$.
    const index_t N = size[0]*size[1];
    // TODO: build sparse matrix A
    return SparseMatrix<double>(N,N);
}

/* \brief Evaluate function f at the coordinates of X and store values in X
 * Computes $(X)_{i,j} = f(i-1,j-1)$ (-1 subracted because indices start at 0).
 * \param[in,out] X Matrix where values of X will be stored, must have correct size
 * \param[in] f A function taking two indices (i,j) and returning a value.
 */
void eval(MatrixXd & X, std::function<double(index_t,index_t)> f) {
    // TODO:implement an evaluation function
}

/* \brief Computes vec(Y) = A * vec(X)
 * Multiply the sparse matrix times vec(X) and store result in vec(Y)
 * \param[in] A a n^2*m^2 sparse marix
 * \param[in] X a n*m "grid" matrix
 * \param[out] Y a n*m "grid" matrix s.t. vec(Y) = A*vec(X)
 */
void mult(const SparseMatrix<double> & A, const MatrixXd & X, MatrixXd & Y) {
    assert( A.rows() == X.rows()*X.cols() &&
            A.cols() == X.rows()*X.cols() &&
            "Inconsistent size of A!");

    Y.resizeLike(X);

    // TODO: compute vec(Y) = A*vec(X)
}

/* \brief Solves vec(Y) = A * vec(X) for X
 * Solve the linear system given by the sparse matrix A and r.h.s Y
 * \param[in] A a n^2*m^2 sparse marix
 * \param[in] Y a n*m "grid" matrix
 * \param[out] X a n*m "grid" matrix s.t. vec(Y) = A*vec(X)
 */
void solve(const SparseMatrix<double> & A, const MatrixXd & Y, MatrixXd & X) {
    assert( A.rows() == Y.rows()*Y.cols() &&
            A.cols() == Y.rows()*Y.cols() &&
            "Inconsistent size of A!");

    X.resizeLike(Y);

    // Sparse LU decomposition object
    SparseLU<SparseMatrix<double>> solver;

    // TODO: solve system vec(Y) = A*vec(X)
}

int main(int argc, char** argv) {
    // Size of the *grid" matrices X and Y
    index_t n = 100, m = 100;
    // If provided, will be read from command line
    if(argc > 1) {
        n = std::stoi(argv[1]);
    }
    if(argc > 2) {
        m = std::stoi(argv[2]);
    }

    // The stencil
    Matrix3d S;
    S << 0, 1, 0,
         1, -4, 1,
         0, 1, 0;

    // Construct S
    std::cout << "--> Building matrix A..." << std::endl;
    SparseMatrix<double> A = build_matrix(S, shape_t{n, m});

    std::cout << A;

    // The grid matrices X and Y
    MatrixXd X(n,m), Y;

    // A function of indices (i,j) of the matrix X
    auto f = [n,m] (index_t i , index_t j) {
        if(i > n/4 && i < n*3/4 && j > m/4 && j < m*3/4 ) return 1.;
        else return 0.;
    };

    // Fill entries of X given f
    std::cout << "--> Evaluating f at the indices of X..." << std::endl;
    eval(X, f);

    // Compute A*vec(X)
    std::cout << "--> Computing A*vec(X)..." << std::endl;
    mult(A, X, Y);

    // Compute A^{-1}*vec(Y)
    std::cout << "--> Solving A*vec(X) = vec(Y)..." << std::endl;
//    solve(A, X, Y);
//    Y = X;

    std::cout << Y;

    // Plot values of X
    mglData Xd(Y.rows(), Y.cols(), Y.data());

    mglGraph gr;
    gr.SubPlot(1,1,0,"<_");
    gr.SetRanges(0,n,0,m);
    gr.SetRange('c', -1, 0);
    gr.Colorbar("bcwyr");
    gr.Title("Visualization of $X$");
    gr.Axis();
    gr.Tile(Xd, "bcwyr");
    gr.WriteEPS("grid.eps");
}
