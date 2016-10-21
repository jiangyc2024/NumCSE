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
// Use this to contain the "size" of our matrix
using shape_t =  std::array<index_t, 2>;

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
    return size[1] * i + j;
}

/* \brief Build sparse matrix constructed from stencil matrix $\mathbf{S}$.
 * The matrix is constructed from triplets, then compressed and returned.
 * \param[in] S The stencil matrix: i.e. an $3 \times 3$ matrix
 * \param[in] size The size of the grid, i.e. the tuple $(n,m)$.
 * \return Sparse matrix $\mathbf{A}$ representing linear operator $L$.
 */
SparseMatrix<double> build_matrix(const Matrix3d & S,
                                 const shape_t & size) {
    // Will be used to store triplets
    std::vector<Triplet<double>> triplets;

    // The matrix $\mathbf{A}$ will have size $N \times N$.
    // N = n*m
    const index_t N = size[0]*size[1];

    // Used to preallocate size of vector "triplets":
    // predicted number of triplets
    const index_t ntriplets =
            (size[0]-2)*(size[1]-2)*3*3 // rows not at boundary:
                // those have 3*3 entries per row
            + 2*size[0] // top/bottom boundary rows
                // 1 entry per row (identity)
            + 2*(size[1] - 2); // left/right rows (minus corners)
                // 1 entry per row (identity)

    // Reserve enough space (more efficient)
    triplets.reserve(ntriplets);

    // Loop over entire  grid
    for(index_t i = 0; i < size[0]; ++i) {
        for(index_t j = 0; j < size[1]; ++j) {
            // The index of $vec(\mathbf{X})$ w.r.t
            // the indices $(i,j)$ of $\mathbf{X}$
            const index_t I = to_vector_index(i, j, size);

            // If stencil falls outside of grid, just use identity
            // (means we are at "boundaries" of $\mathbf{X}$)
            if( i == 0 || i == size[0] - 1 ||
                    j == 0 || j == size[1] - 1 ) {
                // In this case, only entry is diagonal, and has value $1$.
                triplets.push_back(Triplet<double>(I, I, 1.));
                // Do not add furhter entries...
                continue;
            }

            // Loop over the entire stencil matrix $\mathbf{S}$.
            for(index_t k = 0; k < 3; ++k) {
                for(index_t l = 0; l < 3; ++l) {
                    // Index on $vec(\mathbf{X})$ relative to index
                    // $(i+k-2,j+l-2)$
                    // NOTE: C++ indexing starts from $0$, so we change the
                    // offset.
                    const index_t J = to_vector_index(
                                i + k - 1,
                                j + l - 1,
                                size
                                );

                    // Add triplet of corresponging stencil matrix entry
                    triplets.push_back(
                                Triplet<double>(I, J, S(k, l))
                                );
                }
            }
        }
    }

    // Now, construct $\mathbf{A}$.
    SparseMatrix<double> A(N,N);

    // Build from triplets and compress matrix.
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    return A;
}

/* \brief Evaluate function $f$ at the indices of $\mathbf{X}$.
 * Computes $(\mathbf{X})_{i,j} = f(i,j)$
 * ($+1$ added because C++-indices start at $0$).
 * \param[in,out] X Matrix $\mathbf{X}$, must have correct size.
 * \param[in] f A function $\mathbf{N}^2 \rightarrow \mathbf{R}$.
 */
void eval(MatrixXd & X, std::function<double(index_t,index_t)> f) {
    for(index_t i = 0; i < (index_t) X.rows(); ++i) {
        for(index_t j = 0; j < (index_t) X.cols(); ++j) {
            // Notice that we shift indices due to C++ indexing.
            X(i,j) = f(i+1,j+1);
        }
    }
}

/* \brief Compute $vec(\mathbf{Y}) := \mathbf{A} * vec(\mathbf{X})$
 * Multiply the sparse matrix by $vec(\mathbf{X})$
 * and store result in $vec(\mathbf{Y})$.
 * \param[in] A A $nm \times nm$ sparse marix.
 * \param[in] X A $n \times m$ "grid" matrix.
 * \param[out] Y A "grid" matrix s.t.
 * $vec(\mathbf{Y}) := \mathbf{A} * vec(\mathbf{X})$.
 */
void mult(const SparseMatrix<double> & A,
          const MatrixXd & X, MatrixXd & Y) {
    // Size check
    assert( A.rows() == X.rows()*X.cols() &&
            A.cols() == X.rows()*X.cols() &&
            "Inconsistent size of A!");
    // Ensure correct size of $\mathbf{Y}$.
    Y.resizeLike(X);

    // The second "const" here is needed to "promis" X will not be modified.
    const Map<const VectorXd> x(X.data(), X.size());
    // Notice how "Map" can also be written onto.
    Map<VectorXd>(Y.data(), X.size()) = A * x;
}

/* \brief Solve $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$ for $\mathbf{X}$.
 * Solve the linear system given by the sparse matrix $\mathbf{A}$
 * and r.h.s $\mathbf{Y}$.
 * \param[in] A A $nm \times nm$ sparse marix.
 * \param[in] Y A $n \times m$ "grid" matrix.
 * \param[out] X A "grid" matrix s.t.
 * $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$.
 */
void solve(const SparseMatrix<double> & A,
           const MatrixXd & Y, MatrixXd & X) {
    // Size check
    assert( A.rows() == Y.rows()*Y.cols() &&
            A.cols() == Y.rows()*Y.cols() &&
            "Inconsistent size of A!");
    // Ensure correct size of $\mathbf{X}$.
    X.resizeLike(Y);

    // Sparse LU decomposition object
    SparseLU<SparseMatrix<double>> solver;

    // Prepare LU factorization
    solver.analyzePattern(A);
    solver.factorize(A);

    // Check if factorization failed.
    if( solver.info() != Success ) {
        std::cerr << "Fatal: LU decomposition of A failed!!!"
                  << std::endl
                  << "Is the matrix invertible?"
                  << std::endl;
        return;
    }

    // Solve system after factorization.
    // Notice how you can assign to a Map object.
    Map<VectorXd>(X.data(), Y.size()) = solver.solve(
                Map<const VectorXd>(Y.data(), Y.size())
                );
}

int main(int argc, char** argv) {
    // Size of the "grid" matrices $\mathbf{X}$ and $\mathbf{Y}$.
    index_t n = 100, m = 100;
    // If provided, $n$ and $m$ will be read from command line.
    if(argc > 1) {
        n = std::stoi(argv[1]);
    }
    if(argc > 2) {
        m = std::stoi(argv[2]);
    }

    // The stencil matrix $\mathbf{S}$.
    Matrix3d S;
    S << 0, 1,  0,
         1, -4, 1,
         0, 1,  0;

    // Construct $\mathbf{A}$.
    std::cout << "--> Building matrix A..."
              << std::endl;
    SparseMatrix<double> A = build_matrix(S, shape_t{n, m});

    // The grid matrices $\mathbf{X}$ and $\mathbf{Y}$.
    MatrixXd X(n,m), Y;

    // A function of indices $(i,j)$ of the matrix $\mathbf{X}$.
    auto f = [n,m] (index_t i , index_t j) {
        if(i > n/4 && i < n*3/4 && j > m/4 && j < m*3/4 ) return 1.;
        else return 0.;
    };

    // Fill entries of  $\mathbf{X}$ given $f$.
    std::cout << "--> Evaluating f at the indices of X..."
              << std::endl;
    eval(X, f);

    // Compute $\mathbf{A} vec(\mathbf{X})$.
    std::cout << "--> Computing A*vec(X)..."
              << std::endl;
    mult(A, X, Y);

    // Compute $\mathbf{A}^{-1} vec(\mathbf{X})$.
    std::cout << "--> Solving A*vec(X) = vec(Y)..."
              << std::endl;
    solve(A, X, Y);

    // Plot values of $\mathbf{X}$.
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
