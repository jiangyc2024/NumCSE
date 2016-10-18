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
    // TODO: convert grid index $(i,j)$ to index $I$ and return
    return 0;
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
    std::vector<Triplet<double>> triplets;

    // The matrix $\mathbf{A}$ will have size $N \times N$.
    // N = n*m
    const index_t N = size[0]*size[1];
    // TODO: build/compress/return sparse matrix $\mathbf{A}$ from
    // intermediate triplet format.
    return SparseMatrix<double>(N,N);
}
/* SAM_LISTING_END_1 */

/* \brief Evaluate function $f$ at the indices of $\mathbf{X}$.
 * Computes $(\mathbf{X})_{i,j} = f(i,j)$
 * ($+1$ added because C++-indices start at $0$).
 * \param[in,out] X Matrix $\mathbf{X}$, must have correct size.
 * \param[in] f A function $\mathbf{N}^2 \rightarrow \mathbf{R}$.
 */
/* SAM_LISTING_BEGIN_2 */
void eval(MatrixXd & X, std::function<double(index_t,index_t)> f) {
    // TODO: implement an evaluation function computing
    // $(\mathbf{X})_{i,j} = f(i,j)$.
}
/* SAM_LISTING_END_2 */

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

    // TODO: compute $vec(\mathbf{Y}) := \mathbf{A} * vec(\mathbf{X})$.
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

    // Sparse LU decomposition object
    SparseLU<SparseMatrix<double>> solver;

    // TODO: solve system $vec(\mathbf{Y}) = \mathbf{A} * vec(\mathbf{X})$.
}
/* SAM_LISTING_END_4 */

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
    /* SAM_LISTING_BEGIN_0 */
    auto f = [n,m] (index_t i , index_t j) {
        // TODO: implement f
        return 0.;
    };
    /* SAM_LISTING_END_0 */

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
