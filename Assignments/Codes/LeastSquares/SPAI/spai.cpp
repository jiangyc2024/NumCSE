#include <iostream>
#include <cmath>
#include <ctime>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/KroneckerProduct>

#include "timer.h"

using namespace Eigen;

using index_t = int;

/* @brief Compute $B = \argmin_{X \in P(A)} |I-AX|_F$
 * @param[in] A An $n \times n$ matrix
 * @param[out] B The $n \times n$ matrix $= \argmin_{X \in P(A)} |I-AX|_F$
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> spai(SparseMatrix<double> & A) {
    // Size check
    assert(A.rows() == A.cols() &&
    "Matrix must be square!");
    unsigned int N = A.rows();

    // Needed to make sure ___Ptr functions return
    // arrays specified in CRS format
    A.makeCompressed();

    // Obtain pointers to data of A
    double* valPtr = A.valuePtr();
    index_t* innPtr = A.innerIndexPtr();
    index_t* outPtr = A.outerIndexPtr();

    // Create vector for triplets of B and reserve enough space
    std::vector<Triplet<double>> triplets;
    triplets.reserve(A.nonZeros());

#if SOLUTION
    // Loop over each column of $A$
    for(unsigned int i = 0; i < N; ++i) {
        // Number of non zeros in current column of $A$
        index_t nnz_i = outPtr[i+1] - outPtr[i];
        // If column is empty: skip column (matrix is not invertible)
        if(nnz_i == 0) continue;

        // Temporarily create a (small, dense) matrix on which normal formula
        // will be applied. We project the space $\mathcal{P}(A)$
        // onto $\mathcal{P}(a_i)$
         MatrixXd C = MatrixXd::Zero(N, nnz_i);
//        SparseMatrix<double> C(N, nnz_i);
//        std::vector<Triplet<double>> C_triplets;
//        C_triplets.reserve(nnz_i*nnz_i);

        // Need to build matrix $C$. To this end we remove all columns
        // from $A$ for which $a_i = 0$.
        // To this end: loop over all non-zero entries of the $i$-th column
        // This loop has length $n$
        for(unsigned int k = outPtr[i]; k < outPtr[i+1]; ++k) {
            // Row of this non-zero entry
            index_t row_k = innPtr[k];
            // Number of non-zero entries for $row_k$-th column
            index_t nnz_k = outPtr[row_k+1] - outPtr[row_k];
            // Loop over all non-zeros of $row_k$th-column
            // This loop has length complexity $n$
            for(unsigned int l = 0; l < nnz_k; ++l) {
                unsigned int innIdx = outPtr[row_k] + l;
                 C(innPtr[innIdx], k - outPtr[i]) = valPtr[innIdx];
//                C_triplets.emplace_back(Triplet<double>(innPtr[innIdx], k - outPtr[i], valPtr[innIdx]));
            }
        }
//        C.setFromTriplets(C_triplets.begin(), C_triplets.end());
//        C.makeCompressed();
//        SparseMatrix<double> S = C.transpose() * C;
//        MatrixXd M = MatrixXd(S);
//        VectorXd xt = C.row(i).transpose();

        // Compute $C^\top C$ and solve normal equation.
        // Complexity of product:
        // Complexity of solve:
        // Size of $b$ is at most $n$.
//        VectorXd b = M.partialPivLu().solve(xt);
        VectorXd b = (C.transpose() * C).partialPivLu().solve(C.row(i).transpose());
        // Loop as length $n$.
        for(unsigned int k = 0; k < b.size(); ++k) {
            triplets.emplace_back(Triplet<double>(innPtr[outPtr[i] + k], i, b(k)));
        }
    }
#else // TEMPLATE
    // TODO: build $B$ in triplet format, exploiting the sparse format of $A$
#endif // TEMPLATE

    // Build and return SPAI preconditioner
    SparseMatrix<double> B = SparseMatrix<double>(N,N);
    B.setFromTriplets(triplets.begin(), triplets.end());
    B.makeCompressed();
    return B;
}
/* SAM_LISTING_END_1 */

// Run conditionally
const bool small_test = false;
const bool big_test = true;

int main(int argc, char **argv) {
    unsigned int n = 60;
    if(argc > 1) {
        n = std::stoi(argv[1]);
    }
    srand(time(NULL));

    if(small_test)
    {
        SparseMatrix<double> M(5,5);
        M.coeffRef(3,4) = 1;
        M.coeffRef(4,4) = 2;
        M.coeffRef(1,4) = 3;
        M.coeffRef(3,3) = 4;
        M.coeffRef(3,2) = 4;
        M.coeffRef(2,3) = 4;
        M.coeffRef(2,2) = 5;
        M.coeffRef(3,1) = 6;
        M.coeffRef(0,0) = 9;

        SparseMatrix<double> N = spai(M);
        SparseMatrix<double>I(n,n);
        I.setIdentity();

        std::cout << "M = " << std::endl
                  << M
                  << std::endl;
        std::cout << "N = " << std::endl
                  << N
                  << std::endl;

        std::cout << "Error:"
                  << (I - M*N).norm()
                  << std::endl;

    }
    if(big_test)
    {
        SparseMatrix<double> M(n*n,n*n);

        SparseMatrix<double>I(n,n);
        I.setIdentity();
        MatrixXd R = MatrixXd::Random(n,n);
        M = kroneckerProduct(R, I);

        Timer tm;
        tm.start();
        SparseMatrix<double> N = spai(M);
        tm.stop();

        SparseMatrix<double>Ibig(n*n,n*n);
        Ibig.setIdentity();

        std::cout << "Error (n = " << n*n << "):"
                  << (Ibig - M*N).norm()
                  << std::endl
                  << "Elapsed:"
                  << tm.duration() << " s"
                  << std::endl;
    }
}
