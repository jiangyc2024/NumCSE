#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/LU>

#include "timer.h"

using namespace Eigen;

/* @brief Circular shift (downwards) of $b$
 * @param[in,out] b The input $n$-dimensional vector shifted downwards
 */
void shift(VectorXd & b) {
    int n = b.size();

    double temp = b(n-1);
    for(int k = n-2; k >= 0; --k) {
        b(k+1) = b(k);
    }
    b(0) = temp;
}

/* @brief Compute $X = A^{-1}*[b_1,...,b_n],\; b_i = i$-th cyclic shift of $b$.
 * Function with naive implementation.
 * @param[in] A An $n \times n$ matrix
 * @param[in] b An $n$-dimensional vector
 * @param[out] X The $n \times n$ matrix $X = A^{-1}*[b_1,...,b_n]$
 */
/* SAM_LISTING_BEGIN_0 */
void solvpermb(const MatrixXd & A, VectorXd & b, MatrixXd & X) {
    // Size of b, which is the size of A
    int n = b.size();
    assert( n == A.rows() && n == A.cols()
            && "Error: size mismatch!");
    X.resize(n,n);

    // For each loop iteration:
    // 1. solve the linear system $Ax = b$,
    // 2. store the result in a column of $X$
    // 3. and shift $b$ by one element for the next iteration.
    for(int l = 0; l < n; ++l) {
        X.col(l) = A.fullPivLu().solve(b);

        shift(b);
    }
}
/* SAM_LISTING_END_0 */

/* @brief Compute $X = A^{-1}*[b_1,...,b_n],\; b_i = i$-th cyclic shift of $b$,
 * Function has complexity $O(n^3)$
 * @param[in] A An $n \times n$ matrix
 * @param[in] b An $n$-dimensional vector
 * @param[out] X The $n \times n$ matrix $X = A^{-1}*[b_1,...,b_n]$
 */
/* SAM_LISTING_BEGIN_1 */
void solvpermb_on3(const MatrixXd & A, VectorXd & b, MatrixXd & X) {
    // Size of b, which is the size of A
    int n = b.size();
    assert( n == A.cols() && n == A.rows()
            && "Error: size mismatch!");
    X.resize(n,n);

    // Notice that here we reuse the LU factorization
    FullPivLU<MatrixXd> LU = A.fullPivLu();

    for(int l = 0; l < n; ++l) {
        X.col(l) = LU.solve(b);

        shift(b);
    }
}
/* SAM_LISTING_END_1 */

int main() {
    unsigned int n = 9;
    // Compute with both solvers
    std::cout << "--> Check that the solvers are correct" << std::endl;
    MatrixXd A = MatrixXd::Random(n,n);
    VectorXd b = VectorXd::Random(n);
    MatrixXd Xi, Xr, X;

    // std::cout << "b = " << std::endl
    // << b << std::endl;

    solvpermb(A,b,Xi);
    // std::cout << "Direct porting from MATLAB (naive solver): "
    // << std::endl << X << std::endl;
    // std::cout << "A*X = " << std::endl
    // << A*X << std::endl;

    solvpermb_on3(A,b,Xr);
    // std::cout << "Reusing LU: " << std::endl
    // << X << std::endl;
    // std::cout << "A*X = " << std::endl
    // << A*X << std::endl;

    std::cout << "Error = " << (Xi - Xr).norm() << std::endl;

    // Compute runtimes of different solvers
    std::cout << "--> Runtime comparison of naive solver vs reusing LU" << std::endl;
    unsigned int repeats = 3;

    // Header
    std::cout << std::setw(20) << "n"
              << std::setw(20) << "time no reuse [s]"
              << std::setw(20) << "time reuse [s]"
              << std::endl;

    // Loop over matrix size
    for(unsigned int p = 2; p <= 7; ++p) {
        // Timers
        Timer tm_naive, tm_reuseLU;
        unsigned int n = pow(2,p);

        // Repeat test many times
        for(unsigned int r = 0; r < repeats; ++r) {
            A = MatrixXd::Random(n,n);
            b = VectorXd::Random(n);

            // Compute runtime with naive solver
            tm_naive.start();
            solvpermb(A,b,X);
            tm_naive.stop();
            // Compute runtime reusing LU factorisation
            tm_reuseLU.start();
            solvpermb_on3(A,b,X);
            tm_reuseLU.stop();
        }
        // Print runtimes
        std::cout << std::setw(20) << n
                  << std::scientific << std::setprecision(3)
                  << std::setw(20) << tm_naive.min()
                  << std::setw(20) << tm_reuseLU.min()
                  << std::endl;
    }
}
