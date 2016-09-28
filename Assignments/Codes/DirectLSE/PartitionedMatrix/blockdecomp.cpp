#include <Eigen/Dense>
#include <iostream>

#include "timer.h"

using namespace Eigen;

/* \brief Use efficient implementation A*x = bb
 * \param[in] R MatrixXd is nxn and upper triangular
 * \param[in] v VectorXd is nx1
 * \param[in] u VectorXd is nx1
 * \param[in] bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
 * \param[out] x solution A*bb = x
 */
/* SAM_LISTING_BEGIN_1 */
void solvelse(const MatrixXd & R,
              const VectorXd & v, const VectorXd & u,
              const VectorXd & bb, VectorXd & x) {
#if SOLUTION
    // Size of R, which is wize of u, v, and size of bb is n+1
    unsigned int n = R.rows();
    assert( n == R.cols() && n == u.size()
            && n == v.size() && n+1  == bb.size()
            && "Size mismatch!");

    TriangularView<MatrixXd, Upper> triR = R.triangularView<Upper>();

    // $s$ is the Schur complement and, in this case, is a scalar
    // also b_s is a scalar
    // $snv = s^{-1}$, $b_s$ as in lecture notes
    // $sinvbs = s^{-1}*b_s$
    double sinv = - 1. / u.dot(triR.solve(v));
    double bs = (bb(n) - u.dot(triR.solve(bb.head(n))));
    double sinvbs = sinv*bs;

    // Stack the vector (z, \xi)^T =: x
    x = VectorXd::Zero(n+1);
    x << triR.solve(bb.head(n) - v*sinvbs),
         sinvbs;
#else // TEMPLATE
    // TODO: Implement a block-gauss elimination
#endif // TEMPLATE
}
/* SAM_LISTING_END_1 */

//! \brief Use Eigen's LU-solver to solve Ax = y
//! \param[in] R MatrixXd is nxn and upper triangular
//! \param[in] v VectorXd is nx1
//! \param[in] u VectorXd is nx1
//! \param[in] bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
//! \param[out] x solution A*bb = x
/* SAM_LISTING_BEGIN_2 */
void solvelse_lu(const MatrixXd & R,
                 const VectorXd & v, const VectorXd & u,
                 const VectorXd & bb, VectorXd & x) {
#if SOLUTION
    // Size of R, which is wize of u, v, and size of bb is n+1
    unsigned int n = R.rows();
    assert( n == R.cols() && n == u.size()
            && n == v.size() && n+1  == bb.size()
            && "Size mismatch!");

    MatrixXd A(n+1,n+1);
    // Build matrix using blocks
    A << R,             v,
         u.transpose(), 0;

    // Solve system using LU decpomposition
    x = A.partialPivLu().solve(bb);
#else // TEMPLATE
    // TODO: solve the system using LU decomposition
#endif
}
/* SAM_LISTING_END_2 */

int main() {

    ///// CORRECTNESS TESTS

    // Bunch of random vectors/matrices
    int n = 9;

    // Random test vectors
    Eigen::MatrixXd R = Eigen::MatrixXd::Random(n,n).triangularView<Eigen::Upper>();
    Eigen::VectorXd v = Eigen::VectorXd::Random(n);
    Eigen::VectorXd u = Eigen::VectorXd::Random(n);
    Eigen::VectorXd bb = Eigen::VectorXd::Random(n+1);
    Eigen::VectorXd xe, xo;

    // Check that answer is correct
    std::cout << "--> Check correctness." << std::endl;
    solvelse(R, v, u, bb, xo);
    std::cout << "Block Gauss:" << std::endl
              << xo << std::endl;
    // Compare with eigen partial pivot LU-solve
    solvelse_lu(R, v, u, bb, xe);
    std::cout << "Eigen LU:" << std::endl
              << xe << std::endl;
    std::cout << "Error:" << std::endl
              << (xe-xo).norm() << std::endl;

#if SOLUTION
    ///// RUNTIME COMPARISONS

    // Compute runtime of different implementations of kron
    std::cout << "--> Runtime comparisons." << std::endl;
    unsigned int repeats = 3;

    // Header
    std::cout << std::setw(15) << "M"
              << std::setw(15) << "own"
              << std::setw(15) << "eigen"
              << std::endl;
    for(unsigned int p = 2; p <= 10; p++) {
        unsigned int M = pow(2,p);
        Timer tm_own, tm_eigen_lu;

        for(unsigned int r = 0; r < repeats; ++r) {
            // We will use only upper triangular part
            R = MatrixXd::Random(M,M).triangularView<Upper>();
            v = VectorXd::Random(M);
            u = VectorXd::Random(M);
            bb = VectorXd::Random(M+1);

            // Time of own implementation
            tm_own.start();
            solvelse(R, v, u, bb, x);
            tm_own.stop();

            // Time of eigen
            tm_eigen_lu.start();
            solvelse_lu(R, v, u, bb, x);
            tm_eigen_lu.stop();
        }

        // Print table
        std::cout << std::setw(15) << M
                  << std::setw(15) << tm_own.min() / 1000000.
                  << std::setw(15) << tm_eigen_lu.min() / 1000000.
                  << std::endl;
    }
#endif // SOLUTION
}
