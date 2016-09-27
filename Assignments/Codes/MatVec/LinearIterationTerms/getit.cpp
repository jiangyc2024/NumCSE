#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;

#if SOLUTION
/* \brief Performs the computation $y = A^k x$.
 */
#endif // SOLUTION
/* SAM_LISTING_BEGIN_1 */
VectorXcd getit(const MatrixXd &A, const VectorXd &x, unsigned int k) {

# if SOLUTION
    // As said in the problem formulation, we may assume that the
    // following three lines have complexity $O(n^3)$.
#endif // SOLTUION
    EigenSolver<MatrixXd> eig = EigenSolver<MatrixXd>(A);
    const VectorXcd & V = eig.eigenvalues();
    const MatrixXcd & W = eig.eigenvectors();

# if SOLUTION
    // The following operation requires only a loop over
    // the dimension of $cx$, which is $n$.
#endif // SOLTUION
    VectorXcd cx = x.cast<std::complex<double>>();

# if SOLUTION
    // The first operator* is a matrix vector multiplication
    // with complexity $O(n^2)$.
#endif // SOLTUION
    VectorXcd ret = W *
# if SOLUTION
        // The componentwise power has complexity $O(n^2)$.
        // The second operator* is a vector-vector componentwise
        // multiplication, with complexity $O(n)$.
#endif // SOLTUION
        (
          V.array().pow(k) *
# if SOLUTION
          // In the following line, a linear system is solved, operation
          // with complexity $O(n^3)$
#endif // SOLTUION
          (W.partialPivLu().solve(cx)).array()
        );

    return ret;
}
/* SAM_LISTING_END_1 */

int main() {
    // Some arbitrary data to test getit
    MatrixXd A(4,4);
    A << 1,  2,  3,  4,
         5,  6,  7,  8,
         9,  10, 11, 12,
         13, 14, 15, 16;
    VectorXd x(4);
    x << 4,  5,  6,  7;
    unsigned int  k = 9;

    // Testing the implementation with some matrix
    VectorXcd yg = getit(A, x, k);
    std::cout << "getit(A,x, k) = " << std::endl
              << yg << std::endl;

#if SOLUTION
    // Checking that getit works
    VectorXcd yp = A.pow(k)*x;
    std::cout << "A^k x = " << std::endl
              << yp << std::endl;
    double err = (yg - yp).norm();
#endif // SOLUTION
}
