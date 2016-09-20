#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

/* @brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Ustable GS algorithm: output is prone to cancellation issues.
 * @param[in] $\mathbf{A}$ Matrix of linearly independent columns
 * @return Matrix with ONB of $span(a_1, \cdots, a_n)$ as columns
 */
/* SAM_LISTING_BEGIN_1 */
MatrixXd gram_schmidt(const MatrixXd & A) {
    // We create a matrix Q with the same size as A
    // and copied from A
    MatrixXd Q(A);

    // The first vector just gets normalized
    Q.col(0).normalize();

    for(unsigned int j = 1; j < A.cols(); ++j) {
        // Replace inner loop over each previous vector in Q with fast
        // matrix-vector multiplication
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));

        // Normalize vector, if possible
        // (otherwise it means colums of $\mathbf{A}$ are
        // almost linear dependant)
        double eps = std::numeric_limits<double>::denorm_min();
        if( Q.col(j).norm() <= eps * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because "
                      << "A has (almost) linear dependant "
                      << "columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }

    return Q;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
int main(void) {
    // Orthonormality test
    unsigned int n = 9;
    MatrixXd A = MatrixXd::Random(n,n);
    MatrixXd Q = gram_schmidt( A );

    // Compute how far is $\mathbf{Q}*\mathbf{Q}^\top$ from the identity
    // i.e. "How far is Q from being orthonormal?"
    double err = (Q*Q.transpose() - MatrixXd::Identity(n,n))
        .norm();

    // Error has to be small, but not zero (why?)
    std::cout << "Error is: "
              << err
              << std::endl;

    // If error is too big, we exit with error
    double eps = std::numeric_limits<double>::denorm_min();
    exit(err < eps);
}
/* SAM_LISTING_END_2 */
