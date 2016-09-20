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


    return Q;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
int main(void) {
    // Orthonormality test
    unsigned int n = 9;
    MatrixXd A = MatrixXd::Random(n,n);
    // TODO: use gramschmidt to compute orthonormalization of
    // the matrix $\mathbf{A}$.
}
/* SAM_LISTING_END_2 */
