#include <Eigen/Dense>

using namespace Eigen;

/** 
 * @brief Compute the Matrix product $A \times B$ using Strassen's algorithm.
 * @param[in] A Matrix $2^k \times 2^k$
 * @param[in] B Matrix $2^k \times 2^k$
 * @param[out] Matrix product of A and B of dim $2^k \times 2^k$
 */
/* SAM_LISTING_BEGIN_1 */
MatrixXd strassenMatMult(const MatrixXd& A, const MatrixXd& B) {
    const unsigned n = A.rows();

    MatrixXd C(n,n);
    // TODO: implement strassen's algorithm to multiply the two matrices $A$ and $B$ and return the result
    return C;
}
/* SAM_LISTING_END_1 */
