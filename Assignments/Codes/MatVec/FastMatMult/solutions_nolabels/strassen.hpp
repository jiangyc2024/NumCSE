#include <Eigen/Dense>

using namespace Eigen;

/* \brief Compute the Matrix product $A \times B$ using Strassen's algorithm.
 * \param[in] A Matrix $2^k \times 2^k$
 * \param[in] B Matrix $2^k \times 2^k$
 * \param[out] Matrix product of A and B of dim $2^k \times 2^k$
 */
/* SAM_LISTING_BEGIN_1 */
MatrixXd strassenMatMult(const MatrixXd& A, const MatrixXd& B) {
#if SOLUTION
    // Ensure square matrix
    assert(A.rows() == A.cols() && "Matrix A must be square");
    assert(B.rows() == B.cols() && "Matrix B must be square");
    // matrix dimension must be a power of 2
    assert(A.rows() % 2 == 0 && "Matrix dimensions must be a power of two.");

    const unsigned n = A.rows();

    // The function is recursive and acto on bloks of size $n/2 \times n/2$
    // i.e. exploits fast product of 2x2 block matrix
    if ( n==2 ) { // End of recursion
        MatrixXd C(2, 2);
        C << A(0,0)*B(0,0) + A(0,1)*B(1,0),
             A(0,0)*B(0,1) + A(0,1)*B(1,1),
             A(1,0)*B(0,0) + A(1,1)*B(1,0),
             A(1,0)*B(0,1) + A(1,1)*B(1,1);
        return C;
    }

    MatrixXd Q0(n/2, n/2), Q1(n/2, n/2), Q2(n/2, n/2), Q3(n/2, n/2),
             Q4(n/2, n/2), Q5(n/2, n/2), Q6(n/2, n/2);

    MatrixXd A11 = A.topLeftCorner(n/2, n/2);
    MatrixXd A12 = A.topRightCorner(n/2, n/2);
    MatrixXd A21 = A.bottomLeftCorner(n/2, n/2);
    MatrixXd A22 = A.bottomRightCorner(n/2, n/2);

    MatrixXd B11 = B.topLeftCorner(n/2, n/2);
    MatrixXd B12 = B.topRightCorner(n/2, n/2);
    MatrixXd B21 = B.bottomLeftCorner(n/2, n/2);
    MatrixXd B22 = B.bottomRightCorner(n/2, n/2);

    Q0 = strassenMatMult(A11 + A22, B11 + B22);
    Q1 = strassenMatMult(A21 + A22, B11);
    Q2 = strassenMatMult(A11,       B12 - B22);
    Q3 = strassenMatMult(A22,       B21 - B11);
    Q4 = strassenMatMult(A11 + A12, B22);
    Q5 = strassenMatMult(A21 - A11, B11 + B12);
    Q6 = strassenMatMult(A12 - A22, B21 + B22);

    MatrixXd C(n,n);
    C << Q0 + Q3 - Q4 + Q6,
         Q2 + Q4,
         Q1 + Q3,
         Q0 + Q2 - Q1 + Q5;

    return C;
#else
    const unsigned n = A.rows();

    MatrixXd C(n,n);

    // TODO: implement Strassen's algorithm to multiply
    // the two matrices $A$ and $B$ and return the result

    return C;
#endif
}
/* SAM_LISTING_END_1 */
