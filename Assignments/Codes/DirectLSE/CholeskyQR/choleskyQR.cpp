#include <iostream>

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace Eigen;

/* @brief QR decomposition from Cholesky decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
/* SAM_LISTING_BEGIN_0 */
void CholeskyQR(const MatrixXd & A, MatrixXd & R, MatrixXd & Q) {
	
	MatrixXd AtA = A.transpose() * A;
	LLT<MatrixXd> L = AtA.llt();
	R = L.matrixL().transpose();
	Q = R.transpose().partialPivLu().solve(A.transpose()).transpose();
	// LLT only works with symmetric p.d. matrices,
	// but $R^\top$ is at least invertible...
}
/* SAM_LISTING_END_0 */

/* @brief Direct QR decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
/* SAM_LISTING_BEGIN_1 */
void DirectQR(const MatrixXd & A, MatrixXd & R, MatrixXd & Q) {
	
	HouseholderQR<MatrixXd> QR = A.householderQr();
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();
}
/* SAM_LISTING_END_1 */

int main() {
	size_t n = 3;
    MatrixXd A(n,n);
    A << 2, 9, 1, 4, 9, 1, 0, 4, 2;
    std::cout << "A =" << std::endl << A << std::endl;
    
    MatrixXd R, Q;
    CholeskyQR(A, R, Q);
	
	std::cout << "From Cholesky: R =" << std::endl << R << std::endl;
    std::cout << "From Cholesky: Q =" << std::endl << Q << std::endl;
    
    DirectQR(A, R, Q);
    
	std::cout << "Direct QR: R =" << std::endl << R << std::endl;
    std::cout << "Direct QR: Q =" << std::endl << Q << std::endl;
}
