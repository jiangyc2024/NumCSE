#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

using namespace Eigen;

/* @brief factorize X mxn of rank at most k in X=AB';  A mxk, Bnxk
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
/* SAM_LISTING_BEGIN_0 */
void factorizeXAB(const MatrixXd & X, size_t k, MatrixXd & A, MatrixXd & B) {
	
	size_t m = X.rows();
	size_t n = X.cols();
	double tol = 1e-6;
	
	JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
	VectorXd s = svd.singularValues();
	MatrixXd U = svd.matrixU();
	MatrixXd V = svd.matrixV();
	
	// You can ask for only thin $U$ or $V$ to be computed.
	// In case of a rectangular $n \times p$ matrix,
	// letting $m$ be the smaller value among $n$ and $p$,
	// there are only $m$ singular vectors;
	// the remaining columns of $U$ and $V$ do not correspond to actual singular vectors
	// and are not returned in the thin format.
	if(k+1 <= min(m,n) && s(k) > tol) {
		std::cerr << "The rank of the matrix is greater than the required rank" << std::endl;
	}
	
	A = U.leftCols(min(k,U.cols())) * s.head(min(k,s.size()));
	B = V.leftCols(min(k,V.cols()));
}
/* SAM_LISTING_END_0 */

/* @brief Direct QR decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
/* SAM_LISTING_BEGIN_1 */
void svdAB(const MatrixXd & A, const MatrixXd & B,
MatrixXd & U, MatrixXd & S, MatrixXd & V) {
	
	assert(A.rows() == B.rows() && A.cols() == B.cols()
           && "Matrices A and B should have the same dimensions");
	size_t n = A.rows();
	size_t k = A.cols();
	
	// QA,QB: n x k; RA,RB k x k
	HouseholderQR<MatrixXd> QRA = A.householderQr();
    MatrixXd QA = QRA.householderQ();
    MatrixXd RA = QRA.matrixQR().triangularView<Upper>();
    HouseholderQR<MatrixXd> QRB = B.householderQr();
    MatrixXd QB = QRB.householderQ();
    MatrixXd RB = QRB.matrixQR().triangularView<Upper>();
    
    // To return the same "economy-size decomposition" as Matlab
    if(n > k) {
		QA = QA.leftCols(k);
		RA = RA.topRows(k);
		QB = QB.leftCols(k);
		RB = RB.topRows(k);
	}
	
	// U,S,V: k x k
	JacobiSVD<MatrixXd> svd(RA*RB.transpose(), ComputeThinU | ComputeThinV);
	VectorXd s = svd.singularValues();
	S.setZeros(k,k);
	S.diagonal() = s;
	U = svd.matrixU();
	V = svd.matrixV();
	
	// U,V:   n x k
	U = QA*U;
	V = QB*V;
}
/* SAM_LISTING_END_1 */

int main() {
	size_t m = 3;
	size_t n = 2;
    MatrixXd A(m,n);
    double epsilon = std::numeric_limits<double>::denorm_min();
    A << 1, 1, 0.5*epsilon, 0, 0, 0.5*epsilon;
    std::cout << "A =" << std::endl << A << std::endl;
    
    MatrixXd R, Q;
    CholeskyQR(A, R, Q);
	
	std::cout << "From Cholesky: R =" << std::endl << R << std::endl;
    std::cout << "From Cholesky: Q =" << std::endl << Q << std::endl;
    
    DirectQR(A, R, Q);
    
	std::cout << "Direct QR: R =" << std::endl << R << std::endl;
    std::cout << "Direct QR: Q =" << std::endl << Q << std::endl;
}
