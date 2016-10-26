#include <algorithm>
#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

using namespace Eigen;

/* @brief Factorize matrix $X$ as $X = AB'$
 * @param[in] X An $m \times n$ matrix of rank at most $k$
 * @param[in] k Rank of matrix $X$
 * @param[out] A The $m \times k$ matrix from the decomposition of $X$
 * @param[out] B The $n \times k$ matrix from the decomposition of $X$
 */
/* SAM_LISTING_BEGIN_0 */
void factorize_X_AB(const MatrixXd & X, size_t k, MatrixXd & A, MatrixXd & B) {
	
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
	if(k+1 <= std::min(m,n) && s(k) > tol) {
		std::cerr << "The rank of the matrix is greater than the required rank" << std::endl;
	}
	
	A = U.leftCols(std::min(k,(size_t)U.cols())) * s.head(std::min(k,(size_t)s.size()));
	B = V.leftCols(std::min(k,(size_t)V.cols()));
}
/* SAM_LISTING_END_0 */

/* @brief SVD decomposition and ...
 * @param[in] A An $m \times k$ matrix
 * @param[in] B An $n \times k$ matrix
 * @param[out] U The $n \times k$ matrix from ...
 * @param[out] s The $k$-dimensional vector of singular values of ...
 * @param[out] V The $n \times k$ matrix from ...
 */
/* SAM_LISTING_BEGIN_1 */
void svd_AB(const MatrixXd & A, const MatrixXd & B,
MatrixXd & U, VectorXd & s, MatrixXd & V) {
	
	assert(A.cols() == B.cols()
           && "Matrices A and B should have the same column number");
	size_t n = A.rows();
	size_t k = A.cols();
	
	// QA: m x k; QB: n x k; RA,RB k x k
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
	
	// U,V: k x k
	JacobiSVD<MatrixXd> svd(RA*RB.transpose(), ComputeThinU | ComputeThinV);
	s = svd.singularValues();
	U = svd.matrixU();
	V = svd.matrixV();
	
	// U: m x k; V: n x k
	U = QA*U;
	V = QB*V;
}
/* SAM_LISTING_END_1 */

/* @brief Find $Az$ and $Bz$ s.t.
 * $Ax*Bz'$ is the best approximation of $AxBx'+AyBy'$ with rank $k$
 * @param[in] Ax An $m \times k$ matrix
 * @param[in] Ay An $m \times k$ matrix
 * @param[in] Bx An $n \times k$ matrix
 * @param[in] By An $n \times k$ matrix
 * @param[out] Az The $m \times k$ matrix to form $Az*Bz'$
 * @param[out] Bz The $n \times k$ matrix to form $Az*Bz'$
 */
/* SAM_LISTING_BEGIN_2 */
void rank_k_approx(const MatrixXd & Ax, const MatrixXd & Ay,
const MatrixXd & Bx, const MatrixXd & By,
MatrixXd & Az, MatrixXd & Bz) {
	
	assert(Ax.rows() == Ay.rows() && Ax.cols() == Ay.cols()
           && "Matrices Ax and Ay should have the same dimensions");
	assert(Bx.rows() == By.rows() && Bx.cols() == By.cols()
           && "Matrices Bx and By should have the same dimensions");
	
	MatrixXd A(Ax.rows(), Ax.cols()+Ay.cols()); A << Ax, Ay;
	MatrixXd B(Bx.rows(), Bx.cols()+By.cols()); B << Bx, By;
	MatrixXd U, V; VectorXd s;
	svd_AB(A, B, U, s, V);
	// U: m x 2k; S: 2k x 2k; V: n x 2k
	
	size_t k = Ax.cols();
	Az = U.leftCols(k) * s;
	Bz = V.leftCols(k);
}
/* SAM_LISTING_END_2 */

int main() {
	size_t m = 3;
	size_t n = 2;
	size_t k = 2;
    MatrixXd X(m,n);
    X << 5, 0, 2, 1, 7, 4;
    
    MatrixXd A, B;
    factorize_X_AB(X, k, A, B);
    
    std::cout << "A =" << std::endl << A << std::endl;
    std::cout << "B =" << std::endl << B << std::endl;
    
    A.resize(m,k); B.resize(n,k);
	A << 2, 1, 2, 3, 6, 1;
	B << 4, 4, 5, 0;
    MatrixXd U, V; VectorXd s;
    
    svd_AB(A, B, U, s, V);
    
    std::cout << "U =" << std::endl << U << std::endl;
    std::cout << "s =" << std::endl << s << std::endl;
    std::cout << "V =" << std::endl << V << std::endl;
    
    MatrixXd Ax(m,k), Ay(m,k), Bx(n,k), By(n,k), Az, Bz;
	Ax << 1,  0, 9, 2, 6, 3;
	Ay << 8, -2, 3, 4, 5, 8;
	Bx << 2, 1, 2, 3;
	By << 4, 4, 5, 0;
	
	rank_k_approx(Ax, Ay, Bx, By, Az, Bz);
	
	std::cout << "Az =" << std::endl << Az << std::endl;
    std::cout << "Bz =" << std::endl << Bz << std::endl;
}
