//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace Eigen;

/* @brief QR decomposition from Cholesky decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void CholeskyQR(const MatrixXd & A, MatrixXd & R, MatrixXd & Q) {
	
	MatrixXd AtA = A.transpose() * A;
	LLT<MatrixXd> L = AtA.llt();
	R = L.matrixL().transpose();
	Q = R.transpose().partialPivLu().solve(A.transpose()).transpose();
	// LLT only works with symmetric p.d. matrices,
	// but $R^\top$ is at least invertible...
}

/* @brief Direct QR decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void DirectQR(const MatrixXd & A, MatrixXd & R, MatrixXd & Q) {
	
	HouseholderQR<MatrixXd> QR = A.householderQr();
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();
}

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
