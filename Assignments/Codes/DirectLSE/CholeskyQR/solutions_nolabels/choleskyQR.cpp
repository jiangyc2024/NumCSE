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
	Q = R.transpose().triangularView<Lower>().solve(A.transpose()).transpose();
	// .triangularView() template member only accesses the triangular part
	// of a dense matrix and allows to easily solve linear problem
}

/* @brief Direct QR decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void DirectQR(const MatrixXd & A, MatrixXd & R, MatrixXd & Q) {
	
	size_t m = A.rows();
	size_t n = A.cols();
	
	HouseholderQR<MatrixXd> QR = A.householderQr();
    Q = QR.householderQ() * MatrixXd::Identity(m, std::min(m, n));
    R = MatrixXd::Identity(std::min(m, n), m) * QR.matrixQR().triangularView<Upper>();
    // If A: m x n, then Q: m x m and R: m x n.
    // If m > n, however, the extra columns of Q and extra rows of R are not needed.
    // Matlab returns this "economy-size" format calling "qr(A,0)",
    // which does not compute these extra entries.
    // With the code above, Eigen is smart enough to not compute the discarded vectors.
}

int main() {
	size_t m = 3;
	size_t n = 2;
    MatrixXd A(m,n);
    double epsilon = std::numeric_limits<double>::denorm_min();
    A << 3, 5, 1, 9, 7, 1;
  //A << 1, 1, 0.5*epsilon, 0, 0, 0.5*epsilon;
    std::cout << "A =" << std::endl << A << std::endl;
    
    MatrixXd R, Q;
    CholeskyQR(A, R, Q);
	
	std::cout << "From Cholesky: R =" << std::endl << R << std::endl;
    std::cout << "From Cholesky: Q =" << std::endl << Q << std::endl;
    
    DirectQR(A, R, Q);
    
	std::cout << "Direct QR: R =" << std::endl << R << std::endl;
    std::cout << "Direct QR: Q =" << std::endl << Q << std::endl;
}
