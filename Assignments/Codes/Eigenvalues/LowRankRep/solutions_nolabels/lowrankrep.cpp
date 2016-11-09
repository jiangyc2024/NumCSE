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
#include <Eigen/SVD>

using namespace Eigen;

/* @brief Factorize matrix $X$ into $X = AB'$
 * @param[in] X An $m \times n$ matrix of rank at most $k$
 * @param[in] k Rank of matrix $X$
 * @param[out] A The $m \times k$ matrix from the decomposition of $X$
 * @param[out] B The $n \times k$ matrix from the decomposition of $X$
 */
void factorize_X_AB(const MatrixXd & X, size_t k, MatrixXd & A,
						  MatrixXd & B) {
  size_t m = X.rows(), n = X.cols();
  double tol = 1e-6; // Tolerance for numerical rank: \lref{ex:svdrank}
  assert(k <= std::min(m,n)
	 && "Rank k cannot be larger than dimensions of X");
  
  JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
  // With Eigen::svd you can ask for thin $\VU$ or $\VV$ to be computed.
  // In case of a rectangular $m \times n$ matrix,
  // with $j$ the smaller value among $m$ and $n$,
  // there can only be at most $j$ singular values.
  // The remaining columns of $\VU$ and $\VV$ do not correspond
  // to actual singular vectors and are not computed in thin format.
  
  VectorXd s = svd.singularValues();
  MatrixXd S; S.setZero(s.size(),s.size());
  S.diagonal() = s;
  MatrixXd U = svd.matrixU();
  MatrixXd V = svd.matrixV();
  
  if(k+1 <= std::min(m,n) && s(k) > tol*s(0)) {
    // 1.\ condition checks if there is a singular value after $k$.
    // ($\VS$ cannot be larger than the dimensions of $\VX$.)
    // 2.\ condition checks if such singular value
    // (by definition non-negative) is "numerically" different from 0.
    std::cerr << "Rank of matrix X is greater than required rank k"
              << std::endl;
  }

  // Instead of multiplying with $\VS$ we could simply scale
  // the columns of $\VA$
  A = U.leftCols(std::min(k,(size_t)U.cols())) *
    S.topLeftCorner(std::min(k,(size_t)s.size()),
    std::min(k,(size_t)s.size()));
  B = V.leftCols(std::min(k,(size_t)V.cols()));
}

/* @brief Compute the SVD of $AB'$
 * @param[in] A An $m \times k$ matrix
 * @param[in] B An $n \times k$ matrix
 * @param[out] U The $n \times k$ matrix from the SVD of $AB'$
 * @param[out] S The $k \times k$ diagonal matrix of sing. vals of $AB'$
 * @param[out] V The $n \times k$ matrix from the SVD of $AB'$
 */
void svd_AB(const MatrixXd & A, const MatrixXd & B,
MatrixXd & U, MatrixXd & S, MatrixXd & V) {
	
	assert(A.cols() == B.cols()
           && "Matrices A and B should have the same column number");
	size_t m = A.rows();
	size_t n = B.rows();
	size_t k = A.cols();
	
	// QA: m x k; QB: n x k; RA,RB k x k
	HouseholderQR<MatrixXd> QRA = A.householderQr();
    MatrixXd QA = QRA.householderQ() * MatrixXd::Identity(m, std::min(m, k));
    MatrixXd RA = MatrixXd::Identity(std::min(m, k), m) * QRA.matrixQR().triangularView<Upper>();
    HouseholderQR<MatrixXd> QRB = B.householderQr();
    MatrixXd QB = QRB.householderQ() * MatrixXd::Identity(n, std::min(n, k));;
    MatrixXd RB = MatrixXd::Identity(std::min(n, k), n) * QRB.matrixQR().triangularView<Upper>();
    // To compute the same "economy-size" QR decomposition as Matlab -- see choleskyQR.cpp
	
	// U,V: k x k
	JacobiSVD<MatrixXd> svd(RA*RB.transpose(), ComputeFullU | ComputeFullV);
	// Thin matrices are unnecessary here as $RA*RB'$ is
	// a square $k \times k$ matrix!
	// Moreover, if you had computed the thin matrices $U$ and $V$
	// of the direct SVD of $AB'$,
	// you would still have dealt with the minimum between $m$ and $n$,
	// without exploiting $k << m, n$!
		
	VectorXd s = svd.singularValues();
	S.setZero(s.size(),s.size());
	S.diagonal() = s;
	U = svd.matrixU();
	V = svd.matrixV();
	
	// U: m x k; V: n x k
	U = QA*U;
	V = QB*V;
}

/* @brief Find $Az$ and $Bz$ such that
 * $Ax*Bz'$ is the best approximation of $AxBx'+AyBy'$ with rank $k$
 * @param[in] Ax An $m \times k$ matrix
 * @param[in] Ay An $m \times k$ matrix
 * @param[in] Bx An $n \times k$ matrix
 * @param[in] By An $n \times k$ matrix
 * @param[out] Az The $m \times k$ matrix to form $Az*Bz'$
 * @param[out] Bz The $n \times k$ matrix to form $Az*Bz'$
 */
void rank_k_approx(const MatrixXd & Ax, const MatrixXd & Ay,
const MatrixXd & Bx, const MatrixXd & By,
MatrixXd & Az, MatrixXd & Bz) {
	
	assert(Ax.rows() == Ay.rows() && Ax.cols() == Ay.cols()
           && "Matrices Ax and Ay should have the same dimensions");
	assert(Bx.rows() == By.rows() && Bx.cols() == By.cols()
           && "Matrices Bx and By should have the same dimensions");
	
	MatrixXd A(Ax.rows(), Ax.cols()+Ay.cols()); A << Ax, Ay;
	MatrixXd B(Bx.rows(), Bx.cols()+By.cols()); B << Bx, By;
	// U: m x 2k; S: 2k x 2k; V: n x 2k
	MatrixXd U, S, V;
	svd_AB(A, B, U, S, V);
	
	size_t k = Ax.cols();
	Az = U.leftCols(k) * S;
	Bz = V.leftCols(k);
}

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
    MatrixXd U, S, V;
    
    svd_AB(A, B, U, S, V);
    
    std::cout << "U =" << std::endl << U << std::endl;
    std::cout << "S =" << std::endl << S << std::endl;
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
