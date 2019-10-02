
#include "Sylvester.hpp"

using namespace Eigen;

// Kronecker product from kron.hpp.
void kron(const MatrixXd &A, const MatrixXd &B, MatrixXd &C) {
  // Allocate enough space for the matrix
  C = MatrixXd(A.rows() * B.rows(), A.cols() * B.cols());
  for (unsigned int i = 0; i < A.rows(); ++i) {
    for (unsigned int j = 0; j < A.cols(); ++j) {
      // We use eigen block operations to set the values of
      // each $n \times n$ block.
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) =
          A(i, j) * B; // $\in \mathbb{R}^{(n \times n)}$
    }
  }
}

MatrixXd sparseSPD(int n) {
    double rate = 1.0 / n;
    // R is a random sparse matrix.
    ArrayXXd R = ArrayXXd::Random(n,n);    
    R = R * (R.abs() < rate).cast<double>();
    // B is diagonally dominant.
    MatrixXd B = R.matrix();
    B.diagonal() = R.abs().matrix().rowwise().sum();
    // A is symmetric and strictly diagonally dominant.
    MatrixXd A = B + B.transpose() + MatrixXd::Identity(n,n);
    return A;
}

int main() {
  int n = 8;
  std::srand(31);
  MatrixXd A = sparseSPD(n);
  SparseMatrix<double> As = A.sparseView();
  As.makeCompressed();
  
  ArrayXd diagA;
  SparseMatrix<double> Xdiag;
  MatrixXd C, X;
  
  int ans=0;
  std::cin >> ans;
  std::cout << std::setprecision(2);
  switch(ans){
    case 0: // Testing solveDiagSylvesterEq()
            diagA = ArrayXd::LinSpaced(n,1,n);    
            Xdiag = solveDiagSylvesterEq( diagA );
            std::cout << MatrixXd(Xdiag) << std::endl;
            break;
    case 1: // Testing sparseKron()
            kron(A,A,C);
            std::cout << (C - sparseKron(As)).norm() << std::endl;
            break;
    case 2: // Testing solveSpecialSylvesterEq()
            X = solveSpecialSylvesterEq( As );
            std::cout << X << std::endl;
            break;
  }
  return 0;
}
