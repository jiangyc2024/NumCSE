
#include "Sylvester.hpp"

using namespace Eigen;

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

int main() {
    // Initialization
    unsigned int n = 6;
    MatrixXd A(n,n);
    A <<  4, -1,  0,  0,  0,  0,
         -1,  4,  2,  0,  2,  0,
          0,  2,  4,  0,  0, -1,
          0,  0,  0,  4, -1,  0,
          0,  2,  0, -1,  4, -1,
		  0,  0, -1,  0, -1,  4;
    
	// Caste $A$ to 'SparseMatrix' and compress it to store the matrix in CCS format
	SparseMatrix<double> As = A.sparseView();
	As.makeCompressed();
    
    std::cout<< A << std::endl;
    
    
    std::vector<double> stdDiagA(n,3.21);
    VectorXd eigDiagA = VectorXd::Constant(n,1.23);
    //~ std::cout << "solveDiagSylvesterEq(std):\n";
    //~ std::cout<< solveDiagSylvesterEq( stdDiagA ) <<std::endl;
    //~ std::cout << "solveDiagSylvesterEq(eig):\n";
    //~ std::cout<< solveDiagSylvesterEq( eigDiagA ) <<std::endl;
    
    SparseMatrix<double> Bs = sparseKron( As );
    std::cout << "sparseKron():\n";
    //~ std::cout << Bs << std::endl;
    
    MatrixXd C;
    kron(A,A,C);
    
    double err = (C.sparseView() - Bs).norm();
    std::cout << err << std::endl;    
    
    
    std::cout << "solveSpecialSylvesterEq():\n";
    MatrixXd X = solveSpecialSylvesterEq( As );
    std::cout<< X <<std::endl;
    std::cout << "\nError:\n";
    std::cout<< (X + A*X*A - A).norm() <<std::endl;
    
    return 0;
}
