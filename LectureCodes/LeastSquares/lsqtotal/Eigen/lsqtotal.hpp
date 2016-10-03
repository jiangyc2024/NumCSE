#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
using Eigen::VectorXd;
using Eigen::MatrixXd;

// computes only solution \Blue{$\Vx$} of fitted consistent LSE
VectorXd lsqtotal(const MatrixXd& A, const VectorXd& b) {
  const unsigned m = A.rows(), n = A.cols(); // No. of rows and columns

  MatrixXd Ab(m, n + 1);
  Ab << A, b; // Ab = [A,b]
  MatrixXd V = Ab.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).matrixV(); // see \eqref{tlsq:1}

  double s = V(n, n);
  if (s == 0) {
    std::cerr << "No solution!\n"; exit(1);
  }

  return -V.col(n).head(n) / s; // see \eqref{tlsq:3};
}
