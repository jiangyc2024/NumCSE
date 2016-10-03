#include <Eigen/Dense>
#include <Eigen/SVD>
using Eigen::VectorXd;
using Eigen::MatrixXd;

VectorXd lsqsvd(const MatrixXd& A, const VectorXd& b) {
  Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // fast solution:
  // return svd.solve(b);
  
  // 'manual' solution:
  VectorXd sv = svd.singularValues();
  unsigned r = svd.nonzeroSingularValues(); // no. of nonzero singular values
  MatrixXd U = svd.matrixU(), 
           V = svd.matrixV();

  return V.leftCols(r)*( sv.head(r).cwiseInverse().asDiagonal() * (U.leftCols(r).adjoint()*b) );
}
