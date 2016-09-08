# include <Eigen/Dense>
# include <Eigen/SVD>
# include <Eigen/QR>

using Eigen::VectorXd;
using Eigen::MatrixXd;

void qr_solve(const MatrixXd& A, const VectorXd& b, VectorXd& x) {
  Eigen::HouseholderQR<MatrixXd> solver(A);
  x = solver.solve(b);
  // or equivalent:
  // x = A.householderQr().solve(b);
}

void lu_solve(const MatrixXd& A, const VectorXd& b, VectorXd& x) {
  x = A.lu().solve(b); // 'lu()' is short for 'partialPivLu()'
}

void svd_solve(const MatrixXd& A, const VectorXd& b, VectorXd& x) {
  x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}
