// Solves constrained linear least squares problem \eqref{clsq} with \texttt{dim} passing \Blue{$d$}
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

using Eigen::MatrixXd;
using Eigen::VectorXd;

void clsq(const MatrixXd& A, const unsigned dim, double& c, VectorXd& n) {
  unsigned p = A.cols(), m = A.rows();

  if (p < dim + 1) {
    std::cerr << "not enough unknowns\n"; return;
  }
  if (m < dim)  {
    std::cerr << "not enough equations\n"; return;
  }

  m = std::min(m, p);
  // First step: orthogonal transformation, see Code~\ref{mc:qrlsqsolve}
  MatrixXd R = A.householderQr().matrixQR().template triangularView<Eigen::Upper>();
  // compute matrix V from SVD composition of R, solve \eqref{eq:HRmin2}
  MatrixXd V = R.block(p - dim, p - dim, m + dim - p, dim)
                .jacobiSvd(Eigen::ComputeFullV).matrixV();
  n = V.col(dim - 1);

  MatrixXd R_topleft = R.topLeftCorner(p - dim, p - dim);
  c = -(R_topleft.template triangularView<Eigen::Upper>()
        .solve(R.block(0, p - dim, p - dim, dim)) * n)(0);
}
