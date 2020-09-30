// Solves constrained linear least squares problem \eqref{clsq} with \texttt{dim} passing \Blue{$d$}
#include <iostream>
#include <cmath>
#include <exception>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Solves constrained linear least squares problem
// \eqref{clsq} with \texttt{dim} passing \Blue{$d$}
std::pair<Eigen::VectorXd, Eigen::VectorXd> clsq(const MatrixXd& A,
						 const unsigned dim) {
  unsigned p = A.cols(), m = A.rows();
  if (p < dim + 1) throw runtime_error("not enough unknowns");
  if (m < dim)  throw runtime_error("not enough equations"); 
  m = std::min(m, p); // Number of variables 
  // First step: orthogonal transformation, see Code~\ref{mc:qrlsqsolve}
  MatrixXd R = A.householderQr().matrixQR().template triangularView<Eigen::Upper>();
  // compute matrix V from SVD composition of R, solve \eqref{eq:HRmin2}
  MatrixXd V = R.block(p - dim, p - dim, m + dim - p, dim)
                .jacobiSvd(Eigen::ComputeFullV).matrixV();
  VectorXd n = V.col(dim - 1); // Norm-constrained part of solution vector
  // Compute free part of solution vector 
  const auto R_topleft = R.topLeftCorner(p - dim, p - dim);
  // Check for singular matrix
  const auto R_diag = R_topleft.diagonal().cwiseAbs();
  if (R_diag.minCoeff() < (numeric_limits<double>::epsilon())*R_diag.maxCoeff())
    throw runtime_error("Upper left block of R not regular");
  VectorXd c = -(R_topleft.template triangularView<Eigen::Upper>()).
                  solve(R.block(0, p - dim, p - dim, dim) * n);
  return {c,n}; 
}
/* SAM_LISTING_END_0 */
