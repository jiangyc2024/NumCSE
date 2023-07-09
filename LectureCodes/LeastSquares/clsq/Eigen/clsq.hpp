// Solves constrained linear least squares problem \eqref{clsq} with
// \texttt{dim} passing \Blue{$d$}
#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

namespace clsq {


using std::runtime_error;
using std::numeric_limits;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// clang-format off
inline
/* SAM_LISTING_BEGIN_0 */
// Solves constrained linear least squares problem
// \eqref{clsq} with \texttt{dim} passing \Blue{$d$}
std::pair<Eigen::VectorXd, Eigen::VectorXd> clsq(const MatrixXd& A,
						 const unsigned dim) {
  const unsigned p = A.cols();
  unsigned m = A.rows();
  if (p < dim + 1) {
    throw runtime_error("not enough unknowns");
  }
  if (m < dim) {
   throw runtime_error("not enough equations"); 
  }
  m = std::min(m, p); // Number of variables 
  // First step: orthogonal transformation, see Code~\ref{mc:qrlsqsolve}
  MatrixXd R = A.householderQr().matrixQR().template triangularView<Eigen::Upper>();
  // compute matrix V from SVD composition of R, solve \eqref{eq:HRmin2}
  MatrixXd V = R.block(p - dim, p - dim, m + dim - p, dim)
                .jacobiSvd(Eigen::ComputeFullV).matrixV();
  const VectorXd n = V.col(dim - 1); // Norm-constrained part of solution vector
  // Compute free part of solution vector 
  const auto R_topleft = R.topLeftCorner(p - dim, p - dim);
  // Check for singular matrix
  const auto R_diag = R_topleft.diagonal().cwiseAbs();
  if (R_diag.minCoeff() < (numeric_limits<double>::epsilon())*R_diag.maxCoeff()) {
    throw runtime_error("Upper left block of R not regular");
  }
  const VectorXd c = -(R_topleft.template triangularView<Eigen::Upper>()).
                  solve(R.block(0, p - dim, p - dim, dim) * n);
  return {c,n}; 
}
/* SAM_LISTING_END_0 */
// clang-format on

inline
/* SAM_LISTING_BEGIN_9 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> clsq2(const Eigen::MatrixXd &A,
                                                  const unsigned d) {
  const unsigned int n = A.cols();
  const unsigned int m = A.rows();
  assert(m >= n);
  const Eigen::MatrixXd R =
    A.householderQr().matrixQR().template triangularView<Eigen::Upper>(); // \Label[line]{fh:1}
  const auto V = R.block(d, d, m -d , n-d) // \Label[line]{fh:2}
                     .jacobiSvd(Eigen::ComputeFullV) // \Label[line]{fh:3}
                     .matrixV(); // \Label[line]{fh:4}
  const auto y = V.col(n-d - 1); // \Label[line]{fh:5}
  const auto b = R.block(0, d, d, n-d) * y; // \Label[line]{fh:6}
  const auto S = R.topLeftCorner(d, d); // \Label[line]{fh:7}
  const Eigen::VectorXd D = S.diagonal().cwiseAbs(); // \Label[line]{fh:8}
  if (D.minCoeff() < (numeric_limits<double>::epsilon()) * D.maxCoeff()) { // \Label[line]{fh:9}
    throw runtime_error("Upper left block of R not regular"); // \Label[line]{fh:A}
  }
  const auto z = -(S.template triangularView<Eigen::Upper>()).solve(b); // \Label[line]{fh:B}
  return {z, y}; 
}
/* SAM_LISTING_END_9 */


} //namespace clsq